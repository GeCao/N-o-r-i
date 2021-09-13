#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bitmap.h>
#include <nori/bbox.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class EnvEmitter : public Emitter {
public:

    EnvEmitter(const PropertyList& props) {
        m_filename = props.getString("filename");
        m_scaler = props.getFloat("scale");
        Bitmap bitmap(m_filename);
        m_image = bitmap;
        m_isEnv = true;

        m_numRows = m_image.rows();
        m_numCols = m_image.cols();

        m_theta_pdf.resize(m_numRows, 1);
        m_theta_cdf.resize(m_numRows, 1);

        /// Gray version
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> gray;
        gray.resize(m_numRows, m_numCols);
        for (int y = 0; y < m_numRows; ++y)
            for (int x = 0; x < m_numCols; ++x)
                gray.coeffRef(y, x) = (m_image(y, x).r() + m_image(y, x).g() + m_image(y, x).b()) / 3;

        /// Multipled by sinTheta
        for (int y = 0; y < m_numRows; y++)
            gray.row(y) = gray.row(y) * sin(y * M_PI / (m_numRows - 1));
        m_phi_pdf = gray;

        
        /// Create the pdf of theta and pdf of phi
        for (int y = 0; y < m_numRows; y++) {
            m_theta_pdf(y) = gray.row(y).sum();
            if(gray.row(y).sum() != 0)
                m_phi_pdf.row(y) = gray.row(y) / gray.row(y).sum();
        }
        m_theta_pdf = m_theta_pdf / m_theta_pdf.sum();

        /// Create the cdf of theta and cdf of phi
        m_theta_cdf = m_theta_pdf;
        for (int y = 1; y < m_numRows; y++)
            m_theta_cdf(y) += m_theta_cdf(y - 1);
        m_phi_cdf = m_phi_pdf;
        for (int y = 0; y < m_numRows; y++)
            for (int x = 1; x < m_numCols; x++)
                m_phi_cdf(y, x) += m_phi_cdf(y, x - 1);
        
    }

    /// Get the big sphere from the scene
    void activate(const Scene* scene) {
        m_center = scene->getBVH()->getBoundingBox().getCenter();
        m_radius = scene->getBoundingBox().getExtents().norm() / 2 * m_scaler;
    }

    /// Make sure the generated index will not voilate the boundary
    int clamp(const int& x, const int& min, const int& max) const {

        int index = x;
        if (x < min)
            index = min;
        else if (x > max)
            index = max;
        return index;

    }

    /// Find the nearest index in CDF
    int getIndex(const Eigen::VectorXf& cdf, const float& sample) const {

        int index = 0;
        for (int i = 0; i < cdf.size(); i++) {
            if (sample <= cdf(i)) {
                float distance1 = cdf(i) - sample,
                      distance2 = sample - cdf(i - 1);
                if (distance1 > distance2)
                    index = i - 1;
                else
                    index =  i;
                break;
            }
        }
        return index;

    }

    Color3f getEnvRadiance(const Ray3f& ray) const {

        Vector3f d = (m_center - ray.o).normalized();
        float cosTheta = ray.d.dot(d);
        float d1 = (m_center - ray.o).norm(),
              sinTheta = sqrt(1 - cosTheta * cosTheta);

        float d2 = d1 * cosTheta,
              d3 = d1 * sinTheta;
        float d4 = sqrt(m_radius * m_radius - d3 * d3);

        EmitterQueryRecord lRec;
        lRec.p = ray.o + (d2 + d4) *ray.d;

        return eval(lRec);

    }

    virtual std::string toString() const override {
        return tfm::format(
            "EnvironmentLight[\n"
            "  filename = %s,\n"
            "]",
            m_filename);
    }

    void getLocation(const EmitterQueryRecord& lRec, int& row_index, int& col_index) const {

        Vector3f direction = lRec.p - m_center;
        if (direction.y() > m_radius)
            direction.y() = m_radius;
        float Theta = acos(direction.y() / m_radius);
        float Phi = atan2(direction.x(), direction.z());
        if (Phi < 0)
            Phi += 2 * M_PI;
        
        row_index = clamp(floor(Theta * INV_PI * m_numRows), 0, m_numRows - 1);
        col_index = clamp(floor(Phi * 0.5f * INV_PI * m_numCols), 0, m_numCols - 1);

    }

    virtual Color3f eval(const EmitterQueryRecord& lRec) const override {

        int row_index, col_index;
        getLocation(lRec, row_index, col_index);

        return m_image(row_index, col_index);
    }

    virtual Color3f sample(EmitterQueryRecord& lRec, const Point2f& sample) const override {

        int rowIndex, colIndex;
        rowIndex = getIndex(m_theta_cdf, sample.x());
        Eigen::VectorXf selected_row = m_phi_cdf.row(rowIndex);
        colIndex = getIndex(selected_row, sample.y());
        rowIndex = clamp(rowIndex, 0, m_numRows - 1);
        colIndex = clamp(colIndex, 0, m_numCols - 1);
        
        lRec.pdf = m_theta_pdf(rowIndex) * m_phi_pdf(rowIndex, colIndex);

        float theta = rowIndex / (m_numRows - 1) * M_PI,
              phi = colIndex / (m_numCols - 1) * 2 * M_PI;

        lRec.p = Point3f(m_radius * sin(theta) * cos(phi), m_radius * sin(theta) * sin(phi), m_radius * cos(theta));
        lRec.n = (m_center - lRec.p).normalized();
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (lRec.p - lRec.ref).norm() - Epsilon);

        if (eval(lRec).minCoeff() <= 0 || pdf(lRec) <= 0)
            return Color3f(0.0f);
        else
            return eval(lRec) / pdf(lRec);
    }

    virtual float pdf(const EmitterQueryRecord& lRec) const override {

        float cosTheta0 = lRec.n.dot(-lRec.wi),
              sinTheta0 = sqrt(1 - cosTheta0 * cosTheta0);
        if (cosTheta0 <= 0 || sinTheta0 <= 0)
            return 0.0f;

        return m_numCols * m_numRows / (2 * M_PI * M_PI * sinTheta0) * lRec.pdf;

    }

private:
    std::string m_filename;
    Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_image;
    Eigen::VectorXf m_theta_pdf;
    Eigen::VectorXf m_theta_cdf;
    Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_phi_pdf;
    Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_phi_cdf;
    Point3f m_center;
    float m_radius;
    int m_scaler;
    int m_numRows;
    int m_numCols;
};

NORI_REGISTER_CLASS(EnvEmitter, "environment")
NORI_NAMESPACE_END