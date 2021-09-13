#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class SpotEmitter : public Emitter {
public:
    SpotEmitter(const PropertyList& props) {
        m_radiance = props.getColor("radiance");
        m_cosThetaMax = props.getFloat("cosThetaMax", 0.5);
        m_radius = props.getFloat("radius");
        m_center = props.getPoint3("center");
        m_n = props.getVector3("normal");
    }

    virtual std::string toString() const override {
        return tfm::format(
            "SpotLight[\n"
            "  radiance = %s,\n"
            "  cosTheta = %f,\n"
            "]",
            m_radiance.toString(), m_cosThetaMax);
    }

    // A function can be called by the integrators whenever a ray intersects your mesh emitter
    virtual Color3f eval(const EmitterQueryRecord& lRec) const override {
        if (!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        if ((lRec.p - lRec.ref).dot(lRec.n) < 0) { // Front
            return m_radiance;
        }
        else { // Back
            return Color3f(0.0f);
        }
    }

    virtual Color3f sample(EmitterQueryRecord& lRec, const Point2f& sample) const override {
        // Do not need a shape, just use the analytical way.
        Vector3f w = Warp::squareToUniformSphereCap(sample, m_cosThetaMax);

        lRec.ref = lRec.ref;
        lRec.n = Transform_G(m_n).transpose() * w;
        lRec.p = m_radius * lRec.n + m_center;
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.pdf = Warp::squareToUniformSphereCapPdf(w, m_cosThetaMax);
        double cos_theta_0 = -lRec.wi.dot(lRec.n);
        double dis2 = (lRec.p - lRec.ref).squaredNorm();
        if (lRec.pdf == 0 || cos_theta_0 <= 0.0f) {
            return Color3f(0.0f);
        }
        return eval(lRec) / lRec.pdf / dis2 * cos_theta_0;
    }

    virtual float pdf(const EmitterQueryRecord& lRec) const override {
        // Do not need a shape, just use the analytical way.
        double cos_theta_0 = -lRec.wi.dot(lRec.n);
        if (cos_theta_0 >= 0) { // Front
            return Warp::squareToUniformSphereCapPdf(Transform_G(m_n) * lRec.n, m_cosThetaMax);
        }
        else { // Back
            return 0.0f;
        }
    }

    virtual Color3f samplePhoton(Ray3f& ray, const Point2f& sample1, const Point2f& sample2) const override {
        // Do not need a shape, just use the analytical way.
        Normal3f n_local = Warp::squareToUniformSphereCap(sample1, m_cosThetaMax);
        Normal3f n = Transform_G(m_n) * n_local;
        Vector3f w_local = Warp::squareToCosineHemisphere(sample2);
        Vector3f w = Transform_G(n) * w_local;

        ray.d = n; ray.o = m_center + m_radius * n;
        return M_PI * m_radiance / Warp::squareToUniformSphereCapPdf(n_local, m_cosThetaMax);

        ray.d = w; ray.o = m_center + m_radius * n;
        float pdf = Warp::squareToUniformSphereCapPdf(n_local, m_cosThetaMax);
        if (pdf > 0.0f) {
            return M_PI * m_radiance / pdf;
        }
        else {
            return Color3f(0.0f);
        }
    }

    MatrixXf Transform_G(Vector3f n_) const {
        Vector3f n = n_;
        //Nomalize to n = (0, 0, 1)
        MatrixXf G1(3, 3);
        if (n.z() == 0) {
            if (n.y() != 0) {
                G1 << 1.0, 0.0, 0.0,
                    0.0, 0.0, -1.0,
                    0.0, n.y() / std::abs(n.y()), 0.0;
            }
            else if (n.x() != 0) {
                G1 << 0.0, 0.0, -1.0,
                    0.0, 1.0, 0.0,
                    n.x() / std::abs(n.x()), 0.0, 0.0;
            }
        }
        else {
            G1 << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, n.z() / std::abs(n.z());
        }
        n = G1 * n;
        double gamma = n.z() / std::sqrt(n.z() * n.z() + n.y() * n.y()), sigma = n.y() / std::sqrt(n.z() * n.z() + n.y() * n.y());
        MatrixXf G2(3, 3);
        G2 << 1.0, 0.0, 0.0,
            0.0, gamma, -sigma,
            0.0, sigma, gamma;
        n = G2 * n;

        gamma = n.z() / std::sqrt(n.z() * n.z() + n.x() * n.x());
        sigma = n.x() / std::sqrt(n.z() * n.z() + n.x() * n.x());
        MatrixXf G3(3, 3);
        G3 << gamma, 0.0, -sigma,
            0.0, 1.0, 0.0,
            sigma, 0.0, gamma;
        n = G3 * n;
        MatrixXf final_G = G3 * G2 * G1;
        return final_G;
    }

protected:
    Color3f m_radiance;
    Point3f m_center;
    float m_radius;
    float m_cosThetaMax;
    Normal3f m_n;
};

NORI_REGISTER_CLASS(SpotEmitter, "spot_emitter")
NORI_NAMESPACE_END