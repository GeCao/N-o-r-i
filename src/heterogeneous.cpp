#include <nori/medium.h>
#include <nori/scene.h>
#include <fstream>

NORI_NAMESPACE_BEGIN

class Heterogeneous : public Medium {
public:
    Heterogeneous(const PropertyList& props) {
        std::string filename = props.getString("filepath", "D:\\Homework\\computer_graphics\\nori\\scenes\\Validation\\fluid\\smoke.vol");
        m_scale = props.getVector3("scale", Vector3f(1.0f, 1.0f, 1.0f));
        m_translate = props.getVector3("translate", Vector3f(0.0f, 0.0f, 0.0f));
        m_particle_density = props.getFloat("particledensity", 4.5e-5);
        m_density_scale = props.getFloat("densityscale", 1.0f);
        m_integ_res = props.getInteger("resolution", 256);
        Volumetric_g = props.getFloat("volumetric_g", 0.0f);
        //Read the file:
        std::ifstream is(filename, std::ios::binary);

        char* buf = (char*)malloc(4 * sizeof(char));
        is.get(buf, 3 + 1); // VOL;
        char version;
        is.read((char*)&version, sizeof(char)); //char version = 3;

        int value, xres, yres, zres;
        is.read((char*)&value, sizeof(int)); // 1;
        is.read((char*)&xres, sizeof(int));
        is.read((char*)&yres, sizeof(int));
        is.read((char*)&zres, sizeof(int));
        m_res = Vector3i(xres, yres, zres);

        float minX, minY, minZ, maxX, maxY, maxZ;
        is.read((char*)&value, sizeof(int)); // 1;
        is.read((char*)&minX, sizeof(float));
        is.read((char*)&minY, sizeof(float));
        is.read((char*)&minZ, sizeof(float));
        is.read((char*)&maxX, sizeof(float));
        is.read((char*)&maxY, sizeof(float));
        is.read((char*)&maxZ, sizeof(float));
        m_minCorner = Point3f(minX, minY, minZ).cwiseProduct(m_scale);
        m_maxCorner = Point3f(maxX, maxY, maxZ).cwiseProduct(m_scale);
        m_minCorner += m_translate;
        m_maxCorner += m_translate;
        m_girdsize.x() = (m_maxCorner.x() - m_minCorner.x()) / (xres - 1.0f);
        m_girdsize.y() = (m_maxCorner.y() - m_minCorner.y()) / (yres - 1.0f);
        m_girdsize.z() = (m_maxCorner.z() - m_minCorner.z()) / (zres - 1.0f);

        m_density.clear();
        m_maxdensity = 0.0;
        while (true) {
            float density;
            is.read((char*)&density, sizeof(float));
            density *= m_density_scale;
            m_density.push_back(density);
            if (is.eof()) { break; }
            if (density > m_maxdensity) { m_maxdensity = density; }
        }
        m_densitysize = m_density.size() - 1; // Do not know the reason, but here should hava a reduction to make it work.
        cout << "MaxDensity: " << m_maxdensity << std::endl;

        is.close();
        /*
        * A visualization code.
        FILE* fp;
        if (!(fp = fopen("data_density.plt", "w"))) {
            printf("Ð´ÈëÎÄ¼þÊ§°Ü!\n");
            exit(-1);
        }
        fprintf(fp, "TITLE =\"smoke_vol\"\n");
        fprintf(fp, "variables = \"x \" \"y \" \"z\" \"rho\" \n");
        fprintf(fp, "zone i=%d, j=%d, k=%d \n", xres, yres, zres);
        fprintf(fp, "f=point\n");
        int idx = 0;
        for (int k = 0; k < zres; k++) {
            for (int j = 0; j < yres; j++) {
                for (int i = 0; i < xres; i++) {
                    fprintf(fp, "%d %d %d %lf\n", i, j, k, m_density[idx++] / m_maxdensity);
                }
            }
        }
        fclose(fp);
        */
        int idx = 0;
        m_sigma_a = std::vector<float>(m_densitysize, 0.0f);
        m_sigma_s = std::vector<float>(m_densitysize, 0.0f);
        for (int k = 0; k < zres; k++) {
            for (int j = 0; j < yres; j++) {
                for (int i = 0; i < xres; i++) {
                    //lambda = 525 nm, particle diameter = 10 um, eta = 1.33, density = density
                    m_sigma_s[idx] = ComputeSigma_s(5.25e-7, m_particle_density, 1.33, m_density[idx]);

                    //Need to implement by CaoGe: How to get sigma_a by density?
                    //m_sigma_a[idx] = (m_density[idx]) / (maxX - minX);
                    m_sigma_a[idx] = 0.07 * m_sigma_s[idx];
                    idx++;
                }
            }
        }

        std::cout << "Suceessfully loaded hetero-medium" << std::endl;
    }

    int convert_index(int i, int j, int k) const {
        //return i * m_res.y() * m_res.z() + j * m_res.z() + k;
        return k * m_res.y() * m_res.x() + j * m_res.x() + i;
    }

    virtual Vector3f squareToHenyen_Greenstein(const Point2f sample) const {
        float g = Volumetric_g;
        float cos_theta = (1 + g * g - pow((1 - g * g) / (1 - g + 2 * g * sample.x()), 2)) / (2.0 * g);
        float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        float fai = 2 * M_PI * sample.y();
        return Vector3f(sin_theta * cos(fai), sin_theta * sin(fai), cos_theta);
    }

    virtual float squareToHenyen_GreensteinPdf(const Vector3f m) const {
        return ComputePhaseFunction(m.z());
    }

    virtual float ComputePhaseFunction(float costheta) const {
        /*Henyen_Green*/
        float g = Volumetric_g;
        return 1.0 / 4.0 / M_PI * (1 - Volumetric_g * Volumetric_g) / pow(1 + g * g - 2 * g * costheta, 1.5);
        /*Rayleigh*/
        //return 3.0 / 16.0 / M_PI * (1.0 + costheta * costheta);
    }

    virtual bool GetParamsOnMediumPoint(MediumSamplingRecord& MRec) const {
        Point3f xt = MRec.p;
        if (xt.x() < m_minCorner.x() || xt.y() < m_minCorner.y() || xt.z() < m_minCorner.z()) { return false; }
        if (xt.x() > m_maxCorner.x() || xt.y() > m_maxCorner.y() || xt.z() > m_maxCorner.z()) { return false; }
        float x_interpolation = (xt.x() - m_minCorner.x()) / m_girdsize.x();
        float y_interpolation = (xt.y() - m_minCorner.y()) / m_girdsize.y();
        float z_interpolation = (xt.z() - m_minCorner.z()) / m_girdsize.z();
        int i = (int)(x_interpolation), j = (int)(y_interpolation), k = (int)(z_interpolation);
        MRec.sigmaA = 0.0f;  MRec.sigmaS = 0.0f;
        
        for (int left = i; left < 1; left++) {
            for (int forward = j; forward < 1; forward++) {
                for (int down = k; down < 1; down++) {
                    int idx = convert_index(left, forward, down);
                    float weight = (m_girdsize.x() - abs(xt.x() - (left * m_girdsize.x() + m_minCorner.x()))) * \
                        (m_girdsize.y() - abs(xt.y() - forward * m_girdsize.y() + m_minCorner.y())) * \
                        (m_girdsize.z() - abs(xt.z() - down * m_girdsize.z() + m_minCorner.z()));
                    if (idx >= m_res.x() * m_res.y() * m_res.z() || idx < 0) { continue; }
                    MRec.sigmaA += weight * m_sigma_a[idx];
                    MRec.sigmaS += weight * m_sigma_s[idx];
                }
            }
        }
        
        MRec.sigmaA = m_sigma_a[convert_index(i, j, k)];
        MRec.sigmaS = m_sigma_s[convert_index(i, j, k)];
        return true;
    }

    virtual float evalTransmittance(Ray3f& ray) const {
        /*Ray Marching*/
        float tmin = ray.mint, tmax = ray.maxt, Tr = 0.0f;
        if (!get_tmin_tmax_by_AABBIntersect(ray)) { // No medium detected
            return 1.0;
        }
        else { //Medium detected
            tmin = std::max(tmin, ray.mint); tmax = std::min(tmax, ray.maxt);
            if (tmax <= tmin) { return 1.0; }
            int res = m_integ_res;
            float dz = (tmax - tmin) / (res - 1);
            MediumSamplingRecord MRec;
            MRec.p = ray.o + tmin * ray.d;
            for (int k = 1; k < res; k++) {
                MRec.p = MRec.p + (k*dz)* ray.d;
                GetParamsOnMediumPoint(MRec);
                Tr += (dz) * (MRec.sigmaA + MRec.sigmaS);
            }
            Tr = std::exp(-Tr);
            return Tr;
        }
    }

    virtual float freeFlightDistance(Ray3f ray_, float tmax_, float sample) const {
        /*Ray Marching*/
        Ray3f ray = Ray3f(ray_.o, ray_.d);
        ray.mint = ray_.mint; ray.maxt = tmax_;
        float tmin = ray.mint, tmax = tmax_;
        if (!get_tmin_tmax_by_AABBIntersect(ray)) { // No medium detected
            return tmax_;
        }
        else { //Medium detected
            tmin = std::max(tmin, ray.mint); tmax = std::min(tmax, ray.maxt);
            if (tmax <= tmin) { return tmax_; } // Medium detected but not in the routine.
        }

        int res = m_integ_res;
        float dz = (tmax - tmin) / (res - 1);
        float Tr = 1.0f;
        float Tr_withnoexp = 0.0f;
        MediumSamplingRecord MRec;
        MRec.p = ray.o + tmin * ray.d;
        for (int k = 1; k < res; k++) {
            MRec.p = MRec.p + (k * dz) * ray.d;
            GetParamsOnMediumPoint(MRec);
            Tr_withnoexp += -dz * (MRec.sigmaA + MRec.sigmaS);
            Tr = exp(Tr_withnoexp);
            if (Tr < sample) {
                //cout << "Tr: " << Tr;
                //cout << ", Tr component: " << exp(-dz * (MRec.sigmaA + MRec.sigmaS)) << std::endl;
                return tmin + k * dz;
            }
        }
        return tmax_;
    }

    std::string toString() const {
        return "heterogeneous medium";
    }

    virtual float ComputeSigma_s(float lambda, float diameter, float eta, float rho) const {
        /*Rayleigh Scatter*/
        if (rho < 1e-5) { return 0.0; }
        float sigma_s = rho * 2 * pow(M_PI, 5) * pow(diameter, 6) / 3 / pow(lambda, 4) * pow((eta * eta - 1) / (eta * eta + 2), 2);
        return sigma_s;
    }

    virtual bool get_tmin_tmax_by_AABBIntersect(Ray3f& ray) const {
        //1st: Test that if the ray.o is inside the AABB.
        bool isinteresect = false;

        //2nd: iterate the 6 faces:
        float dist, d_coeff;
        Point3f its_p;
        std::vector<Point3f> its_p_vec;
        ray.maxt = Epsilon; ray.mint = std::numeric_limits<float>::infinity();

        //Face1: x = Xmin
        dist = m_minCorner.x() - ray.o.x();
        if (ray.d.x() != 0.0f) {
            d_coeff = dist / ray.d.x();
            if (d_coeff >= 0) {
                its_p = ray.o + d_coeff * ray.d;
                its_p_vec.push_back(its_p);
                if (PointOnAABB(its_p) || PointAtAABBLine(its_p)) { ray.mint = std::min(ray.mint, d_coeff); ray.maxt = std::max(ray.maxt, d_coeff); isinteresect = true; }
            }
        }

        //Face2: y = Ymin
        dist = m_minCorner.y() - ray.o.y();
        if (ray.d.y() != 0.0f) {
            d_coeff = dist / ray.d.y();
            if (d_coeff >= 0) {
                its_p = ray.o + d_coeff * ray.d;
                its_p_vec.push_back(its_p);
                if (PointOnAABB(its_p) || PointAtAABBLine(its_p)) { ray.mint = std::min(ray.mint, d_coeff); ray.maxt = std::max(ray.maxt, d_coeff); isinteresect = true; }
            }
        }

        //Face3: z = Zmin
        dist = m_minCorner.z() - ray.o.z();
        if (ray.d.z() != 0.0f) {
            d_coeff = dist / ray.d.z();
            if (d_coeff >= 0) {
                its_p = ray.o + d_coeff * ray.d;
                its_p_vec.push_back(its_p);
                if (PointOnAABB(its_p) || PointAtAABBLine(its_p)) { ray.mint = std::min(ray.mint, d_coeff); ray.maxt = std::max(ray.maxt, d_coeff); isinteresect = true; }
            }
        }

        //Face4: x = Xmax
        dist = m_maxCorner.x() - ray.o.x();
        if (ray.d.x() != 0.0f) {
            d_coeff = dist / ray.d.x();
            if (d_coeff >= 0) {
                its_p = ray.o + d_coeff * ray.d;
                its_p_vec.push_back(its_p);
                if (PointOnAABB(its_p) || PointAtAABBLine(its_p)) { ray.mint = std::min(ray.mint, d_coeff); ray.maxt = std::max(ray.maxt, d_coeff); isinteresect = true; }
            }
        }

        //Face5: y = Ymax
        dist = m_maxCorner.y() - ray.o.y();
        if (ray.d.y() != 0.0f) {
            d_coeff = dist / ray.d.y();
            if (d_coeff >= 0) {
                its_p = ray.o + d_coeff * ray.d;
                its_p_vec.push_back(its_p);
                if (PointOnAABB(its_p) || PointAtAABBLine(its_p)) { ray.mint = std::min(ray.mint, d_coeff); ray.maxt = std::max(ray.maxt, d_coeff); isinteresect = true; }
            }
        }

        //Face6: z = Zmax
        dist = m_maxCorner.z() - ray.o.z();
        if (ray.d.z() != 0.0f) {
            d_coeff = dist / ray.d.z();
            if (d_coeff >= 0) {
                its_p = ray.o + d_coeff * ray.d;
                its_p_vec.push_back(its_p);
                if (PointOnAABB(its_p) || PointAtAABBLine(its_p)) { ray.mint = std::min(ray.mint, d_coeff); ray.maxt = std::max(ray.maxt, d_coeff); isinteresect = true; }
            }
        }

        if (ray.mint == ray.maxt) {
            //The ray.o might be inside the Bounding Box
            if (!PointInsideAABB(ray.o) && !PointOnAABB(ray.o)) {
                //If no, Then the only reason is that this only intersection is on the line or point of bounding box
                its_p = ray.o + ray.maxt * ray.d;
                if (!PointAtAABBLine(its_p)) {
                    //If still no, some wrong with the implement: Check the reason!
                    std::cout << "AABB implement is not correct!" << std::endl;
                    std::cout << "The ray origin point: " << ray.o.x() << ", " << ray.o.y() << ", " << ray.o.z() << std::endl;
                    std::cout << "The ray d: : " << ray.d.x() << ", " << ray.d.y() << ", " << ray.d.z() << std::endl;
                    std::cout << "The ray tmax" << ray.maxt << std::endl;
                    std::cout << "The ray intersect: " << std::endl;
                    for (int i = 0; i < its_p_vec.size(); i++) {
                        std::cout << its_p_vec[i].x() << ", " << its_p_vec[i].y() << ", " << its_p_vec[i].z() << std::endl;
                    }
                }
                else {
                    ray.mint = Epsilon; ray.maxt = std::numeric_limits<float>::infinity(); // The line collision not count.
                    return false;
                }
            }
            else {
                ray.mint = Epsilon; ray.maxt = ray.maxt;
            }
        }

        return isinteresect;
    }

private:
    bool PointInsideAABB(Point3f xt) const {
        if (xt.x() < m_minCorner.x() - 1e-4 || xt.y() < m_minCorner.y() - 1e-4 || xt.z() < m_minCorner.z() - 1e-4) { return false; }
        if (xt.x() > m_maxCorner.x() + 1e-4 || xt.y() > m_maxCorner.y() + 1e-4 || xt.z() > m_maxCorner.z() + 1e-4) { return false; }
        return true;
    }

    bool PointOnAABB(Point3f xt) const {
        if (abs(xt.x() - m_minCorner.x()) < 1e-4 || abs(xt.x() - m_maxCorner.x()) < 1e-4) {
            if (xt.y() >= m_minCorner.y() && xt.y() <= m_maxCorner.y() &&
                xt.z() >= m_minCorner.z() && xt.z() <= m_maxCorner.z()) {
                return true;
            }
        }
        if (abs(xt.y() - m_minCorner.y()) < 1e-4 || abs(xt.y() - m_maxCorner.y()) < 1e-4) {
            if (xt.x() >= m_minCorner.x() && xt.x() <= m_maxCorner.x() &&
                xt.z() >= m_minCorner.z() && xt.z() <= m_maxCorner.z()) {
                return true;
            }
        }
        if (abs(xt.z() - m_minCorner.z()) < 1e-4 || abs(xt.z() - m_maxCorner.z()) < 1e-4) {
            if (xt.y() >= m_minCorner.y() && xt.y() <= m_maxCorner.y() &&
                xt.x() >= m_minCorner.x() && xt.x() <= m_maxCorner.x()) {
                return true;
            }
        }
        return false;
    }

    bool PointAtAABBLine(Point3f xt) const {
        if (abs(xt.x() - m_minCorner.x()) < 1e-4 || abs(xt.x() - m_maxCorner.x()) < 1e-4) {
            if (abs(xt.y() - m_minCorner.y()) < 1e-4 || abs(xt.y() - m_maxCorner.y()) < 1e-4) {
                if (xt.z() >= m_minCorner.z() && xt.z() <= m_maxCorner.z()) {
                    return true;
                }
            }
            else if (abs(xt.z() - m_minCorner.z()) < 1e-4 || abs(xt.z() - m_maxCorner.z()) < 1e-4) {
                if (xt.y() >= m_minCorner.y() && xt.y() <= m_maxCorner.y()) {
                    return true;
                }
            }
        }
        else if (xt.x() >= m_minCorner.x() && xt.x() <= m_maxCorner.x()) {
            if (abs(xt.y() - m_minCorner.y()) < 1e-4 || abs(xt.y() - m_maxCorner.y()) < 1e-4) {
                if (abs(xt.z() - m_minCorner.z()) < 1e-4 || abs(xt.z() - m_maxCorner.z()) < 1e-4) {
                    return true;
                }
            }
        }
        return false;
    }

protected:
    std::vector<float> m_density;
    std::vector<float> m_sigma_s;
    std::vector<float> m_sigma_a;
    float m_maxdensity;
    int m_densitysize;
    Point3f m_minCorner, m_maxCorner, m_girdsize;
    Vector3i m_res;
    Vector3f m_scale, m_translate;

    float m_particle_density;
    float m_density_scale;
    int m_integ_res;
    float Volumetric_g;
};

NORI_REGISTER_CLASS(Heterogeneous, "heterogeneous");
NORI_NAMESPACE_END