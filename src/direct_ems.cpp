#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectEMS : public Integrator {
public:
    DirectEMS(const PropertyList& props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        Normal3f n = its.shFrame.n;
        MatrixXf final_G = Transform_G(n); //Nomalize to n = (0, 0, 1)

        EmitterQueryRecord lRec(its.p);
        auto light_list = scene->getLights();
        Color3f light_L(0.0);
        //Test if this intersect area is an emit light.
        Color3f Le(0.0f);
        EmitterQueryRecord new_lRec(ray.o, its.p, n); // ref, p, n
        if (its.mesh->isEmitter()) {
            Le = its.mesh->getEmitter()->eval(new_lRec);
        }

        std::vector<Color3f> Lr;
        std::vector<float> light_pdf;
        for (auto Light : light_list) {
            // L: light sample
            light_L = Light->sample(lRec, sampler->next2D()); // area light.
            light_pdf.push_back(Light->pdf(lRec));
            
            //VISIBILITY: try to find the visibility.
            float V = 1;
            float cos_theta_i = lRec.wi.dot(n);
            float cos_theta_0 = lRec.n.dot(-lRec.wi);
            Ray3f new_ray = Ray3f(its.p, lRec.wi);
            new_ray.maxt = (lRec.p - new_ray.o).norm() - Epsilon;
            new_ray.mint = Epsilon;
            //Test if there is any occlusion.
            if (cos_theta_i <= 0 || cos_theta_0 <= 0 || scene->rayIntersect(new_ray)) {
                V = 0;
                Lr.push_back(Color3f(0.0f));
                continue;
            }
            float G = V * cos_theta_i;

            //BRDF: Nomalize to n = (0, 0, 1)
            BSDFQueryRecord bRec(-final_G * ray.d, final_G * lRec.wi, ESolidAngle); //wi(p->camera), wo(p->light), measure
            bRec.uv = its.uv;
            Color3f brdf_eval = its.mesh->getBSDF()->eval(bRec);

            Lr.push_back(G * light_L * brdf_eval);
        }

        Color3f Lr_sum(0.0f);
        float pdf_sum = 0.0;
        for (int i = 0; i < Lr.size(); i++) {
            pdf_sum += light_pdf[i];
            Lr_sum += Lr[i] * light_pdf[i];
        }
        if (pdf_sum == 0.0) {
            return Le;
        }
        Lr_sum /= pdf_sum;
        return Le + Lr_sum;
    }

    std::string toString() const {
        return "DirectEMS";
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
};

NORI_REGISTER_CLASS(DirectEMS, "direct_ems");
NORI_NAMESPACE_END