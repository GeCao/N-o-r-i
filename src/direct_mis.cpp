#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectMIS : public Integrator {
public:
    DirectMIS(const PropertyList& props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Point2f sample_p = sampler->next2D();
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

        Color3f Lr_ems(0.0f);
        float pdf_ems_wie = 0.0;
        float pdf_mat_wie = 0.0;
        for (auto Light : light_list) {
            // L: light sample
            light_L = Light->sample(lRec, sample_p); // area light.
            pdf_ems_wie += Light->pdf(lRec);

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
                continue;
            }
            float G = V * cos_theta_i;

            //BRDF: Nomalize to n = (0, 0, 1)
            BSDFQueryRecord bRec(-final_G * ray.d, final_G * lRec.wi, ESolidAngle); //wi(p->camera), wo(p->light), measure
            bRec.uv = its.uv;
            Color3f brdf_eval = its.mesh->getBSDF()->eval(bRec);
            pdf_mat_wie += its.mesh->getBSDF()->pdf(bRec);

            Lr_ems += G * light_L * brdf_eval;
        }
        pdf_mat_wie /= light_list.size();
        pdf_ems_wie /= light_list.size();
        Lr_ems /= light_list.size();

        //BRDF: Nomalize to n = (0, 0, 1)
        //    wi(p -> light), wo(p -> camera), measure
        BSDFQueryRecord bRec(-final_G * ray.d, Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
        bRec.uv = its.uv;
        Color3f fr_with_pdf = its.mesh->getBSDF()->sample(bRec, sample_p);
        float pdf_mat_wim = its.mesh->getBSDF()->pdf(bRec);
        float pdf_ems_wim = 0.0;

        lRec.ref = its.p;
        lRec.wi = final_G.transpose() * bRec.wo;
        light_L = Color3f(0.0);
        float V = 1;
        Intersection emit_its;
        Ray3f new_ray = Ray3f(its.p, lRec.wi);
        new_ray.mint = Epsilon;
        //->Test if there is any intersection.
        if (scene->rayIntersect(new_ray, emit_its)) {
            if (emit_its.mesh->isEmitter()) {
                V = 1;
                lRec.p = emit_its.p;
                lRec.n = emit_its.shFrame.n;
                light_L = emit_its.mesh->getEmitter()->eval(lRec);
                pdf_ems_wim = emit_its.mesh->getEmitter()->pdf(lRec);
            }
            else {
                //Intersect but it is not light
                V = 0;
                light_L = Color3f(0.0f);
            }
        }
        else {
            //No intersect
            V = 0;
            light_L = Color3f(0.0f);
        }
        float cos_theta_i = lRec.wi.dot(n);
        float cos_theta_0 = lRec.n.dot(-lRec.wi);
        if (cos_theta_i <= 0 || cos_theta_0 <= 0) {
            V = 0;
            light_L = Color3f(0.0f);
        }

        Color3f Lr_mat =  fr_with_pdf * light_L * V;

        Color3f res_Lr = Le;

        if (pdf_ems_wie + pdf_mat_wie > 0) {
            res_Lr += Lr_ems * pdf_ems_wie / (pdf_ems_wie + pdf_mat_wie);
        }
        if (pdf_ems_wim + pdf_mat_wim > 0) {
            res_Lr += Lr_mat * pdf_mat_wim / (pdf_ems_wim + pdf_mat_wim);
        }

        return res_Lr;
    }

    std::string toString() const {
        return "DirectMIS";
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

NORI_REGISTER_CLASS(DirectMIS, "direct_mis");
NORI_NAMESPACE_END