#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMIS : public Integrator {
public:
    PathMIS(const PropertyList& props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f Lr(0.0f);
        Color3f mult_t(1.0);
        float pdf_ems_wi_e = 0.0f, pdf_ems_wi_m = 0.0f;
        float pdf_mats_wi_e = 0.0f, pdf_mats_wi_m = 0.0f;
        Ray3f recur_ray = ray;

        bool hit_specular = true; //For the first step, it is specular

        while (true) {
            double weight_mat = 0.0f, weight_ems = 0.0f;
            // Find the surface that is visible in the requested direction
            Intersection its;
            if (!scene->rayIntersect(recur_ray, its))
                break;
            Normal3f n = its.shFrame.n;

            //Le Term:
            EmitterQueryRecord new_lRec(recur_ray.o, its.p, its.shFrame.n); // ref, p, n
            Color3f Le(0.0f);
            pdf_ems_wi_m = 0.0f;
            if (its.mesh->isEmitter()) {
                Le = its.mesh->getEmitter()->eval(new_lRec);
                pdf_ems_wi_m = its.mesh->getEmitter()->pdf(new_lRec);
            }

            //Compute Lr
            if (pdf_ems_wi_m + pdf_mats_wi_m > 0) {
                weight_mat = pdf_mats_wi_m / (pdf_ems_wi_m + pdf_mats_wi_m);
            }
            if (hit_specular) {
                weight_mat = 1.0;
            }
            Lr += weight_mat * mult_t * Le;

            //russianRoulette: Termination
            float succes_probibility = std::min(mult_t.maxCoeff(), 0.99f);
            if (sampler->next1D() <= succes_probibility) {
                //recursion
                mult_t /= succes_probibility;
            }
            else {
                //Termination
                break;
            }

            //BRDF Term:
            BSDFQueryRecord bRec(-its.toLocal(recur_ray.d), Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
            bRec.uv = its.uv;
            Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            pdf_mats_wi_m = its.mesh->getBSDF()->pdf(bRec);

            if (bRec.measure != EDiscrete)
            {
                hit_specular = false;

                //ems part:
                EmitterQueryRecord lRec(its.p);
                auto light_list = scene->getLights();
                int idx_Light = int(sampler->next1D() * light_list.size());
                auto Light = light_list[idx_Light];
                Color3f light_L_with_pdf(0.0);
                Color3f Lr_ems(0.0f);

                // L: light sample
                ShapeQueryRecord ems_Shape;
                lRec.pdf = its.mesh->pdfSurface(ems_Shape);
                light_L_with_pdf = Light->sample(lRec, sampler->next2D()); // area light.

                //VISIBILITY: try to find the visibility.
                float cos_theta_i = lRec.wi.dot(n), cos_theta_0 = lRec.n.dot(-lRec.wi);
                Ray3f new_ray = Ray3f(its.p, lRec.wi);
                new_ray.maxt = (lRec.p - new_ray.o).norm() - Epsilon;
                new_ray.mint = Epsilon;
                //Test if there is any occlusion.
                if (cos_theta_0 <= 0 || scene->rayIntersect(new_ray)) {
                    Lr_ems = Color3f(0.0f);
                }
                else {
                    pdf_ems_wi_e = Light->pdf(lRec);

                    //BRDF:
                    BSDFQueryRecord bRec_ems(-its.toLocal(recur_ray.d), its.toLocal(lRec.wi), ESolidAngle); //wi(p->camera), wo(p->light), measure
                    bRec_ems.uv = its.uv;
                    Color3f brdf_eval = its.mesh->getBSDF()->eval(bRec_ems);
                    pdf_mats_wi_e = its.mesh->getBSDF()->pdf(bRec_ems);

                    Lr_ems = std::abs(cos_theta_i) * light_L_with_pdf * brdf_eval;
                }

                //Compute Lr
                if (pdf_ems_wi_e + pdf_mats_wi_e > 0) {
                    weight_ems = pdf_ems_wi_e / (pdf_ems_wi_e + pdf_mats_wi_e);
                }
                Lr += weight_ems * mult_t * Lr_ems;
            }
            else {
                hit_specular = true;
            }

            mult_t *= fr_with_pdf_and_cos;

            //update ray.
            recur_ray = Ray3f(its.p, its.toWorld(bRec.wo));
            recur_ray.mint = Epsilon;
        }
        return Lr;
    }

    std::string toString() const {
        return "DirectMIS";
    }
};

NORI_REGISTER_CLASS(PathMIS, "path_mis");
NORI_NAMESPACE_END