#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMATS : public Integrator {
public:
    PathMATS(const PropertyList& props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f res_Li(0.0f);
        Color3f mult_t(1.0f);
        Ray3f recur_ray(ray.o, ray.d);
        while (true) {
            // Find the surface that is visible in the requested direction
            Intersection its;
            if (!scene->rayIntersect(recur_ray, its))
                break;
            Normal3f n = its.shFrame.n;

            //Le Term:
            Color3f Le(0.0f);
            EmitterQueryRecord new_lRec(recur_ray.o, its.p, its.shFrame.n); // ref, p, n
            if (its.mesh->isEmitter()) {
                Le = its.mesh->getEmitter()->eval(new_lRec);
            }
            res_Li += mult_t * Le; // Contribute from material sampling.

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

            //BRDF Term: Nomalize to n = (0, 0, 1)
            BSDFQueryRecord bRec(-its.toLocal(recur_ray.d), Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
            bRec.uv = its.uv;
            Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

            mult_t *= fr_with_pdf_and_cos;

            //update ray.
            //recur_ray = Ray3f(its.p, final_G.transpose() * bRec.wo);
            recur_ray = Ray3f(its.p, its.toWorld(bRec.wo));
            recur_ray.mint = Epsilon;
        }
        return res_Li;
    }

    std::string toString() const {
        return "PATHMAT";
    }
};

NORI_REGISTER_CLASS(PathMATS, "path_mats");
NORI_NAMESPACE_END