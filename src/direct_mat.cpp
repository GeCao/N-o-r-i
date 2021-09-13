#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectMAT : public Integrator {
public:
    DirectMAT(const PropertyList& props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        Normal3f n = its.shFrame.n;

        //Le Term:
        Color3f Le(0.0f);
        EmitterQueryRecord new_lRec(ray.o, its.p, n); // ref, p, n
        if (its.mesh->isEmitter()) {
            Le = its.mesh->getEmitter()->eval(new_lRec);
        }

        //BRDF Term: Nomalize to n = (0, 0, 1)
        MatrixXf final_G = Transform_G(n);
        //    wi(p -> light), wo(p -> camera), measure
        BSDFQueryRecord bRec(-final_G * ray.d, Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
        bRec.uv = its.uv;
        Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

        //L Term:
        EmitterQueryRecord lRec(its.p);
        lRec.wi = final_G.transpose() * bRec.wo;
        Color3f light_L(0.0f);
        Intersection emit_its;
        Ray3f new_ray = Ray3f(its.p, lRec.wi);
        new_ray.mint = Epsilon;
        //->Test if there is any intersection.
        if (scene->rayIntersect(new_ray, emit_its)) {
            if (emit_its.mesh->isEmitter()) {
                //Intersect and it is a light
                lRec.p = emit_its.p;
                lRec.n = emit_its.shFrame.n;
                light_L = emit_its.mesh->getEmitter()->eval(lRec);
            }
            else {
                //Intersect but it is not light
                return Le;
            }
        }
        else {
            //No intersect
            return Le;
        }

        //Test if the ray is under the intersection plane.
        float cos_theta_i = lRec.wi.dot(n);
        float cos_theta_0 = lRec.n.dot(-lRec.wi);
        if (cos_theta_i <= 0 || cos_theta_0 <= 0) {
            return Le;
        }

        return Le + fr_with_pdf_and_cos * light_L;
    }

    std::string toString() const {
        return "DirectMAT";
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

NORI_REGISTER_CLASS(DirectMAT, "direct_mats");
NORI_NAMESPACE_END