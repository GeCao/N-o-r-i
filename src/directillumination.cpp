#include <nori/integrator.h>
#include <nori/warp.h>
#include <nori/scene.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectIllusionIntegrator : public Integrator {
public:
    DirectIllusionIntegrator(const PropertyList& props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its1;
        if (!scene->rayIntersect(ray, its1))
            return Color3f(0.0f);

        /* Return the component-wise absolute
           value of the shading normal as a color */
        Normal3f n = its1.shFrame.n.cwiseAbs();
        EmitterQueryRecord lRec(its1.p);
        Point2f sample_get = sampler->next2D();
        Color3f light_L = scene->getLights()[0]->sample(lRec, sample_get);

        //new generated ray and try to find the visibility.
        Ray3f new_ray = ray;
        new_ray.o = its1.p;
        new_ray.d = lRec.wi;
        new_ray.maxt = (lRec.p - new_ray.o).norm();
        new_ray.mint = 0;
        float V = 1;
        float cos_theta_i_k = new_ray.d.dot(n);
        /*
        if (scene->rayIntersect(new_ray)) {
            V = 0;
        }
        if (cos_theta_i_k <= 0) {
            V = 0;
        }
        */
        float Lr = V * cos_theta_i_k;

        BSDFQueryRecord bRec(-ray.d, new_ray.d, ESolidAngle);
        bRec.uv = its1.uv;

        //Lr *= scene->getShapes()[0]->getBSDF()->pdf(bRec);
        return Lr * light_L * its1.mesh->getBSDF()->eval(bRec);
    }

    std::string toString() const {
        return "DirectIllusionIntegrator[]";
    }
};

NORI_REGISTER_CLASS(DirectIllusionIntegrator, "direct");
NORI_NAMESPACE_END