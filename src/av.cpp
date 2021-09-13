#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AverageVisibilityIntegrator : public Integrator {
public:
    AverageVisibilityIntegrator(const PropertyList& props) {
        m_length = props.getFloat("length");
        std::cout << "Length value was : " << m_length << std::endl;
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(1.0f);

        Normal3f n = its.shFrame.n;

        Vector3f v = Warp::sampleUniformHemisphere(sampler, n); //Ray direction.
        Ray3f new_ray(its.p, v);
        Intersection new_its;
        if (scene->rayIntersect(new_ray, new_its)) {
            if (new_its.t < m_length) {
                return Color3f(0.0f);
            }
        }
        return Color3f(1.0f);
    }

    std::string toString() const {
        return "NormalIntegrator[]";
    }
protected:
    float m_length;
};

NORI_REGISTER_CLASS(AverageVisibilityIntegrator, "av");
NORI_NAMESPACE_END