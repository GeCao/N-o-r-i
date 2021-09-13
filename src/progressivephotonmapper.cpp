#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class ProgressivePhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    ProgressivePhotonMapper(const PropertyList& props) {
        /* Lookup parameters */
        m_photonCount = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
        m_alpha = props.getFloat("photonAlpha", 2.0 / 3.0f /* Default: automatic */);
        m_itr = 0;

        /* Create a sample generator for the preprocess step */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }

    virtual void preprocess(const Scene* scene) override {
        if (++m_itr > 1) {
            //update
            m_photonRadius *= sqrt((m_itr - 1.0f + m_alpha) / (m_itr - 1.0f + 1.0f));
            std::cout << "Radius: " << m_photonRadius << ", iteration: " << m_itr << std::endl;
        }
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

        /* Estimate a default photon radius */
        if (m_photonRadius == 0)
            m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;
        if (m_alpha < 0.0f || m_alpha > 1.0f) { m_alpha = 2.0f / 3.0f; }



        /* How to add a photon?
         * m_photonMap->push_back(Photon(
         *	Point3f(0, 0, 0),  // Position
         *	Vector3f(0, 0, 1), // Direction
         *	Color3f(1, 2, 3)   // Power
         * ));
         */

         // put your code to trace photons here
        int count_photons = 0;
        auto light_list = scene->getLights();
        actual_phtons = 0;
        m_photonMap->clear();
        while (true) {
            for (auto single_light : light_list) {
                actual_phtons += 1;
                Ray3f ray;
                Color3f single_power = single_light->samplePhoton(ray, m_sampler->next2D(), m_sampler->next2D());

                Ray3f recur_ray(ray.o, ray.d);
                while (true) {
                    Intersection its;
                    if (!scene->rayIntersect(recur_ray, its)) { break; }

                    if (its.mesh->getBSDF()->isDiffuse()) {
                        //deposit photon
                        m_photonMap->push_back(Photon(its.p, -recur_ray.d, single_power));
                        if (++count_photons == m_photonCount) { break; }
                    }

                    //russianRoulette: Termination
                    float succes_probibility = std::min(single_power.maxCoeff(), 0.99f);
                    if (m_sampler->next1D() <= succes_probibility) { single_power /= succes_probibility; } //recursion
                    else { break; } //Termination

                    //BRDF Term: Nomalize to n = (0, 0, 1)
                    BSDFQueryRecord bRec(-its.toLocal(recur_ray.d), Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
                    bRec.uv = its.uv;
                    Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, m_sampler->next2D());

                    single_power *= fr_with_pdf_and_cos;

                    //update ray.
                    recur_ray = Ray3f(its.p, its.toWorld(bRec.wo));
                    recur_ray.mint = Epsilon;
                }
                //Trace Photons
                if (count_photons == m_photonCount) { break; }
            }
            //Trace Photons
            if (count_photons == m_photonCount) { break; }
        }
        actual_phtons /= light_list.size();

        /* Build the photon map */
        m_photonMap->build();
    }

    virtual Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& _ray) const override {

        /* How to find photons?
         * std::vector<uint32_t> results;
         * m_photonMap->search(Point3f(0, 0, 0), // lookup position
         *                     m_photonRadius,   // search radius
         *                     results);
         *
         * for (uint32_t i : results) {
         *    const Photon &photon = (*m_photonMap)[i];
         *    cout << "Found photon!" << endl;
         *    cout << " Position  : " << photon.getPosition().toString() << endl;
         *    cout << " Power     : " << photon.getPower().toString() << endl;
         *    cout << " Direction : " << photon.getDirection().toString() << endl;
         * }
         */

         // put your code for path tracing with photon gathering here
        Color3f res_Li(0.0f);
        Color3f mult_t(1.0f);
        Ray3f recur_ray(_ray.o, _ray.d);
        while (true) {
            // Find the surface that is visible in the requested direction
            Intersection its;
            if (!scene->rayIntersect(recur_ray, its)) { break; }

            //Le Term:
            res_Li += mult_t * RequestLe_OnDiffuseSurface(recur_ray, &its);

            if (its.mesh->getBSDF()->isDiffuse()) {
                //Make Radiance estimate.
                res_Li += mult_t * PhotonDistribution_OnDiffuseSurface(recur_ray, &its);
                //is diffuse, terminate.
                break;
            }

            //russianRoulette: Termination
            float succes_probibility = std::min(mult_t.maxCoeff(), 0.99f);
            if (sampler->next1D() <= succes_probibility) { mult_t /= succes_probibility; } //recursion
            else { break; } //Termination

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

    virtual std::string toString() const override {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }

private:
    Color3f RequestLe_OnDiffuseSurface(Ray3f ray, Intersection* its) const {
        //Le Term:
        Color3f Le(0.0f);
        EmitterQueryRecord lRec(ray.o, its->p, its->shFrame.n); // ref, p, n

        if (its->mesh->isEmitter()) {
            Le = its->mesh->getEmitter()->eval(lRec);
        }
        return Le; // Contribute from material sampling.
    }

    Color3f PhotonDistribution_OnDiffuseSurface(Ray3f ray, Intersection* its) const {
        //Make Radiance estimate.
        std::vector<uint32_t> results;
        m_photonMap->search(its->p, m_photonRadius, results);
        Color3f photon_contribution(0.0f);
        for (uint32_t i : results) {
            const Photon& photon = (*m_photonMap)[i];

            BSDFQueryRecord bRec(-its->toLocal(ray.d), its->toLocal(photon.getDirection()), ESolidAngle); //wo can be any value, since we will reset it in sample function.
            bRec.uv = its->uv;
            float sigma2 = m_photonRadius * m_photonRadius / 9.0f;
            photon_contribution += its->mesh->getBSDF()->eval(bRec) * photon.getPower() / \
                (2.0 * M_PI * sigma2) * exp(-0.5 * (photon.getPosition() - its->p).squaredNorm() / sigma2);
        }
        //photon_contribution /= M_PI * m_photonRadius * m_photonRadius;
        photon_contribution /= actual_phtons;
        return photon_contribution;
    }

private:
    /*
     * Important: m_photonCount is the total number of photons deposited in the photon map,
     * NOT the number of emitted photons. You will need to keep track of those yourself.
     */
    int m_photonCount;
    float m_photonRadius;
    float m_alpha;
    int m_itr;
    int actual_phtons;
    std::unique_ptr<PhotonMap> m_photonMap;

    Sampler* m_sampler;
};

NORI_REGISTER_CLASS(ProgressivePhotonMapper, "progressivephotonmapper");
NORI_NAMESPACE_END
