#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/medium.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class VolumetricPhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    VolumetricPhotonMapper(const PropertyList& props) {
        /* Lookup parameters */
        m_photonCount = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);

        actual_phtons = 0;
    }

    virtual void preprocess(const Scene* scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler* sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

        /* Estimate a default photon radius */
        if (m_photonRadius == 0)
            m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;


        /* How to add a photon?
         * m_photonMap->push_back(Photon(
         *	Point3f(0, 0, 0),  // Position
         *	Vector3f(0, 0, 1), // Direction
         *	Color3f(1, 2, 3)   // Power
         * ));
         */

         // put your code to trace photons here
        m_medium = scene->getMedium();

        int count_photons = 0;
        auto light_list = scene->getLights();
        actual_phtons = 0;
        int medium_photons = 0;
        while (true) {
            for (auto single_light : light_list) { // Iteration the whole emitter list.
                actual_phtons += 1;
                Ray3f ray;
                Color3f single_power = single_light->samplePhoton(ray, sampler->next2D(), sampler->next2D());

                Ray3f recur_ray(ray.o, ray.d);
                while (true) {
                    Intersection its;
                    float tmax;
                    bool Intersected = scene->rayIntersect(recur_ray, its);
                    if (!Intersected) {
                        if (m_medium->get_tmin_tmax_by_AABBIntersect(recur_ray)) {
                            its.p = (recur_ray.maxt) * recur_ray.d;
                            tmax = recur_ray.maxt;
                        }
                        else { break; }
                    }
                    else {
                        tmax = its.t;
                    }
                    //Normal3f n = its.shFrame.n;

                    float t = m_medium->freeFlightDistance(recur_ray, tmax, sampler->next1D());
                    if (t < tmax) {//Media scattering
                        medium_photons += 1;
                        //propagate photon
                        Point3f xt = recur_ray.o + t * recur_ray.d;

                        //deposit photon
                        m_photonMap->push_back(Photon(xt, -recur_ray.d, single_power));
                        if (++count_photons == m_photonCount) { break; }

                        MediumSamplingRecord MRec;
                        MRec.p = xt;
                        m_medium->GetParamsOnMediumPoint(MRec);
                        float sigma_t = MRec.sigmaS + MRec.sigmaA;
                        if (sigma_t != 0.0f) { single_power *= MRec.sigmaS / sigma_t; }

                        //update ray.
                        recur_ray = Ray3f(xt, SamplePF(recur_ray.d, sampler->next2D())); // Need to implement SamplePF by CaoGe
                        recur_ray.mint = Epsilon;
                    }
                    else {
                        if (!Intersected) { break; }
                        if (its.mesh->getBSDF()->isDiffuse()) { //deposit photon
                            m_photonMap->push_back(Photon(its.p, -recur_ray.d, single_power));
                            if (++count_photons == m_photonCount) { break; }
                        }

                        //russianRoulette: Termination
                        float succes_probibility = std::min(single_power.maxCoeff(), 0.99f);
                        if (sampler->next1D() <= succes_probibility) { single_power /= succes_probibility; } //recursion
                        else { break; } //Termination

                        //BRDF Term: Nomalize to n = (0, 0, 1)
                        BSDFQueryRecord bRec(-its.toLocal(recur_ray.d), Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
                        bRec.uv = its.uv;
                        Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                        single_power *= fr_with_pdf_and_cos;

                        //update ray.
                        recur_ray = Ray3f(its.p, its.toWorld(bRec.wo));
                        recur_ray.mint = Epsilon;
                    }
                }
                //Trace Photons
                if (count_photons == m_photonCount) {
                    break;
                }
            }
            //Trace Photons
            if (count_photons == m_photonCount) {
                break;
            }
        }
        actual_phtons /= light_list.size();
        std::cout << "Medium Photons: " << medium_photons << std::endl;

        /* Build the photon map */
        m_photonMap->build();
    }

    Vector3f SamplePF(Normal3f n, Point2f sample) const {
        Frame shFrame(n);
        Vector3f local_wi = m_medium->squareToHenyen_Greenstein(sample);
        float pdf = m_medium->squareToHenyen_GreensteinPdf(local_wi);
        Vector3f world_wi = shFrame.toWorld(local_wi);
        return world_wi;
    }

    Color3f PhotonDistribution_OnVolume(Ray3f ray, Point3f its_p, float& Tr) const {
        //Make Radiance estimate.
        Color3f photon_contribution(0.0f);
        ray.mint = Epsilon; ray.maxt = (its_p - ray.o).norm(); // Costrain the routine.
        float tmin = ray.mint, tmax = ray.maxt;
        if (!m_medium->get_tmin_tmax_by_AABBIntersect(ray)) { return Color3f(0.0f); } // Medium not detected
        tmin = std::max(tmin, ray.mint); tmax = std::min(tmax, ray.maxt);
        if (tmax <= tmin) { return Color3f(0.0f); } // Medium Detected but it is not in out routine.

        //Medium Detected
        int res = 128;
        float dt = (tmax - tmin) / (float)(res - 1), sigma_s = 0.0f, phasefunc = 0.0f;
        MediumSamplingRecord MRec;
        MRec.p = ray.o + tmin * ray.d;
        Tr = 1.0f;
        for (int i = 0; i < res; i++) {
            MRec.p = MRec.p + (i * dt) * ray.d;
            m_medium->GetParamsOnMediumPoint(MRec);
            Tr *= exp(-dt * (MRec.sigmaA + MRec.sigmaS));
            sigma_s = MRec.sigmaS;

            //Find the photons.
            std::vector<uint32_t> results;
            m_photonMap->search(MRec.p, m_photonRadius, results);
            Color3f res_Ls(0.0f);
            if (results.size() != 0) {
                for (uint32_t k : results) {
                    const Photon& photon = (*m_photonMap)[k];
                    Point3f xt_i = photon.getPosition();
                    float costheta = ray.d.dot(-photon.getDirection());
                    phasefunc = m_medium->ComputePhaseFunction(costheta);
                    res_Ls += phasefunc * photon.getPower();
                }
                res_Ls /= 4.0f / 3.0f * M_PI * pow(m_photonRadius, 3);
                res_Ls /= results.size();
            }

            photon_contribution += Tr * sigma_s * res_Ls * dt;
        }
        
        return photon_contribution;
    }

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

            BSDFQueryRecord photon_bRec(-its->toLocal(ray.d), its->toLocal(photon.getDirection()), ESolidAngle); //wo can be any value, since we will reset it in sample function.
            photon_bRec.uv = its->uv;
            photon_contribution += its->mesh->getBSDF()->eval(photon_bRec) * photon.getPower();
        }
        photon_contribution /= M_PI * m_photonRadius * m_photonRadius;
        photon_contribution /= actual_phtons;
        return photon_contribution;
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
            bool Intersected = scene->rayIntersect(recur_ray, its);
            if (!Intersected) {
                if (m_medium->get_tmin_tmax_by_AABBIntersect(recur_ray)) {
                    its.p = (recur_ray.maxt) * recur_ray.d;
                }
                else { break; }
            }
            Normal3f n = its.shFrame.n;

            float Tr = 1.0;
            res_Li += mult_t * PhotonDistribution_OnVolume(recur_ray, its.p, Tr);
            mult_t *= Tr;

            if (!Intersected) {
                break;
            }

            //Le Term:
            res_Li += mult_t * RequestLe_OnDiffuseSurface(recur_ray, &its);

            if (its.mesh->getBSDF()->isDiffuse()) {
                //Make Radiance estimate.
                res_Li += mult_t * PhotonDistribution_OnDiffuseSurface(recur_ray, &its);
                //is diffuse, terminate.
                break;
            }
            else {
                //russianRoulette: Termination
                float succes_probibility = std::min(mult_t.maxCoeff(), 0.99f);
                if (sampler->next1D() <= succes_probibility) { mult_t /= succes_probibility; } //recursion
                else { break; } //Termination
            }

            //BRDF Term: Nomalize to n = (0, 0, 1)
            BSDFQueryRecord bRec(-its.toLocal(recur_ray.d), Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
            bRec.uv = its.uv;
            Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

            mult_t *= fr_with_pdf_and_cos;

            //update ray.
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
    /*
     * Important: m_photonCount is the total number of photons deposited in the photon map,
     * NOT the number of emitted photons. You will need to keep track of those yourself.
     */
    int m_photonCount;
    float m_photonRadius;
    int actual_phtons;
    std::unique_ptr<PhotonMap> m_photonMap;

    const Medium* m_medium = nullptr;
};

NORI_REGISTER_CLASS(VolumetricPhotonMapper, "volumetricphotonmapper");
NORI_NAMESPACE_END