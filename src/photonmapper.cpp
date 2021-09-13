/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class PhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    PhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
    }

    virtual void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
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
        int count_photons = 0;
        auto light_list = scene->getLights();
        actual_phtons = 0;
        while (true) {
            for (auto single_light : light_list) {
                actual_phtons += 1;
                Ray3f ray;
                Color3f single_power = single_light->samplePhoton(ray, sampler->next2D(), sampler->next2D());

                Ray3f recur_ray(ray.o, ray.d);
                while (true) {
                    Intersection its;
                    if (!scene->rayIntersect(recur_ray, its))
                        break;
                    Normal3f n = its.shFrame.n;

                    if (its.mesh->getBSDF()->isDiffuse()) {
                        //deposit photon
                        m_photonMap->push_back(Photon(its.p, -recur_ray.d, single_power));
                        count_photons += 1;
                        if (count_photons == m_photonCount) {
                            break;
                        }
                    }

                    //russianRoulette: Termination
                    float succes_probibility = std::min(single_power.maxCoeff(), 0.99f);
                    if (sampler->next1D() <= succes_probibility) {
                        //recursion
                        single_power /= succes_probibility;
                    }
                    else {
                        //Termination
                        break;
                    }

                    //BRDF Term: Nomalize to n = (0, 0, 1)
                    BSDFQueryRecord bRec(-its.toLocal(recur_ray.d), Vector3f(0.0f), ESolidAngle); //wo can be any value, since we will reset it in sample function.
                    bRec.uv = its.uv;
                    Color3f fr_with_pdf_and_cos = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                    single_power *= fr_with_pdf_and_cos;

                    //update ray.
                    recur_ray = Ray3f(its.p, its.toWorld(bRec.wo));
                    recur_ray.mint = Epsilon;
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

		/* Build the photon map */
        m_photonMap->build();
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
    	
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
            if (!scene->rayIntersect(recur_ray, its))
                break;
            Normal3f n = its.shFrame.n;
            //MatrixXf final_G = Transform_G(n);

            //Le Term:
            Color3f Le(0.0f);
            EmitterQueryRecord new_lRec(recur_ray.o, its.p, its.shFrame.n); // ref, p, n
            if (its.mesh->isEmitter()) {
                Le = its.mesh->getEmitter()->eval(new_lRec);
            }
            res_Li += mult_t * Le; // Contribute from material sampling.

            if (its.mesh->getBSDF()->isDiffuse()) {
                //Make Radiance estimate.
                std::vector<uint32_t> results;
                m_photonMap->search(its.p, m_photonRadius, results);
                Color3f photon_contribution(0.0f);
                for (uint32_t i : results) {
                    const Photon& photon = (*m_photonMap)[i];

                    BSDFQueryRecord photon_bRec(-its.toLocal(recur_ray.d), its.toLocal(photon.getDirection()), ESolidAngle); //wo can be any value, since we will reset it in sample function.
                    photon_bRec.uv = its.uv;
                    photon_contribution += its.mesh->getBSDF()->eval(photon_bRec) * photon.getPower();
                }
                photon_contribution /= M_PI * m_photonRadius * m_photonRadius;
                photon_contribution /= actual_phtons;
                res_Li += mult_t * photon_contribution;

                //is diffuse, terminate.
                break;
            }

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
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
