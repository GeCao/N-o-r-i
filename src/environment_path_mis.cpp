#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class EnvironmentPathtMis : public Integrator {
public:
    EnvironmentPathtMis(const PropertyList& props) {

    }

    virtual void preprocess(const Scene* scene) override {

        std::vector<Emitter*> light;
        light = scene->getLights();
        int numLight = light.size();
        for (int i = 0; i < numLight; i++) {
            if (light[i]->isEnvEmitter())
                light[i]->activate(scene);
        }

    }


    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {

        /*Get all the light sources*/
        std::vector<Emitter*> light = scene->getLights();
        int numLight = light.size();

        /*Compute indirect illumation: Li*/
        Color3f Li(0.0f);
        Color3f throughput(1.0f);
        int bounce = 0;

        /*Initialize the first ray*/
        Point3f xo = ray.o;
        Vector3f wi = ray.d;

        /*Outside the loop, store the information of last iteration*/
        EMeasure isDiscrete = EUnknownMeasure;
        float pdf_mat = 0.0f;

        /*Initialize the record used to sample direction*/
        BSDFQueryRecord bRec(-wi);
        bRec.measure = EDiscrete;

        while (true) {
            bounce++;
            Intersection test;
            Ray3f rayI(xo, wi);
            float w_mat, w_ems;
            if (!scene->rayIntersect(rayI, test)) {
                float pdf;
                Color3f env = scene->getEnv(rayI);
                if (bounce == 0)
                    return env;
                Li += env * throughput;
                break;
            }
                
            xo = test.p;

            /*Add material sample contribution*/
            if (test.mesh->isEmitter()) {
                EmitterQueryRecord eRec(rayI.o, xo, test.shFrame.n);
                ShapeQueryRecord emtShape;
                eRec.pdf = test.mesh->pdfSurface(emtShape);
                if (bRec.measure == EDiscrete)
                    w_mat = 1.0f;
                else
                    w_mat = pdf_mat / (pdf_mat + test.mesh->getEmitter()->pdf(eRec) / numLight);
                Li += w_mat * test.mesh->getEmitter()->eval(eRec) * throughput;
            }

            /*Russian Roulette*/
            float q = (throughput.maxCoeff() < 0.99 ? throughput.maxCoeff() : 0.99);
            if (bounce >= 3) {
                if (sampler->next1D() >= q)
                    break;
                else
                    throughput /= q;
            }

            /*Sample bsdf on hemisphere on x0*/
            bRec.wi = test.toLocal(-rayI.d);
            Color3f textureI = test.mesh->getBSDF()->sample(bRec, sampler->next2D());
            isDiscrete = bRec.measure;
            pdf_mat = test.mesh->getBSDF()->pdf(bRec);
            wi = test.toWorld(bRec.wo);

            /*Add emitter sample contribution*/
            if (bRec.measure != EDiscrete) {
                EmitterQueryRecord lRec(xo);
                int randomIndex = sampler->next1D() * numLight;
                Color3f Le = light[randomIndex]->sample(lRec, sampler->next2D());
                Intersection emitterSample;
                bool flag = scene->rayIntersect(lRec.shadowRay, emitterSample);
                if (!flag && test.shFrame.n.dot(lRec.wi) > 0 && lRec.n.dot(-lRec.wi) > 0) {
                    BSDFQueryRecord bRecE(test.toLocal(-rayI.d), test.toLocal(lRec.wi), ESolidAngle);
                    w_ems = light[randomIndex]->pdf(lRec) / numLight / (light[randomIndex]->pdf(lRec) / numLight + test.mesh->getBSDF()->pdf(bRecE));
                    Li += w_ems * throughput * Le * test.mesh->getBSDF()->eval(bRecE) * test.shFrame.n.dot(lRec.wi) * numLight;
                }
            }

            throughput *= textureI;
        }
        return Li;
    }

    std::string toString() const {
        return "PathMis[]";
    }


};

NORI_REGISTER_CLASS(EnvironmentPathtMis, "environment_path_mis");
NORI_NAMESPACE_END