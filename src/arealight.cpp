/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

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

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
    AreaEmitter(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }

    virtual std::string toString() const override {
        return tfm::format(
                "AreaLight[\n"
                "  radiance = %s,\n"
                "]",
                m_radiance.toString());
    }

    // A function can be called by the integrators whenever a ray intersects your mesh emitter
    virtual Color3f eval(const EmitterQueryRecord & lRec) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        if ((lRec.p - lRec.ref).dot(lRec.n) < 0) { // Front
            return m_radiance;
        }
        else { // Back
            return Color3f(0.0f);
        }
    }

    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        //Create enquiry.
        ShapeQueryRecord sRec;
        //Sample the pdf into sRec.pdf, Sample the n into SRec.n, Sample the p into sRec.p;
        m_shape->sampleSurface(sRec, sample);
        //From sRec -> lRec
        lRec.p = sRec.p; lRec.n = sRec.n;
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.pdf = pdf(lRec);

        double cos_theta_0 = -lRec.wi.dot(lRec.n);
        double dis2 = (lRec.p - lRec.ref).squaredNorm();
        if (lRec.pdf == 0) {
            return Color3f(0.0f);
        }
        return eval(lRec) / lRec.pdf / dis2 * cos_theta_0;
    }

    virtual float pdf(const EmitterQueryRecord &lRec) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        //Create enquiry.
        ShapeQueryRecord sRec;
        //get pdf.
        double cos_theta_0 = -lRec.wi.dot(lRec.n);
        if (cos_theta_0 > 0) { // Front
            double dis2 = (lRec.p - lRec.ref).squaredNorm();
            return m_shape->pdfSurface(sRec);
        }
        else { // Back
            return 0.0f;
        }
    }


    virtual Color3f samplePhoton(Ray3f &ray, const Point2f &sample1, const Point2f &sample2) const override {
        if (!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        //Create enquiry.
        ShapeQueryRecord sRec;
        //Sample the pdf into sRec.pdf = 1/A, Sample the n into SRec.n, Sample the p into Srec.p;
        m_shape->sampleSurface(sRec, sample1);
        Vector3f w = Warp::squareToCosineHemisphere(sample2);
        MatrixXf final_G = Transform_G(sRec.n);
        w = final_G.transpose() * w;
        ray.d = w; ray.o = sRec.p;
        return M_PI * m_radiance / sRec.pdf;
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

protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END