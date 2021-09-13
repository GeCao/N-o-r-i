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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DisneyBRDF : public BSDF {
public:
    DisneyBRDF(const PropertyList& propList) {
        /* Base color */
        m_baseColor = propList.getColor("baseColor", Color3f(0.5f));

        /* Subsurface */
        m_subsurface = propList.getFloat("subsurface", 0.0f);

        /* Metallic */
        m_metallic = propList.getFloat("metallic", 0.0f);

        /* Specular */
        m_specular = propList.getFloat("specular", 0.5f);

        /* SpecularTint */
        m_specularTint = propList.getFloat("specularTint", 0.0f);

        /* Roughness */
        m_roughness = propList.getFloat("roughness", 0.5f);

        /* Anisotropic */
        m_anisotropic = propList.getFloat("anisotropic", 0.0f);

        /* Sheen */
        m_sheen = propList.getFloat("sheen", 0.0f);

        /* SheenTint */
        m_sheenTint = propList.getFloat("sheenTint", 0.5f);

        /* Clearcoat */
        m_clearcoat = propList.getFloat("clearcoat", 0.0f);

        /* ClearcoatGloss */
        m_clearcoatGloss = propList.getFloat("clearcoatGloss", 1.0f);

        /* Parameters to discribe roughness and anisotropy */
        m_a = sqrt(1 - 0.9 * m_anisotropic);
        m_alpha_x = fmax(Epsilon, m_roughness * m_roughness / m_a);
        m_alpha_y = fmax(Epsilon, m_roughness * m_roughness * m_a);
        m_alpha = 0.1 + (0.01 - 0.1) * m_clearcoatGloss;

        /* Normalization */
        m_Cdlin = Color3f(pow(m_baseColor.r(), 2.2),
                          pow(m_baseColor.g(), 2.2),
                          pow(m_baseColor.b(), 2.2));
        Clum = 0.3 * m_Cdlin.r() + 0.6 * m_Cdlin.g() + 0.1 * m_Cdlin.b();
        if (Clum > 0)
            m_Ctint = m_Cdlin / Clum;
        else
            m_Ctint = Color3f(1.0f);
    }

    /// Square: pow(_, 2)
    float sqr(const float& x) const { return x * x; }

    /// Interpolation
    float mix(const float& a, const float& b, const float& ratio) const { return a + (b - a) * ratio; }
    Color3f mix(const Color3f& a, const Color3f& b, const float& ratio) const { return a + (b - a) * ratio; }

    /// SchlickFresnel: (1 - cosTheta)^5
    float SchlickFresnel(const float& u)const {
        float m = 1 - u;
        float m2 = m * m;
        return m2 * m2 * m;
    }

    /// GGX
    float smith_GGX(const Vector3f& w) const {
        float tanTheta = Frame::tanTheta(w);
        float A = -0.5 + 0.5 * sqrt(1 + sqr(0.25 * tanTheta));
        return 1/(1 + A);
    }

    /// GGX_aniso
    float smith_GGX_aniso(const Vector3f& w) const {
        float tanTheta = Frame::tanTheta(w),
              cosPhi = Frame::cosPhi(w),
              sinPhi = Frame::sinPhi(w);
        float A = -0.5 + 0.5 * sqrt(1 + (sqr(m_alpha_x * cosPhi) + sqr(m_alpha_y * sinPhi)) * sqr(tanTheta));
        return 1 / (1 + A);
    }
    
    /// GTR2
    Vector3f gtr2Sample(const Point2f& sample) const{
        float x = m_alpha_x * sqrt(sample.x() / (1 - sample.x())) * cos(2 * M_PI * sample.y()),
              y = m_alpha_y * sqrt(sample.x() / (1 - sample.x())) * sin(2 * M_PI * sample.y());
        Vector3f h(-x, -y, 1);
        return h.normalized();
    }

    /// GTR2_pdf (microfacet normal distribution)
    float gtr2Pdf(const Vector3f& wh) const {
        float cosTheta = Frame::cosTheta(wh),
              sinTheta = Frame::sinTheta(wh),
              cosPhi = Frame::cosPhi(wh),
              sinPhi = Frame::sinPhi(wh);
        float aniso = sqr(cosPhi / m_alpha_x) + sqr(sinPhi / m_alpha_y);
        float pdf = 1 / (M_PI * m_alpha_x * m_alpha_y * sqr(sqr(sinTheta) * aniso + sqr(cosTheta)));
        return pdf;
    }

    /// GTR1
    Vector3f gtr1Sample(const Point2f& sample) const {
        Vector3f wh;
        float cosTheta = sqrt((pow(m_alpha, 2 - 2 * sample.x()) - 1) / (m_alpha * m_alpha - 1));
        float sinTheta = sqrt(1 - cosTheta * cosTheta);
        wh.x() = sinTheta * cos(2 * M_PI * sample.y());
        wh.y() = sinTheta * sin(2 * M_PI * sample.y());
        wh.z() = cosTheta;
        return wh;
    }

    /// GTR1_pdf (microfacet normal distribution)
    float gtr1Pdf(const Vector3f& wh) const {
        float cosTheta = Frame::cosTheta(wh);
        float sinTheta = Frame::sinTheta(wh);
        float pdf = (m_alpha * m_alpha - 1) / (2 * M_PI * log(m_alpha) * (sqr(m_alpha * cosTheta) + sqr(sinTheta)));
        return pdf;
    }

    /// f_diffuse
    float fDiffuse(const float& cosThetaI, const float& cosThetaO, const float& cosThetaD) const {
        float fD_90 = 0.5 + 2 * sqr(cosThetaD) * m_roughness;
        return mix(1, fD_90, SchlickFresnel(cosThetaI)) * mix(1, fD_90, SchlickFresnel(cosThetaO));
    }

    /// f_subsurface
    float fSubsurface(const float& cosThetaI, const float& cosThetaO, const float& cosThetaD) const {
        float fSS_90 = sqr(cosThetaD) * m_roughness;
        float fSS = mix(1, fSS_90, SchlickFresnel(cosThetaI)) * mix(1, fSS_90, SchlickFresnel(cosThetaO));
        return 1.25 * (fSS * (1 / (cosThetaI + cosThetaO) - 0.5) + 0.5);
    }

    /// f_sheen
    Color3f fSheen(const float& cosThetaD) const {
        return mix(Color3f(1.0f), m_Ctint, m_sheenTint) * m_sheen * SchlickFresnel(cosThetaD);
    }

    /// Specular
    Color3f specularContrib(const Vector3f& wi, const Vector3f& wo, const Vector3f& wh, const float& cosThetaI, const float& cosThetaO, const float& cosThetaD) const {
        
        /* Frensel */
        Color3f C_specular = mix(Color3f(1.0f), m_Ctint, m_specularTint) * m_specular * 0.08;
        Color3f Cs = mix(C_specular, m_baseColor, m_metallic);
        Color3f Fs = mix(Cs, Color3f(1.0f), SchlickFresnel(cosThetaD));

        /* Shadow factor (Smith GGX2 aniso) */
        float Gs = smith_GGX_aniso(wi) * smith_GGX_aniso(wo);

        /* Macrofacet distribution (GTR2) */
        float Ds = gtr2Pdf(wh);
        return Fs * Gs * Ds / (4 * cosThetaI * cosThetaO);
    }

    /// Clearcoat
    Color3f clearcoatContrib(const Vector3f& wi, const Vector3f& wo, const Vector3f& wh, const float& cosThetaI, const float& cosThetaO, const float& cosThetaD) const {
        
        /* Frensel */
        Color3f Fc = Color3f(0.04f) + Color3f(0.96f) * SchlickFresnel(cosThetaD);

        /* Shadow factor (Smith GGX1) */
        float Gc = smith_GGX(wi) * smith_GGX(wo);

        /* Macrofacet distribution (GTR1) */
        float Dc = gtr1Pdf(wh);
        return Fc * Gc * Dc / (4 * cosThetaI * cosThetaO);
    }


    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord& bRec) const override {
        Normal3f wh = (bRec.wo + bRec.wi).normalized();

        float cosThetaI = Frame::cosTheta(bRec.wi),
              cosThetaO = Frame::cosTheta(bRec.wo),
              cosThetaD = wh.dot(bRec.wi);

        if (cosThetaO <= 0 || cosThetaI <= 0)
            return Color3f(0.0f);

        /* Contrib from diffuse */
        float fD = fDiffuse(cosThetaI, cosThetaO, cosThetaD);
        float fSS = fSubsurface(cosThetaI, cosThetaO, cosThetaD);
        /* Mix the diffuse and subsurface */
        float mixD_SS = mix(fD, fSS, m_subsurface);

        Color3f fSh = fSheen(cosThetaD);

        /* Contrib of specular */
        Color3f sF = specularContrib(bRec.wi, bRec.wo, wh, cosThetaI, cosThetaO, cosThetaD);

        /* Contrib of clearcoat */
        Color3f cF = clearcoatContrib(bRec.wi, bRec.wo, wh, cosThetaI, cosThetaO, cosThetaD);

        Color3f DisneyF = (1 - m_metallic) * (m_baseColor * INV_PI * mixD_SS + fSh) + sF + m_clearcoat / 4 * cF;
        if (DisneyF.minCoeff() <= 0)
            return Color3f(0.0f);
        else
            return DisneyF;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord& bRec) const override {
        Normal3f wh = (bRec.wo + bRec.wi).normalized();
        float Jh = 1 / (4 * wh.dot(bRec.wi));
        if (Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;
        
        /* Diffuse */
        float pDiffuse = Frame::cosTheta(bRec.wo) * INV_PI;
        
        /* Specular(GTR2) */
        float pSpecular = gtr2Pdf(wh) * Frame::cosTheta(wh) * Jh;

        /* Clearcoat(GTR1) */
        float pClearcoat = gtr1Pdf(wh) * Frame::cosTheta(wh) * Jh;

        float prob1 = fmin(0.8, 1 - m_metallic),
              prob2 = 1 / (1 + m_clearcoat);

        float _pdf = prob1 * pDiffuse + (1 - prob1) * (prob2 * pSpecular + (1 - prob2) * pClearcoat);
        if (_pdf <= 0)
            return 0.0f;
        else
            return _pdf;
    }

    /// Sample the BRDF
    virtual Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const override {
        if (bRec.wi.z() <= 0)
            return Color3f(0.0f);
        bRec.measure = ESolidAngle;
        Vector3f wo, wh;
        float prob1 = fmin(0.8, 1 - m_metallic),
              prob2 = 1 / (1 + m_clearcoat);
        float thresh1 = prob1,
              thresh2 = prob1 + prob2 * (1 - prob1);
        Point2f sample = _sample;
        if (sample.x() < thresh1) {
            /* Sample the diffuse */
            sample.x() /= thresh1;
            wo = Warp::squareToCosineHemisphere(sample);
            bRec.wo = wo;
        }
        else if (sample.x() < thresh2) {
            /* Sample the specular */
            sample.x() = (sample.x() - thresh1) / (thresh2 - thresh1);
            wh = gtr2Sample(sample);
            wo = 2 * bRec.wi.dot(wh) * wh - bRec.wi;
            bRec.wo = wo;
        }
        else {
            /* Sample the clearcoat */
            sample.x() = (sample.x() - thresh2) / (1 - thresh2);
            wh = gtr1Sample(sample);
            wo = 2 * bRec.wi.dot(wh) * wh - bRec.wi;
            bRec.wo = wo;
        }

        if (Frame::cosTheta(wo) <= 0 || pdf(bRec) <= 0 || eval(bRec).minCoeff() <= 0) {
            bRec.wo *= 0;
            return Color3f(0.0f);
        }
        else
            return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);

    }

    virtual std::string toString() const override {
        return tfm::format(
            "Disney BRDF[\n"
            "  subsurface = %f,\n"
            "  metallic = %f,\n"
            "  specular = %f,\n"
            "  specularTint = %s,\n"
            "  roughness = %f\n"
            "  anisotropic = %f\n"
            "  sheen = %f\n"
            "  sheenTint = %f\n"
            "  clearcoat = %f\n"
            "  clearcoatGloss = %f\n"
            "]",
            m_subsurface,
            m_metallic,
            m_specular,
            m_specularTint,
            m_roughness,
            m_anisotropic,
            m_sheen,
            m_sheenTint,
            m_clearcoat,
            m_clearcoatGloss
        );
    }
private:
    float m_subsurface;
    float m_metallic;
    float m_specular;
    float m_specularTint;
    float m_roughness;
    float m_anisotropic;
    float m_sheen;
    float m_sheenTint;
    float m_clearcoat;
    float m_clearcoatGloss;
    Color3f m_baseColor;

    float m_a;
    float m_alpha_x;
    float m_alpha_y;
    float m_alpha;
    float Clum;
    Color3f m_Ctint;
    Color3f m_Cdlin;
};

NORI_REGISTER_CLASS(DisneyBRDF, "disney_brdf");
NORI_NAMESPACE_END
