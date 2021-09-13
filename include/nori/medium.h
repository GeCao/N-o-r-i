#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct MediumSamplingRecord {
public:

    /// Traveled distance
    float t;

    /// Location of the scattering interaction
    Point3f p;

    /// Time value associated with the medium scattering event
    float time;

    /// Local particle orientation at \ref p
    Vector3f orientation;

    /**
     * \brief Specifies the transmittance along the segment [mint, t]
     *
     * When sampling a distance fails, this contains the
     * transmittance along the whole ray segment [mint, maxDist].
     */
    float transmittance;

    /// The medium's absorption coefficient at \ref p
    float sigmaA;

    /// The medium's scattering coefficient at \ref p
    float sigmaS;

    /// Records the probability density of sampling a medium interaction at p
    float pdfSuccess;

    /**
     * \brief Records the probability density of sampling a medium
     * interaction in the reverse direction
     *
     * This is essentially the density of obtained by calling \ref sampleDistance,
     * but starting at \c p and stopping at \c ray.o. These probabilities
     * are important for bidirectional methods.
     */
    float pdfSuccessRev;

    /**
     * When the \ref Medium::sampleDistance() is successful, this function
     * returns the probability of \a not having generated a medium interaction
     * until \ref t. Otherwise, it records the probability of
     * not generating any interactions in the whole interval [mint, maxt].
     * This probability is assumed to be symmetric with respect to
     * sampling from the other direction, which is why there is no
     * \c pdfFailureRev field.
     */
    float pdfFailure;

public:
    inline MediumSamplingRecord() { }

    /// Return a string representation
    std::string toString() const;
};

/**
 * \brief Superclass of all mediums
 */
class Medium : public NoriObject {
public:

    virtual std::string toString() const = 0;

    /// Virtual destructor
    virtual ~Medium() { }

    virtual EClassType getClassType() const override { return EMedium; }

    virtual float ComputeSigma_s(float lambda, float diameter, float eta, float rho) const = 0;

    virtual bool GetParamsOnMediumPoint(MediumSamplingRecord& MRec) const = 0;

    virtual Vector3f squareToHenyen_Greenstein(const Point2f sample) const = 0;

    virtual float squareToHenyen_GreensteinPdf(const Vector3f m) const = 0;

    virtual float ComputePhaseFunction(float costheta) const = 0;

    virtual float evalTransmittance(Ray3f& ray) const = 0;

    /*ONLY NEED the parameters: ray.o & ray.d*/
    virtual bool get_tmin_tmax_by_AABBIntersect(Ray3f& ray) const = 0;

    /*
    * tmax_ is the distance between ray.o and its.p
    * sample is used to samle on the volume.
    */
    virtual float freeFlightDistance(Ray3f ray_, float tmax_, float sample) const = 0;
protected:
};

NORI_NAMESPACE_END

#endif /* __NORI_MEDIUM_H */
