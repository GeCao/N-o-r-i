#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter{
public:
	PointLight(const PropertyList& propList) {
		m_power = propList.getColor("power", Color3f());
		m_position = propList.getPoint3("position", Point3f());
	}

	virtual Color3f sample(EmitterQueryRecord& lRec, const Point2f& sample) const override {
		lRec.wi = (m_position - lRec.ref).normalized();
		lRec.p = m_position;
		lRec.n = -lRec.wi;
		return eval(lRec);
	}

	virtual Color3f eval(const EmitterQueryRecord& lRec) const override {
		Color3f I = m_power / 4.0 / M_PI;
		return I / (m_position - lRec.ref).squaredNorm();
	}

	virtual float pdf(const EmitterQueryRecord& lRec) const override {
		return 1.0f;
	}

	virtual std::string toString() const {
		return "Point Light";
	}

protected:
	Color3f m_power; //Watt
	Point3f m_position; //World space position.
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END