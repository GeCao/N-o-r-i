#include <nori/camera.h>
#include <nori/rfilter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

/**
 * \brief Perspective camera with depth of field
 *
 * This class implements a simple perspective camera model. It uses an
 * infinitesimally small aperture, creating an infinite depth of field.
 */
    class PerspectiveDepthCamera : public Camera {
    public:
        PerspectiveDepthCamera(const PropertyList& propList) {
            /* Width and height in pixels. Default: 720p */
            m_outputSize.x() = propList.getInteger("width", 1280);
            m_outputSize.y() = propList.getInteger("height", 720);
            m_invOutputSize = m_outputSize.cast<float>().cwiseInverse();

            /* Specifies an optional camera-to-world transformation. Default: none */
            m_cameraToWorld = propList.getTransform("toWorld", Transform());

            /* Horizontal field of view in degrees */
            m_fov = propList.getFloat("fov", 30.0f);

            /* Near and far clipping planes in world-space units */
            m_nearClip = propList.getFloat("nearClip", 1e-4f);
            m_farClip = propList.getFloat("farClip", 1e4f);

            /* Depth of field parameters */
            m_focal_distance = propList.getFloat("focaldistance", 0.0f); //光圈与film之间的距离：0 for default
            m_aperture_radius = propList.getFloat("apertureradius", 0.2f);//光圈的大小: 

            m_focal_length = propList.getFloat("focallength", 0.050f); //焦距: 50mm for default.

            m_rfilter = NULL;
        }

        virtual void activate() override {
            float aspect = m_outputSize.x() / (float)m_outputSize.y();

            /* Project vectors in camera space onto a plane at z=1:
             *
             *  xProj = cot * x / z
             *  yProj = cot * y / z
             *  zProj = (far * (z - near)) / (z * (far-near))
             *  The cotangent factor ensures that the field of view is
             *  mapped to the interval [-1, 1].
             */
            float recip = 1.0f / (m_farClip - m_nearClip),
                cot = 1.0f / std::tan(degToRad(m_fov / 2.0f));

            Eigen::Matrix4f perspective;
            perspective <<
                cot, 0, 0, 0,
                0, cot, 0, 0,
                0, 0, m_farClip* recip, -m_nearClip * m_farClip * recip,
                0, 0, 1, 0;

            /**
             * Translation and scaling to shift the clip coordinates into the
             * range from zero to one. Also takes the aspect ratio into account.
             */
            m_sampleToCamera = Transform(
                Eigen::DiagonalMatrix<float, 3>(Vector3f(0.5f, -0.5f * aspect, 1.0f)) *
                Eigen::Translation<float, 3>(1.0f, -1.0f / aspect, 0.0f) * perspective).inverse();

            /* If no reconstruction filter was assigned, instantiate a Gaussian filter */
            if (!m_rfilter) {
                m_rfilter = static_cast<ReconstructionFilter*>(
                    NoriObjectFactory::createInstance("gaussian", PropertyList()));
                m_rfilter->activate();
            }
        }

        Color3f sampleRay(Ray3f& ray,
            const Point2f& samplePosition,
            const Point2f& apertureSample) const {
            /* Compute the corresponding position on the
               near plane (in local camera space) */
            Point3f nearP = m_sampleToCamera * Point3f(
                samplePosition.x() * m_invOutputSize.x(),
                samplePosition.y() * m_invOutputSize.y(), 0.0f);

            /* Turn into a normalized ray direction, and
               adjust the ray interval accordingly */
            Vector3f d = nearP.normalized();
            float invZ = 1.0f / d.z();

            ray.o = m_cameraToWorld * Point3f(0, 0, 0);
            ray.d = m_cameraToWorld * d;
            ray.mint = m_nearClip * invZ;
            ray.maxt = m_farClip * invZ;

            //Make depth of filed activate:
            /*
            Vector3f n_in_CameraSpace = Vector3f(0.0, 0.0, 1.0f);
            if (nearP.dot(n_in_CameraSpace) < 0) { n_in_CameraSpace = -n_in_CameraSpace; }
            Vector3f n_in_WorldSpace = (m_cameraToWorld * n_in_CameraSpace).normalized();
            Vector3f s, t;
            coordinateSystem(n_in_WorldSpace, s, t);
            Point3f center_aperture = ray.o + m_focal_distance * n_in_WorldSpace;
            Point3f sample_aperture = center_aperture + m_aperture_radius * (apertureSample.x() * s + apertureSample.y() * t);
            Point3f focus_point = ray.o + (nearP.dot(n_in_CameraSpace)) * ray.d;
            ray.o = sample_aperture;
            ray.d = (focus_point - ray.o).normalized();
            */
            //Sample points on lens:
            Point3f pLens(m_aperture_radius * apertureSample.x(), m_aperture_radius * apertureSample.y(), 0.0);
            float ft = m_focal_distance / d.z();
            Point3f pFocus = ft * d;
            ray.o = m_cameraToWorld * pLens;
            Vector3f d_norm = (pFocus - pLens).normalized();
            ray.d = m_cameraToWorld * d_norm;

            ray.update();

            return Color3f(1.0f);
        }

        virtual void addChild(NoriObject* obj) override {
            switch (obj->getClassType()) {
            case EReconstructionFilter:
                if (m_rfilter)
                    throw NoriException("Camera: tried to register multiple reconstruction filters!");
                m_rfilter = static_cast<ReconstructionFilter*>(obj);
                break;

            default:
                throw NoriException("Camera::addChild(<%s>) is not supported!",
                    classTypeName(obj->getClassType()));
            }
        }

        /// Return a human-readable summary
        virtual std::string toString() const override {
            return tfm::format(
                "PerspectiveDepthCamera[\n"
                "  cameraToWorld = %s,\n"
                "  outputSize = %s,\n"
                "  fov = %f,\n"
                "  clip = [%f, %f],\n"
                "  rfilter = %s\n"
                "]",
                indent(m_cameraToWorld.toString(), 18),
                m_outputSize.toString(),
                m_fov,
                m_nearClip,
                m_farClip,
                indent(m_rfilter->toString())
            );
        }
    private:
        Vector2f m_invOutputSize;
        Transform m_sampleToCamera;
        Transform m_cameraToWorld;
        float m_fov;
        float m_nearClip;
        float m_farClip;

        float m_focal_distance;
        float m_aperture_radius;
        float m_focal_length;
};

NORI_REGISTER_CLASS(PerspectiveDepthCamera, "perspective_depth");
NORI_NAMESPACE_END