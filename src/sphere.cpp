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

#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Sphere : public Shape {
public:
    Sphere(const PropertyList & propList) {
        m_position = propList.getPoint3("center", Point3f());
        m_radius = propList.getFloat("radius", 1.f);

        m_bbox.expandBy(m_position - Vector3f(m_radius));
        m_bbox.expandBy(m_position + Vector3f(m_radius));
    }

    virtual BoundingBox3f getBoundingBox(uint32_t index) const override { return m_bbox; }

    virtual Point3f getCentroid(uint32_t index) const override { return m_position; }

    /*
    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {
        Vector3f vector_o_c = m_position - ray.o;
        float closet_ray_dist = vector_o_c.dot(ray.d);
        float D_2 = vector_o_c.squaredNorm() - closet_ray_dist * closet_ray_dist;
        if (D_2 < 0 || D_2 > m_radius*m_radius) { return false; }
        
        //if (vector_o_c.norm() > m_radius) {
            //t = closet_ray_dist - sqrt(m_radius * m_radius - D_2); // The ray.o is outside the sphere
        //}
        //else {
            //t = closet_ray_dist + sqrt(m_radius * m_radius - D_2); // The ray.o is inside the sphere
        //}
        
        t = closet_ray_dist - sqrt(m_radius * m_radius - D_2); // The ray.o is outside the sphere
        if (t >= ray.mint && t <= ray.maxt) { return true; }
        t = closet_ray_dist + sqrt(m_radius * m_radius - D_2); // The ray.o is inside the sphere
        if (t >= ray.mint && t <= ray.maxt) { return true; }
        return false;
    }

    virtual void setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const override {
        float t, u, v;
        if (!rayIntersect(index, ray, u, v, t)) {
            return;
        }
        its.p = ray.o + ray.d * t;
        Vector3f n = (its.p - m_position).normalized();
        its.shFrame.n = n;
        its.geoFrame.n = n;
        
        float theta = acos(n.z()); //get [0, pi]
        float fai = acos(n.x() / sin(theta)); //get [0, pi]
        if ((n.y() / sin(theta)) < 0) { fai = -fai + 2 * M_PI; } //sin(fai) < 0
        v = theta / M_PI;
        u = fai / (2 * M_PI);
        its.uv = Point2f(u, v);
        
        Vector3f Vec_s, Vec_t;
        coordinateSystem(its.shFrame.n, Vec_s, Vec_t);
        its.shFrame.t = Vec_t;
        its.geoFrame.t = Vec_t;
        its.shFrame.s = Vec_s;
        its.geoFrame.s = Vec_s;
    }
    */
    virtual bool rayIntersect(uint32_t index, const Ray3f& ray, float& u, float& v, float& t) const override {

        float a, b, c, delta, temp1, temp2;
        a = pow(ray.d(0), 2) + pow(ray.d(1), 2) + pow(ray.d(2), 2);
        b = 2 * ((ray.o(0) - m_position(0)) * ray.d(0) + (ray.o(1) - m_position(1)) * ray.d(1) + (ray.o(2) - m_position(2)) * ray.d(2));
        c = pow((ray.o(0) - m_position(0)), 2) + pow((ray.o(1) - m_position(1)), 2)
            + pow((ray.o(2) - m_position(2)), 2) - pow(m_radius, 2);
        delta = pow(b, 2) - 4 * a * c;
        if (delta < 0)
            return false;
        else {
            temp1 = (-b - sqrt(delta)) / (2 * a);
            temp2 = (-b + sqrt(delta)) / (2 * a);
        }

        bool temp1b, temp2b;
        temp1b = false;
        temp2b = false;


        if (temp1 >= ray.mint && temp1 < ray.maxt)
            temp1b = true;
        if (temp2 >= ray.mint && temp2 < ray.maxt)
            temp2b = true;

        if (temp1b == false && temp2b == false)
            return false;
        else if (temp1b == false && temp2b == true)
            t = temp2;
        else
            t = temp1;
        return true;

    }

    virtual void setHitInformation(uint32_t index, const Ray3f& ray, Intersection& its) const override {

        its.p = ray.o + ray.d * its.t;
        its.shFrame.n = (its.p - m_position).normalized();
        its.geoFrame.n = its.shFrame.n;
        Vector3f t, s;
        coordinateSystem(its.shFrame.n, s, t);
        its.shFrame.t = t;
        its.geoFrame.t = t;
        its.shFrame.s = s;
        its.geoFrame.s = s;
        float u, v;
        Point3f localPoint = its.p - m_position;

        float Theta = atan2(localPoint(1), localPoint(0));
        if (Theta >= 0)
            u = Theta / 2 * INV_PI;
        else
            u = (2 * M_PI + Theta) / 2 * INV_PI;
        v = acos(localPoint(2) / m_radius) * INV_PI;
        its.uv = Point2f(1 - u, v);
    }

    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override {
        Vector3f q = Warp::squareToUniformSphere(sample);
        sRec.p = m_position + m_radius * q;
        sRec.n = q;
        sRec.pdf = std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override {
        return std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }


    virtual std::string toString() const override {
        return tfm::format(
                "Sphere[\n"
                "  center = %s,\n"
                "  radius = %f,\n"
                "  bsdf = %s,\n"
                "  emitter = %s\n"
                "]",
                m_position.toString(),
                m_radius,
                m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
                m_emitter ? indent(m_emitter->toString()) : std::string("null"));
    }

protected:
    Point3f m_position;
    float m_radius;
};

NORI_REGISTER_CLASS(Sphere, "sphere");
NORI_NAMESPACE_END
