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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    double r = std::sqrt(sample.x()), theta = 2 * M_PI * sample.y();
    return Point2f(r * cos(theta), r * sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if (p.norm() < 1.0) {
        return 1.0 / M_PI;
    }
    else {
        return 0.0;
    }
}

Vector3f Warp::squareToUniformCylinder(const Point2f& sample) {
    double fai = 2.0 * M_PI * sample.y();
    return Vector3f(std::cos(fai), std::sin(fai), 2 * sample.x() - 1);
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    Vector3f Cylinder_res = squareToUniformCylinder(sample);
    double z = Cylinder_res.z() * (1.0 - cosThetaMax) / 2.0 + (1.0 + cosThetaMax) / 2.0;
    double r = std::sqrt(1 - z * z);
    return Vector3f(r * Cylinder_res.x(), r * Cylinder_res.y(), z);
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    if (v.z() > cosThetaMax) {
        return 1.0 / 2.0 / M_PI / (1.0 - cosThetaMax);
    }
    else {
        return 0.0;
    }
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    Vector3f Cylinder_res = squareToUniformCylinder(sample);
    double r = std::sqrt(1 - Cylinder_res.z() * Cylinder_res.z());
    return Vector3f(r * Cylinder_res.x(), r * Cylinder_res.y(), Cylinder_res.z());
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return 1.0 / 4.0 / M_PI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    Vector3f Cylinder_res = squareToUniformCylinder(sample);
    double z = Cylinder_res.z() * 0.5 + 0.5;
    double r = std::sqrt(1 - z * z);
    return Vector3f(r * Cylinder_res.x(), r * Cylinder_res.y(), z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    if (v.z() > 0) {
        return 1.0 / 2.0 / M_PI;
    }
    else {
        return 0.0;
    }
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Vector3f Cylinder_res = squareToUniformCylinder(sample);
    double cos_theta = std::sqrt(0.5 * (1.0 + Cylinder_res.z()));
    double sin_theta = std::sqrt(0.5 * (1.0 - Cylinder_res.z()));
    return Vector3f(sin_theta * Cylinder_res.x(), sin_theta * Cylinder_res.y(), cos_theta);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    //cout << "PDF: " <<  v.x() << " , " << v.y() << " , " << v.z() << endl;
    if (v.z() > 0.0) {
        return v.z() / M_PI;
    }
    else {
        return 0.0;
    }
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    Vector3f Cylinder_res = squareToUniformCylinder(sample);
    //double theta = std::atan(alpha * std::sqrt(-std::log((1.0 - Cylinder_res.z()) * 0.5)));
    double tan_theta = alpha * std::sqrt(-std::log((1.0 - Cylinder_res.z()) * 0.5));
    double cos_theta = (1.0 / std::sqrt(1 + tan_theta * tan_theta)) * tan_theta / std::abs(tan_theta);
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    //return Vector3f(std::sin(theta) * Cylinder_res.x(), std::sin(theta) * Cylinder_res.y(), std::cos(theta));
    return Vector3f(sin_theta * Cylinder_res.x(), sin_theta * Cylinder_res.y(), cos_theta);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (m.z() > 0.0) {
        return std::exp(-(1.0 / m.z() / m.z() - 1) / alpha / alpha) / pow(m.z(), 3) / alpha / alpha / M_PI;
    }
    else{
        return 0.0;
    }
}

Vector3f Warp::squareToUniformTriangle(const Point2f &sample) {
    float su1 = sqrtf(sample.x());
    float u = 1.f - su1, v = sample.y() * su1;
    return Vector3f(u,v,1.f-u-v);
}

NORI_NAMESPACE_END
