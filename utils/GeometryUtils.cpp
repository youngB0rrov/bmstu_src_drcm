#include "GeometryUtils.h"
#include <cmath>

#define _USE_MATH_DEFINES

namespace geom
{
    double sphereVolume(double r)
    {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    double ellipsoidVolume(double rx, double ry, double rz)
    {
        return (4.0 / 3.0) * M_PI * rx * ry * rz;
    }

    double cubeVolume(double a)
    {
        return a * a * a;
    }

    double sphereVolumeFraction(double r, double a)
    {
        return sphereVolume(r) / cubeVolume(a);
    }
    double ellipsoidVolumeFraction(double rx, double ry, double rz, double a)
    {
        return ellipsoidVolume(rx, ry, rz) / cubeVolume(a);
    }
};