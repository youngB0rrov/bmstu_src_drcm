#pragma once

namespace geom 
{
    double sphereVolume(double r);
    double ellipsoidVolume(double rx, double ry, double rz);
    double cubeVolume(double a);
    double sphereVolumeFraction(double r, double a);
    double ellipsoidVolumeFraction(double rx, double ry, double rz, double a);
};