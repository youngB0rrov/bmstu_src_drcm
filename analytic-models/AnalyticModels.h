#pragma once
#include <Eigen/Dense>

using M6 = Eigen::Matrix<double, 6, 6>;
using V6 = Eigen::Matrix<double, 6, 1>;

namespace analytical
{
    struct HashinShtrikmanCTEBounds
    {
        double alpha_lower; // нижняя граница эффективного КЛТР
        double alpha_upper; // верхняя граница эффективного КЛТР
    };

    struct MTResult
    {
        double alpha_x;
        double alpha_y;
        double alpha_z;
        Eigen::Matrix<double, 6, 6> Ceff;
    };

    double alphaIsotropicVoigt(double alpha_m, double alpha_f, double phi_f);
    double alphaIsotropicTurner(double alpha_m, double alpha_f, double E_m, double nu_m, double E_f, double nu_f, double phi_f);
    double shearModulus(double E, double nu);
    double bulkModulus(double E, double nu);
    double alphaIsotropicKerner(double alpha_m, double alpha_f, double E_m, double nu_m, double E_f, double nu_f, double V_f);
    double hashinShtrikmanBulkLower(double K_m, double G_m, double K_f, double phi_f);
    double hashinShtrikmanBulkUpper(double K_m, double G_m, double K_f, double G_f, double phi_f);
    double levinAlphaFromBulk(double alpha_m, double alpha_f, double K_m, double K_f, double K_eff);
    HashinShtrikmanCTEBounds hashinShtrikmanAlphaBounds(double E_m, double nu_m, double alpha_m, double E_f, double nu_f, double alpha_f, double phi_f);
    MTResult moriTanakaCTE_spheroid_axisZ(double E_m, double nu_m, double alpha_m, double E_f, double nu_f, double alpha_f, double vf, double aspect);
    M6 isotropicStiffnessKelvin(double E, double nu);
    V6 isotropicAlphaKelvin(double alpha);
};