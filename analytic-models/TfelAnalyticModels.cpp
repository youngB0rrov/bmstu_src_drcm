#include "TfelAnalyticModels.h"
#include "TFEL/Material/IsotropicEshelbyTensor.hxx"
#include "TFEL/Material/LocalisationTensor.hxx"
#include "TFEL/Material/LinearHomogenizationSchemes.hxx"
#include "TFEL/Material/IsotropicModuli.hxx"
#include "TFEL/Math/tvector.hxx"
#include "AnalyticModels.h"

namespace analytical
{
    MTResult moriTanakaCTE_ellipsoid_tfel(double E_m, double nu_m, double alpha_m, double E_f, double nu_f, double alpha_f, double vf, double rx, double ry, double rz)
    {
        using M6 = Eigen::Matrix<double, 6, 6>;
        using V6 = Eigen::Matrix<double, 6, 1>;
        using TVec3 = tfel::math::tvector<3u, double>;

        const double vm = 1.0 - vf;

        const M6 Cm = analytical::isotropicStiffnessKelvin(E_m, nu_m);
        const M6 Cf = analytical::isotropicStiffnessKelvin(E_f, nu_f);

        const V6 am = isotropicAlphaKelvin(alpha_m);
        const V6 af = isotropicAlphaKelvin(alpha_f);

        const auto IM0 = tfel::material::YoungNuModuli<double>(E_m, nu_m);
        const auto IMi = tfel::material::YoungNuModuli<double>(E_f, nu_f);

        TVec3 n_a;
        n_a[0] = 1.0; n_a[1] = 0.0; n_a[2] = 0.0;

        TVec3 n_b;
        n_b[0] = 0.0; n_b[1] = 1.0; n_b[2] = 0.0;

        // Dilute strain localisation tensor for the inclusion
        const auto A_dil_tfel =
            tfel::material::homogenization::elasticity
                ::computeEllipsoidLocalisationTensor<double>(IM0, IMi, n_a, rx, n_b, ry, rz);

        M6 A_dil = tfelToEigen6x6(A_dil_tfel);

        // Effective stiffness from TFEL MT
        const auto Ceff_tfel =
            tfel::material::homogenization::elasticity
                ::computeOrientedMoriTanakaScheme<double>(IM0, vf, IMi, n_a, rx, n_b, ry, rz);

        M6 Ceff = tfelToEigen6x6(Ceff_tfel);

        // Mechanical strain concentration tensors in MT
        const M6 A_m = (vm * M6::Identity() + vf * A_dil).inverse();
        const M6 A_f = A_dil * A_m;

        // beta = C * alpha
        const V6 beta_m = Cm * am;
        const V6 beta_f = Cf * af;

        // Mandel–Levin / Benveniste formula:
        // alpha_eff = Ceff^{-1} * [ vm * A_m^T * beta_m + vf * A_f^T * beta_f ]
        const V6 beta_eff =
            vm * A_m.transpose() * beta_m +
            vf * A_f.transpose() * beta_f;

        const V6 aeff = Ceff.fullPivLu().solve(beta_eff);

        MTResult out{};
        out.alpha_x = aeff(0);
        out.alpha_y = aeff(1);
        out.alpha_z = aeff(2);
        out.Ceff    = Ceff;
        return out;
        // using M6 = Eigen::Matrix<double, 6, 6>;
        // using V6 = Eigen::Matrix<double, 6, 1>;
        // using TVec3 = tfel::math::tvector<3u, double>;

        // const double vm = 1.0 - vf;
        // const M6 I = M6::Identity();

        // const M6 Cm = analytical::isotropicStiffnessKelvin(E_m, nu_m);
        // const M6 Cf = analytical::isotropicStiffnessKelvin(E_f, nu_f);

        // const V6 am = isotropicAlphaKelvin(alpha_m);
        // const V6 af = isotropicAlphaKelvin(alpha_f);

        // const auto IM0 = tfel::material::YoungNuModuli<double>(E_m, nu_m);
        // const auto IMi = tfel::material::YoungNuModuli<double>(E_f, nu_f);

        // TVec3 n_a;
        // n_a[0] = 1.0; n_a[1] = 0.0; n_a[2] = 0.0;

        // TVec3 n_b;
        // n_b[0] = 0.0; n_b[1] = 1.0; n_b[2] = 0.0;

        // // dilute localisation tensor
        // const auto A_dil_tfel =
        //     tfel::material::homogenization::elasticity
        //         ::computeEllipsoidLocalisationTensor<double>(IM0, IMi, n_a, rx, n_b, ry, rz);

        // const M6 Adil = tfelToEigen6x6(A_dil_tfel);

        // const auto Ceff_tfel =
        //     tfel::material::homogenization::elasticity
        //         ::computeOrientedMoriTanakaScheme<double>(IM0, vf, IMi, n_a, rx, n_b, ry, rz);

        // const M6 Ceff = tfelToEigen6x6(Ceff_tfel);

        // // MT-Benveniste: сначала beta = C * alpha
        // const V6 beta_m = Cm * am;
        // const V6 beta_f = Cf * af;

        // const M6 Abar = vm * I + vf * Adil;

        // // важная часть: считаем effective (C:alpha), а не alpha напрямую
        // const V6 beta_eff =
        //     Abar.inverse().transpose() *
        //     (vm * beta_m + vf * (Cf * Adil * af));

        // const V6 aeff = Ceff.inverse() * beta_eff;

        // MTResult out{};
        // out.alpha_x = aeff(0);
        // out.alpha_y = aeff(1);
        // out.alpha_z = aeff(2);
        // out.Ceff = Ceff;
        // return out;
    }
}