#include "AnalyticModels.h"
#include <stdexcept>
#include <Eigen/Dense>

namespace analytical {

    double alphaIsotropicVoigt(double alpha_m, double alpha_f, double phi_f)
    {
        double phi_m = 1.0 - phi_f;

        return (alpha_f * phi_f + alpha_m * phi_m);
    }

    double alphaIsotropicTurner(double alpha_m, double alpha_f, double E_m, double nu_m, double E_f, double nu_f, double phi_f)
    {
        double phi_m = 1.0 - phi_f;
        double K_m = E_m / (3.0 * (1 - 2.0 * nu_m));
        double K_f = E_f / (3.0 * (1 - 2.0 * nu_f));

        return (alpha_m * K_m * phi_m + alpha_f * K_f * phi_f) / (K_m * phi_m + K_f * phi_f);
    }

    double shearModulus(double E, double nu)
    {
        return E / (2.0 * (1.0 + nu));
    }

    double bulkModulus(double E, double nu)
    {
        return E / (3.0 * (1.0 - 2.0 * nu));
    }

    double alphaIsotropicKerner(double alpha_m, double alpha_f, double E_m, double nu_m, double E_f, double nu_f, double V_f)
    {
        double V_m = 1.0 - V_f;

        double K_m = bulkModulus(E_m, nu_m);
        double K_f = bulkModulus(E_f, nu_f);
        double G_m = shearModulus(E_m, nu_m);

        double denom = V_m * K_m + V_f * K_f + (3.0 * K_m * K_f) / (4.0 * G_m);

        return alpha_m * V_m + alpha_f * V_f + V_m * V_f * (alpha_f - alpha_m) * (K_f - K_m) / denom;
    }

    double hashinShtrikmanBulkLower(double K_m, double G_m, double K_f, double phi_f)
    {
        if (phi_f < 0.0 || phi_f > 1.0)
            throw std::invalid_argument("hashinShtrikmanBulkLower: phi_f must be in [0, 1]");

        if (K_m <= 0.0 || G_m <= 0.0 || K_f <= 0.0)
            throw std::invalid_argument("hashinShtrikmanBulkLower: moduli must be > 0");

        if (std::abs(K_f - K_m) < 1e-30)
            return K_m;

        const double denom = 1.0 / (K_f - K_m) + 3.0 * (1.0 - phi_f) / (3.0 * K_m + 4.0 * G_m);

        return K_m + phi_f / denom;
    }

    double hashinShtrikmanBulkUpper(double K_m, double G_m, double K_f, double G_f, double phi_f)
    {
        if (phi_f < 0.0 || phi_f > 1.0)
            throw std::invalid_argument("hashinShtrikmanBulkUpper: phi_f must be in [0, 1]");

        if (K_m <= 0.0 || G_m <= 0.0 || K_f <= 0.0 || G_f <= 0.0)
            throw std::invalid_argument("hashinShtrikmanBulkUpper: moduli must be > 0");

        if (std::abs(K_f - K_m) < 1e-30)
            return K_m;

        const double denom = 1.0 / (K_m - K_f) + 3.0 * phi_f / (3.0 * K_f + 4.0 * G_f);

        return K_f + (1.0 - phi_f) / denom;
    }

    double levinAlphaFromBulk(double alpha_m, double alpha_f, double K_m, double K_f, double K_eff)
    {
        if (K_m <= 0.0 || K_f <= 0.0 || K_eff <= 0.0)
            throw std::invalid_argument("levinAlphaFromBulk: bulk moduli must be > 0");

        if (std::abs(K_f - K_m) < 1e-30)
            return alpha_m;

        return alpha_m + (K_f / K_eff) * ((K_m - K_eff) / (K_m - K_f)) * (alpha_f - alpha_m);
    }

    HashinShtrikmanCTEBounds hashinShtrikmanAlphaBounds(double E_m, double nu_m, double alpha_m, double E_f, double nu_f, double alpha_f, double phi_f)
    {
        if (phi_f < 0.0 || phi_f > 1.0)
            throw std::invalid_argument("hashinShtrikmanAlphaBounds: phi_f must be in [0, 1]");

        const double K_m = bulkModulus(E_m, nu_m);
        const double G_m = shearModulus(E_m, nu_m);
        const double K_f = bulkModulus(E_f, nu_f);
        const double G_f = shearModulus(E_f, nu_f);

        const double K_hs_lower = hashinShtrikmanBulkLower(K_m, G_m, K_f, phi_f);
        const double K_hs_upper = hashinShtrikmanBulkUpper(K_m, G_m, K_f, G_f, phi_f);

        const double alpha_lower = levinAlphaFromBulk(alpha_m, alpha_f, K_m, K_f, K_hs_upper);
        const double alpha_upper = levinAlphaFromBulk(alpha_m, alpha_f, K_m, K_f, K_hs_lower);

        return { alpha_lower, alpha_upper };
    }

    using M6 = Eigen::Matrix<double, 6, 6>;
    using V6 = Eigen::Matrix<double, 6, 1>;

    // isotropic stiffness in Kelvin/Mandel notation
    M6 isotropicStiffnessKelvin(double E, double nu)
    {
        const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double mu     = E / (2.0 * (1.0 + nu));

        M6 C = M6::Zero();
        C(0,0) = C(1,1) = C(2,2) = lambda + 2.0 * mu;
        C(0,1) = C(0,2) = C(1,0) = C(1,2) = C(2,0) = C(2,1) = lambda;

        // Kelvin/Mandel shear entries
        C(3,3) = C(4,4) = C(5,5) = 2.0 * mu;
        return C;
    }

    V6 isotropicAlphaKelvin(double alpha)
    {
        V6 a = V6::Zero();
        a(0) = a(1) = a(2) = alpha;
        return a;
    }

    // Eshelby tensor for a PROLATE spheroid, axis of symmetry along z.
    // aspect = c/a = rz/rx, with rx = ry = a, rz = c, and aspect >= 1.
    // Returned in Kelvin/Mandel 6x6 notation.
    M6 eshelbySpheroidKelvinAxisZ(double nu_m, double aspect)
    {
        const double eps = 1e-12;
        M6 S = M6::Zero();

        // sphere limit
        if (std::abs(aspect - 1.0) < eps)
        {
            const double den = 15.0 * (1.0 - nu_m);
            const double S1111 = (7.0 - 5.0 * nu_m) / den;
            const double S1122 = (5.0 * nu_m - 1.0) / den;
            const double S1212 = (4.0 - 5.0 * nu_m) / den;

            S(0,0) = S(1,1) = S(2,2) = S1111;
            S(0,1) = S(0,2) = S(1,0) = S(1,2) = S(2,0) = S(2,1) = S1122;
            S(3,3) = S(4,4) = S(5,5) = 2.0 * S1212;
            return S;
        }

        // prolate branch: aspect > 1
        const double a  = aspect;
        const double a2 = a * a;
        const double root = std::sqrt(a2 - 1.0);
        const double g = a / std::pow(a2 - 1.0, 1.5) * (a * root - std::acosh(a));

        // These are the standard closed-form spheroidal Eshelby components
        // in regular tensor notation, originally written for symmetry axis along x.
        // We remap them here to symmetry axis along z.
        const double A = 1.0 / (2.0 * (1.0 - nu_m)) *
                        (1.0 - 2.0 * nu_m + (3.0 * a2 - 1.0) / (a2 - 1.0)
                        - (1.0 - 2.0 * nu_m + 3.0 * a2 / (a2 - 1.0)) * g);

        const double B = 3.0 / (8.0 * (1.0 - nu_m)) * a2 / (a2 - 1.0)
                        + 1.0 / (4.0 * (1.0 - nu_m))
                        * (1.0 - 2.0 * nu_m - 9.0 / (4.0 * (a2 - 1.0))) * g;

        const double C = 1.0 / (4.0 * (1.0 - nu_m))
                        * (a2 / (2.0 * (a2 - 1.0))
                            - (1.0 - 2.0 * nu_m + 3.0 / (4.0 * (a2 - 1.0))) * g);

        const double D = -1.0 / (2.0 * (1.0 - nu_m)) * a2 / (a2 - 1.0)
                        + 1.0 / (4.0 * (1.0 - nu_m))
                        * (3.0 * a2 / (a2 - 1.0) - (1.0 - 2.0 * nu_m)) * g;

        const double E = -1.0 / (2.0 * (1.0 - nu_m))
                        * (1.0 - 2.0 * nu_m + 1.0 / (a2 - 1.0))
                        + 1.0 / (2.0 * (1.0 - nu_m))
                        * (1.0 - 2.0 * nu_m + 3.0 / (2.0 * (a2 - 1.0))) * g;

        const double F = 1.0 / (4.0 * (1.0 - nu_m))
                        * (a2 / (2.0 * (a2 - 1.0))
                            + (1.0 - 2.0 * nu_m - 3.0 / (4.0 * (a2 - 1.0))) * g);

        const double G = 1.0 / (4.0 * (1.0 - nu_m))
                        * (1.0 - 2.0 * nu_m - (a2 + 1.0) / (a2 - 1.0)
                            - 0.5 * (1.0 - 2.0 * nu_m
                                    - 3.0 * (a2 + 1.0) / (a2 - 1.0)) * g);

        // Remapped to symmetry axis along z
        // normal-normal block
        S(0,0) = B;  S(1,1) = B;  S(2,2) = A;
        S(0,1) = C;  S(1,0) = C;

        // first pair transverse, second pair axial
        S(0,2) = D;  S(1,2) = D;

        // first pair axial, second pair transverse
        S(2,0) = E;  S(2,1) = E;

        // shear block in Kelvin/Mandel has factor 2 relative to tensor shear components
        S(3,3) = 2.0 * G; // yz
        S(4,4) = 2.0 * G; // xz
        S(5,5) = 2.0 * F; // xy

        return S;
    }

    MTResult moriTanakaCTE_spheroid_axisZ(
        double E_m, double nu_m, double alpha_m,
        double E_f, double nu_f, double alpha_f,
        double vf, double aspect // aspect = rz/rx, rx=ry
    )
    {
        const double vm = 1.0 - vf;

        const M6 I  = M6::Identity();
        const M6 Cm = isotropicStiffnessKelvin(E_m, nu_m);
        const M6 Cf = isotropicStiffnessKelvin(E_f, nu_f);
        const M6 Dm = Cm.inverse();
        const M6 Df = Cf.inverse();

        const V6 am = isotropicAlphaKelvin(alpha_m);
        const V6 af = isotropicAlphaKelvin(alpha_f);

        const M6 S = eshelbySpheroidKelvinAxisZ(nu_m, aspect);

        // elastic MT
        const M6 Adil = (I + S * Dm * (Cf - Cm)).inverse();
        const M6 Amt  = Adil * (vm * I + vf * Adil).inverse();
        const M6 Ceff = Cm + vf * (Cf - Cm) * Amt;

        // thermal-expansion MT
        const M6 Bdil = (I + Cm * (I - S) * (Df - Dm)).inverse();
        const M6 Bmt  = Bdil * (vm * I + vf * Bdil).inverse();

        const V6 aeff = am + vf * (Bmt * (af - am));

        MTResult out;
        out.alpha_x = aeff(0);
        out.alpha_y = aeff(1);
        out.alpha_z = aeff(2);
        out.Ceff    = Ceff;
        return out;
    }
};
