#include "TFEL/Material/IsotropicEshelbyTensor.hxx"
#include "TFEL/Material/LocalisationTensor.hxx"
#include "TFEL/Material/LinearHomogenizationSchemes.hxx"
#include "TFEL/Material/IsotropicModuli.hxx"
#include "TFEL/Math/tvector.hxx"
#include <Eigen/Dense>
#include "AnalyticModels.h"

namespace analytical
{
    using M6 = Eigen::Matrix<double, 6, 6>;
    using V6 = Eigen::Matrix<double, 6, 1>;
    using TVec3 = tfel::math::tvector<3u, double>;

    template <typename T4>
    M6 tfelToEigen6x6(const T4& A)
    {
        M6 out = M6::Zero();
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                out(i, j) = static_cast<double>(A(i, j));
            }
        }
        return out;
    }

    MTResult moriTanakaCTE_ellipsoid_tfel(
        double E_m, double nu_m, double alpha_m,
        double E_f, double nu_f, double alpha_f,
        double vf,
        double rx, double ry, double rz);
}