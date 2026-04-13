#include "ThermoElasticHomProblem.h"
#include <fstream>
#include "FEReaders.h"
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include "GeometryUtils.h"
#include "AnalyticModels.h"
#include "PlotUtils.h"
#include <vector>
#include <algorithm>

#define _USE_MATH_DEFINES

using namespace nla3d;
using namespace hom3d;

namespace expirements::spheroid_center
{
    class Problem: public hom3d::ThermoElasticHomProblem
    {
        virtual void mesh() override
        {
            std::ostringstream geo;
            geo << "algebraic3d\n"
                << "solid box1= plane(0,0,0;0,0,-1)-bc=16;\n"
                << "solid box2= plane(0,0,0;0,-1,0)-bc=14;\n"
                << "solid box3= plane(0,0,0;-1,0,0)-bc=18;\n"
                << "solid box4= plane(1,1,1;1,0,0)-bc=22;\n"
                << "solid box5= plane(1,1,1;0,1,0)-bc=15;\n"
                << "solid box6= plane(1,1,1;0,0,1)-bc=19;\n"
                << "solid YAP = box1 and box2 and box3 and box4 and box5 and box6;\n"
                << "solid spheroid0 = ellipsoid(0.5,0.5,0.5; " << _rx << ",0,0; 0," << _ry << ",0; 0,0," << _rz << ")-bc=100;\n"
                << "solid matrix = YAP and not spheroid0;\n"
                << "solid arm = spheroid0 and YAP;\n"
                << "tlo matrix -col=[0,0,1] -transparent;\n"
                << "tlo arm -col=[1,1,0];\n";

            std::ofstream f("disp_spheroid.geo");
            f << geo.str();
            f.close();

            system("netgen -geofile=disp_spheroid.geo -meshfile=disp_spheroid.neu -meshfiletype=\"Neutral Format\" -batchmode -veryfine >mesh_spheroid.out");
            fromNeuFile("disp_spheroid.neu");
            setRVEsize(_a, _b, _c);
        }

        virtual void phases()
        {
            // Фаза матрицы (эпоксидная смола)
            addElasticOrthotropicPhase(3.07e9, 3.07e9, 3.07e9, 0.35, 0.35, 0.35, 63e-6, 63e-6, 63e-6);

            // Фаза включения (оксид кремния)
            addElasticOrthotropicPhase(76e9, 76e9, 76e9, 0.17, 0.17, 0.17, 0.5e-6, 0.5e-6, 0.5e-6);
        }

    private:
        double _rx, _ry, _rz;
        float _a, _b, _c;

    public:
        Problem(double rx = 0.1, double ry = 0.1, double rz = 0.1, float a = 1.0): _rx(rx), _ry(ry), _rz(rz), _a(a), _b(a), _c(a) {}
    };

    int run()
    {
        double E_m = 3.07e9;
        double nu_m = 0.35;
        double alpha_m = 63e-6;

        double E_f = 76e9;
        double nu_f = 0.17;
        double alpha_f = 0.5e-6;

        double RVE_a = 1.0;
        double ratio_aspect = 1.3;
        double ratio_aspect_y = 1.0;
        
        std::ofstream data("cte_vs_vf_spheroid.dat");

        if (!data)
        {
            std::cerr << "Не удалось создать cte_vs_vf_spheroid.dat\n";
            return 1;
        }

        data << "vf alpha_x alpha_y alpha_z alpha_ort alpha_parallel alpha_mean delta_alpha anis_ratio alpha_x_mt alpha_y_mt alpha_z_mt\n";
        data << std::scientific << std::setprecision(10);

        for (double rxy = 0.1; rxy < 0.5 / std::max(1.0, ratio_aspect); rxy+=0.05)
        {
            double rx = rxy;
            double ry = rxy * ratio_aspect_y;
            double rz = ratio_aspect * rxy;

            double phi_f = geom::ellipsoidVolumeFraction(rx, ry, rz, RVE_a);
            double fiberVolumePercent = phi_f * 100.0;

            LOG(INFO) << "Fiber vol.%: " << fiberVolumePercent << "%" << std::endl;

            Problem problem(rx, ry, rz, RVE_a);
            problem.setHomMethod(AHM_SYM_CONST);
            problem.setName("ThermoElasticSpheroid");
            problem.init();
            problem.solve();

            // Анизотропия
            double alpha_perp = 0.5 * (problem.alphaX() + problem.alphaY());
            double alpha_parallel = problem.alphaZ();
            double alpha_mean = (2.0 * alpha_perp + alpha_parallel) / 3.0;
            double delta_alpha = alpha_perp - alpha_parallel;
            double anis_ratio = delta_alpha / alpha_mean;

            // Мори-Танака
            analytical::MTResult mtResult = analytical::moriTanakaCTE_spheroid_axisZ(E_m, nu_m, alpha_m, E_f, nu_f, alpha_f, phi_f, ratio_aspect);

            double ax_mt = mtResult.alpha_x;
            double ay_mt = mtResult.alpha_y;
            double az_mt = mtResult.alpha_z;

            data << fiberVolumePercent << ' '
                << problem.alphaX() * 1e6 << ' '
                << problem.alphaY() * 1e6 << ' '
                << problem.alphaZ() * 1e6 << ' '
                << alpha_perp * 1e6 << ' '
                << alpha_parallel * 1e6 << ' '
                << alpha_mean * 1e6 << ' '
                << delta_alpha * 1e6 << ' '
                << anis_ratio << ' '
                << ax_mt * 1e6 << ' '
                << ay_mt * 1e6 << ' '
                << az_mt * 1e6 << ' '
                << '\n';
        }

        data.close();

        std::vector<plot::PlotSpec> plots =
        {
            {"results/vf_vs_cte_spheroid.png", "Зависимость КЛТР от объемной доли включений для сфероида", plot::PlotType::Basic, {}},
            {"results/vf_vs_cte_spheroid_anisotropy.png", "Зависимость абсолютной анизотропии от объемной доли включений для сфероида", plot::PlotType::Comparison,
            {
                {8, "delta_{alpha}", "linespoints", 2, 7}
            }, "Объемная доля включений, %", "Абсолютная анизотропия, 10^{-6}/K"},
            {"results/vf_vs_cte_spheroid_mori_tanaka.png", "Сравнение численного решения с моделью Мори-Танака", plot::PlotType::Comparison,
            {
                {2, "alpha_x", "linespoints", 2, 5},
                {3, "alpha_y", "linespoints", 2, 7},
                {4, "alpha_z", "linespoints", 2, 9},
                {10, "alpha_x (MT)", "linespoints", 2, 11, 2},
                {11, "alpha_y (MT)", "linespoints", 2, 13, 2},
                {12, "alpha_z (MT)", "linespoints", 2, 3, 2}
            }},
        };

        if (!plot::generatePlots("cte_vs_vf_spheroid.dat", plots))
        {
            return 1;
        }

        LOG(INFO) << "Графики сохранены в директории results\n";

        return 0;
    }
}