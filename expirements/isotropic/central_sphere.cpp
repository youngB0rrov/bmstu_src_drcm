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

#define _USE_MATH_DEFINES

using namespace nla3d;
using namespace hom3d;

namespace expirements::sphere_center
{
    class Problem: public hom3d::ThermoElasticHomProblem
    {
        virtual void mesh()
        {
            std::string geo_text = "algebraic3d \
                    solid box1= plane(0,0,0;0,0,-1)-bc=16; \
                    solid box2= plane(0,0,0;0,-1,0)-bc=14; \
                    solid box3= plane(0,0,0;-1,0,0)-bc=18; \
                    solid box4= plane(1,1,1;1,0,0)-bc=22; \
                    solid box5= plane(1,1,1;0,1,0)-bc=15; \
                    solid box6= plane(1,1,1;0,0,1)-bc=19; \
                    solid YAP = box1 and box2 and box3 and box4 and box5 and box6; \
                    solid sphere0 = sphere(0.5,0.5,0.5; radius)-bc=100; \
                    solid matrix = YAP and not sphere0; \
                    solid arm = sphere0 and YAP; \
                    tlo matrix -col=[0,0,1] -transparent; \
                    tlo arm  -col=[1,1,0] ; \
                    ";
            std::string subs = "radius";
            size_t start_pos = geo_text.find(subs);
            geo_text.replace(start_pos, subs.length(), std::to_string(_radius));
            std::ofstream f("disp.geo");
            f << geo_text;
            f.close();
            system("netgen -geofile=disp.geo -meshfile=disp.neu -meshfiletype=\"Neutral Format\" -batchmode -veryfine >mesh.out");
            fromNeuFile("disp.neu");
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
        double _radius;
        float _a, _b, _c;

    public:
        Problem(double radius = 0.1, float a = 1.0): _radius(radius), _a(a), _b(a), _c(a) {};
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

        std::ofstream data("cte_vs_vf.dat");

        if (!data)
        {
            std::cerr << "Не удалось создать cte_vs_vf.dat\n";
            return 1;
        }

        data << "vf alpha_x alpha_y alpha_z alpha_voigt alpha_turner alpha_kerner alpha_hs_lower alpha_hs_upper\n";
        data << std::scientific << std::setprecision(10);

        for (double r = 0.05; r < RVE_a / 2.0; r+=0.05)
        {
            double phi_f = geom::sphereVolumeFraction(r, RVE_a);
            double fiberVolumePercent = phi_f * 100.0;

            LOG(INFO) << "Fiber vol.%: " << fiberVolumePercent << "%" << std::endl;

            Problem problem(r, RVE_a);
            problem.setHomMethod(AHM_SYM_CONST);
            problem.setName("ThermoElastic");
            problem.init();
            problem.solve();
            problem.printResults();

            double alpha_voigt = analytical::alphaIsotropicVoigt(alpha_m, alpha_f, phi_f);
            double alpha_turner = analytical::alphaIsotropicTurner(alpha_m, alpha_f, E_m, nu_m, E_f, nu_f, phi_f);
            double alpha_kerner = analytical::alphaIsotropicKerner(alpha_m, alpha_f, E_m, nu_m, E_f, nu_f, phi_f);
            auto hs = analytical::hashinShtrikmanAlphaBounds(E_m, nu_m, alpha_m, E_f, nu_f, alpha_f, phi_f);

            double alpha_hs_lower = hs.alpha_lower;
            double alpha_hs_upper = hs.alpha_upper;

            data << fiberVolumePercent << ' ' << problem.alphaX() * 1e6 << ' ' << problem.alphaY() * 1e6 << ' ' << problem.alphaZ() * 1e6 << ' ' << alpha_voigt * 1e6 << ' ' << alpha_turner * 1e6 << ' ' << alpha_kerner * 1e6 << ' ' << alpha_hs_lower * 1e6 << ' ' << alpha_hs_upper * 1e6 << '\n';
        }

        data.close();

        std::vector<plot::PlotSpec> plots =
        {
            {"results/vf_vs_cte.png", "Зависимость КЛТР от объемной доли включений", plot::PlotType::Basic, {}},
            {"results/vf_vs_cte_voigt.png", "Сравнение численного решения с оценкой Фойгта", plot::PlotType::Comparison, 
            {
                {2, "alpha_x", "linespoints", 2, 7},
                {3, "alpha_y", "linespoints", 2, 5},
                {4, "alpha_z", "linespoints", 2, 9},
                {5, "alpha_v",   "linespoints", 2, 11},
            }},
            {"results/vf_vs_cte_turner.png", "Сравнение численного решения с оценкой Тёрнера", plot::PlotType::Comparison,
            {
                {2, "alpha_x", "linespoints", 2, 7},
                {3, "alpha_y", "linespoints", 2, 5},
                {4, "alpha_z", "linespoints", 2, 9},
                {6, "alpha_t",   "linespoints", 2, 11},
            }},
            {"results/vf_vs_cte_kerner.png", "Сравнение численного решения с оценкой Кернера", plot::PlotType::Comparison,
            {
                {2, "alpha_x", "linespoints", 2, 7},
                {3, "alpha_y", "linespoints", 2, 5},
                {4, "alpha_z", "linespoints", 2, 9},
                {7, "alpha_k",   "linespoints", 2, 11},
            }},
            {"results/vf_vs_cte_hashin_strikman.png", "Сравнение численного решения с оценкой Хашина-Штрикмана", plot::PlotType::Comparison,
            {
                {2, "alpha_x", "linespoints", 2, 7},
                {3, "alpha_y", "linespoints", 2, 5},
                {4, "alpha_z", "linespoints", 2, 9},
                {8, "alpha_{lower}",   "linespoints", 2, 11},
                {9, "alpha_{upper}",   "linespoints", 2, 13},
            }}
        };

        if (!plot::generatePlots("cte_vs_vf.dat", plots))
        {
            return 1;
        }

        LOG(INFO) << "Графики сохранены в директории results\n";

        return 0;
    }
}