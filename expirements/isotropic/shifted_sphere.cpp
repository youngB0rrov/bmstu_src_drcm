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

namespace expirements::sphere_shifted
{
    static void replaceAll(std::string& text, const std::string& from, const std::string& to)
    {
        if (from.empty()) return;

        size_t start_pos = 0;
        while ((start_pos = text.find(from, start_pos)) != std::string::npos)
        {
            text.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }
    }

    static std::string toStringPrec(double value)
    {
        std::ostringstream ss;
        ss << std::setprecision(16) << value;
        return ss.str();
    }

    class Problem : public hom3d::ThermoElasticHomProblem
    {
        virtual void mesh() override
        {
            std::string geo_text = R"(
                algebraic3d
                solid box1 = plane(0,0,0;0,0,-1)-bc=16;
                solid box2 = plane(0,0,0;0,-1,0)-bc=14;
                solid box3 = plane(0,0,0;-1,0,0)-bc=18;
                solid box4 = plane(1,1,1;1,0,0)-bc=22;
                solid box5 = plane(1,1,1;0,1,0)-bc=15;
                solid box6 = plane(1,1,1;0,0,1)-bc=19;

                solid YAP = box1 and box2 and box3 and box4 and box5 and box6;
                solid sphere0 = sphere(CX,CY,CZ; RADIUS)-bc=20;

                solid matrix = YAP and not sphere0;
                solid arm = sphere0 and YAP;

                tlo matrix -col=[0,0,1] -transparent;
                tlo arm -col=[1,1,0];
                )";

            replaceAll(geo_text, "CX", toStringPrec(_cx));
            replaceAll(geo_text, "CY", toStringPrec(_cy));
            replaceAll(geo_text, "CZ", toStringPrec(_cz));
            replaceAll(geo_text, "RADIUS", toStringPrec(_radius));

            std::ofstream f("disp_shifted.geo");
            f << geo_text;
            f.close();

            system("netgen -geofile=disp_shifted.geo -meshfile=disp_shifted.neu -meshfiletype=\"Neutral Format\" -batchmode -veryfine > shifted_mesh.out");

            fromNeuFile("disp_shifted.neu");
            setRVEsize(_a, _b, _c);
        }

        virtual void phases() override
        {
            // Фаза матрицы (эпоксидная смола)
            addElasticOrthotropicPhase(
                3.07e9, 3.07e9, 3.07e9,
                0.35, 0.35, 0.35,
                63e-6, 63e-6, 63e-6
            );

            // Фаза включения (оксид кремния)
            addElasticOrthotropicPhase(
                76e9, 76e9, 76e9,
                0.17, 0.17, 0.17,
                0.5e-6, 0.5e-6, 0.5e-6
            );
        }

    private:
        double _radius;
        double _cx, _cy, _cz;
        float _a, _b, _c;

    public:
        Problem(double radius, double cx, double cy, double cz, float a = 1.0)
            : _radius(radius), _cx(cx), _cy(cy), _cz(cz), _a(a), _b(a), _c(a)
        {
            if (_cx < _radius || _cx > _a - _radius ||
                _cy < _radius || _cy > _b - _radius ||
                _cz < _radius || _cz > _c - _radius)
            {
                throw std::runtime_error("Сфера выходит за границы RVE");
            }
        }
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

        // Фиксируем радиус, чтобы объемная доля оставалась постоянной
        double r = 0.18;

        double phi_f = geom::sphereVolumeFraction(r, RVE_a);
        double fiberVolumePercent = phi_f * 100.0;

        // Максимально допустимое смещение от центра:
        // центр куба = 0.5, сфера должна остаться внутри [r, 1-r]
        double maxShift = 0.5 - r;

        std::ofstream data("cte_vs_shifted.dat");
        if (!data)
        {
            std::cerr << "Не удалось создать cte_vs_shifted.dat\n";
            return 1;
        }

        data << "eta alpha_x alpha_y alpha_z cx cy cz vf shift\n";
        data << std::scientific << std::setprecision(10);

        for (double eta = 0.0; eta <= 1.0 + 1e-12; eta += 0.1)
        {
            double shift = eta * maxShift;

            double cx = 0.5;
            double cy = 0.5;
            double cz = 0.5 + shift;

            LOG(INFO) << "Volume fraction = " << fiberVolumePercent << "%" << std::endl;
            LOG(INFO) << "eta = " << eta
                        << ", shift = " << shift
                        << ", center = (" << cx << ", " << cy << ", " << cz << ")"
                        << std::endl;

            Problem problem(r, cx, cy, cz, RVE_a);
            problem.setHomMethod(AHM_SYM_CONST);
            problem.setName("ThermoElasticShifted");
            problem.init();
            problem.solve();

            data << eta << ' '
                    << problem.alphaX() * 1e6 << ' '
                    << problem.alphaY() * 1e6 << ' '
                    << problem.alphaZ() * 1e6 << ' '
                    << cx << ' '
                    << cy << ' '
                    << cz << ' '
                    << fiberVolumePercent << ' '
                    << shift << ' '
                    << '\n';
        }

        data.close();

        std::vector<plot::PlotSpec> plots =
        {
            {"results/vf_vs_cte_shifted_z.png", "Зависимость КЛТР от смещения сферы от центра вдоль оси OZ", plot::PlotType::Basic, {}, "Cмещение сферы относительно центра", "Эффективный КЛТР, 10^{-6}/K"}
        };

        if (!plot::generatePlots("cte_vs_shifted.dat", plots))
        {
            return 1;
        }

        LOG(INFO) << "Графики сохранены в директории results\n";

        return 0;
    }
}