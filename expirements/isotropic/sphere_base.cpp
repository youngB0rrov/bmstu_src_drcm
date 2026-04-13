#include "ThermoElasticHomProblem.h"
#include <fstream>
#include "FEReaders.h"
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "GeometryUtils.h"
#include "AnalyticModels.h"

#define _USE_MATH_DEFINES

using namespace nla3d;
using namespace hom3d;

namespace expirements::base
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
		double RVE_a = 1.0;

		std::ofstream data("cte_vs_vf.dat");

		if (!data)
		{
			std::cerr << "Не удалось создать cte_vs_vf.dat\n";
			return 1;
		}

		data << "vf alpha_x alpha_y alpha_z\n";
		data << std::scientific << std::setprecision(10);

		for (double r = 0.05; r < RVE_a / 2.0; r+=0.05)
		{
			double fiberVolumeFrac = geom::sphereVolumeFraction(r, RVE_a);
			double fiberVolumePercent = fiberVolumeFrac * 100;

			LOG(INFO) << "Fiber vol.%: " << fiberVolumePercent << "%" << std::endl;

			Problem problem(r, RVE_a);
			problem.setHomMethod(AHM_SYM_CONST);
			problem.setName("ThermoElastic");
			problem.init();
			problem.solve();
			problem.printResults();

			data << fiberVolumePercent << ' ' << problem.alphaX() * 1e6 << ' ' << problem.alphaY() * 1e6 << ' ' << problem.alphaZ() * 1e6 << '\n';
		}

		data.close();

		std::ofstream gp("plot.gp");
		if (!gp) {
			std::cerr << "Не удалось создать plot.gp\n";
			return 1;
		}

		gp <<
			"set term pngcairo size 1200,800\n"
			"set output 'vf_vs_cte.png'\n"
			"set encoding utf8\n"
			"set datafile columnheaders\n"
			"set title 'Зависимость КЛТР от объемной доли включений'\n"
			"set xlabel 'Объемная доля включений, %'\n"
			"set ylabel 'Эффективный КЛТР, 1 * 10^{-6}/K'\n"
			"set grid\n"
			"set key outside right top\n"
			"\n"
			"plot \\\n"
			"  'cte_vs_vf.dat' using ($1):2 with linespoints lw 2 pt 7 title 'alpha_x', \\\n"
			"  'cte_vs_vf.dat' using ($1):3 with linespoints lw 2 pt 5 title 'alpha_y', \\\n"
			"  'cte_vs_vf.dat' using ($1):4 with linespoints lw 2 pt 9 title 'alpha_z'\n";

		gp.close();

		int rc = std::system("gnuplot plot.gp");
		if (rc != 0) {
			std::cerr << "Ошибка запуска gnuplot.\n";
			return 1;
		}

		LOG(INFO) << "Cоздан файл vf_vs_cte.png с результатами расчета\n";

		return 0;
	}
}