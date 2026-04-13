#include "PlotUtils.h"
#include "FEReaders.h"
#include <filesystem>
#include <vector>

namespace plot
{
    void writeComparisonPlot(std::ofstream &gp, const std::string &dataFile, const PlotSpec &spec)
    {
        gp << "set term pngcairo size 1200,800\n";
        gp << "set output '" << spec.outputFile << "'\n";
        gp << "set encoding utf8\n";
        gp << "set datafile columnheaders\n";
        gp << "set title '" << spec.title << "'\n";
        gp << "set xlabel '" << spec.xTitle << "'\n";
        gp << "set ylabel '" << spec.yTitle << "'\n";
        gp << "set grid\n";
        gp << "set key outside right top\n\n";

        gp << "plot \\\n";

        for (size_t i = 0; i < spec.curves.size(); ++i)
        {
            const auto &c = spec.curves[i];

            gp << "  '" << dataFile << "' using 1:" << c.column
            << " with " << c.style
            << " lw " << c.lineWidth
            << " pt " << c.pointType
            << " dt " << c.dashType
            << " title '" << c.label << "'";

            if (i + 1 != spec.curves.size())
                gp << ", \\\n";
            else
                gp << "\n";
        }

        gp << "\n";
    }

    void writePlotBase(std::ofstream &gp, const std::string &dataFile, const PlotSpec &spec)
    {
        gp <<
            "set term pngcairo size 1200,800\n"
            "set output '" << spec.outputFile << "'\n"
            "set encoding utf8\n"
            "set datafile columnheaders\n"
            "set title '" << spec.title << "'\n"
            "set xlabel '" << spec.xTitle << "'\n"
            "set ylabel '" << spec.yTitle << "'\n"
            "set grid\n"
            "set key outside right top\n"
            "\n"
            "plot \\\n"
            "  '" << dataFile << "' using 1:2 with linespoints lw 2 pt 7 title 'alpha_x', \\\n"
            "  '" << dataFile << "' using 1:3 with linespoints lw 2 pt 5 title 'alpha_y', \\\n"
            "  '" << dataFile << "' using 1:4 with linespoints lw 2 pt 9 title 'alpha_z'\n";
    }

    bool generatePlots(const std::string &dataFile, const std::vector<PlotSpec>& plots)
    {
        std::filesystem::create_directories("results");

        std::ofstream gp("plot.gp");
        if (!gp)
        {
            LOG(ERROR) << "Не удалось создать plot.gp\n";
            return false;
        }

        for (const auto& spec : plots)
        {
            if (spec.plotType == PlotType::Basic)
            {
                writePlotBase(gp, dataFile, spec);
            }
            else
            {
                writeComparisonPlot(gp, dataFile, spec);
            }
        }

        gp.close();

        int rc = std::system("gnuplot plot.gp");
        if (rc != 0)
        {
            LOG(ERROR) << "Ошибка запуска gnuplot.\n";
            return false;
        }

        return true;
    }
}