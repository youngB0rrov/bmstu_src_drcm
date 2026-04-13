#pragma once

#include <string>
#include <fstream>
#include <vector>

namespace plot
{
    enum class PlotType
    {
        Basic,
        Comparison
    };

    struct CurveSpec
    {
        int column;
        std::string label;
        std::string style = "linespoints";
        int lineWidth = 2;
        int pointType = 7;
        int dashType = 1;
    };

    struct PlotSpec
    {
        std::string outputFile;
        std::string title;
        PlotType plotType;
        std::vector<CurveSpec> curves;
        std::string xTitle = "Объемная доля включений, %";
        std::string yTitle = "Эффективный КЛТР, 10^{-6}/K";
    };

    void writeComparisonPlot(std::ofstream& gp, const std::string& dataFile, const PlotSpec& spec);
    void writePlotBase(std::ofstream &gp, const std::string &dataFile, const PlotSpec &spec);
    bool generatePlots(const std::string& dataFile, const std::vector<PlotSpec>& plots);
}

