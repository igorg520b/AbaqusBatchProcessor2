#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include "mesh2d.h"
#include <Eigen/Core>

#include <QDir>
#include <QFileInfo>

#include <cmath>
#include <iomanip>

constexpr double tekscanRowToAngle(int row) { return (3.14159265358979323846/2 - 1.27713*(38-row)/44.); }


class Generator
{
public:

    constexpr static double indenterRadius = 0.161925;
    constexpr static double indentationDepth = 0.1016; // 4" in meters
    constexpr static double forceMagnitude = 1e6;
    constexpr static double pi = 3.14159265358979323846;
    constexpr static double cropStart = tekscanRowToAngle(10);
    constexpr static double cropEnd = tekscanRowToAngle(20);
    constexpr static double forceAngle = 40*pi/180; // (cropStart+cropEnd)/2;

    constexpr static bool useNormalDistribution = false;
    constexpr static bool attachSides = false; // generate blocks with boundary condition on the sides

    constexpr static double timeToRun = 0.25;
    constexpr static int nFrames = 560;

    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 4e6;
    constexpr static double czElasticity = 1e11;
    constexpr static double czEnergy = 100;

    double indenterX, indenterY; // initialized in LoadFromFileWithCrop
    int setupType; // initialized in LoadFromFileWithCrop
    double fX, fY;  // force components for load

    icy::Mesh2D mesh2d;

    std::string AbaqusExecutablePath = "G:\\SIMULIA\\Commands\\abaqus.bat";

    void LoadFromFileWithCrop(std::string MSHFileName, double indenterX, double indenterY, int cropType);
    void AddEntryToJobGeneratorBat(std::ofstream &s);
    void AddEntryToJobExecutionBat(std::ofstream &s);
    void AddEntryToBinaryExportBat(std::ofstream &s);
private:
    std::string fileName;
    void CreatePy2D();
    void CreateExportScript();

    static double uniformDistribution(const double s, const double e, const double x);
    static double normalDistribution(const double sigma, const double mu, const double x);
    std::vector<int> forcedNodes; // a list of the nodes to which the force is applied (nd.group==5)
};

#endif // GENERATOR_H

