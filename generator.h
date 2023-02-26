#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include "mesh2d.h"
#include <Eigen/Core>

#include <QDir>
#include <QFileInfo>

class Generator
{
public:

    double indenterRadius = 0.161925;
    double indenterOffset = 0;

    double interactionRadius = 0.01; // for exponential interaction property

    double timeToRun = 0.1;
    int nFrames = 400;

    bool loadWithIndenter = true;   // if false -> static load

    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 4e6;
    constexpr static double czElasticity = 1e11;
    constexpr static double czEnergy = 100;

    icy::Mesh2D mesh2d;

    void LoadFromFileWithCrop(std::string MSHFileName, double indenterX, double indenterY, int cropType);
private:
    void CreatePy2D(std::string outputFileName);
};

#endif // GENERATOR_H

