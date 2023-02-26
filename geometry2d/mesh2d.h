#ifndef MESH2D_H
#define MESH2D_H

#include <vector>
#include <string>
#include <algorithm>

#include <gmsh.h>
#include <Eigen/Core>

#include <QString>
#include <QDir>
#include <QDebug>

#include "node2d.h"

namespace icy { class Mesh2D; struct Node2D; struct Element2D; struct CohesiveZone2D;}

class icy::Mesh2D
{
public:

    std::vector<icy::Node2D> nodes;
    std::vector<icy::Element2D> elems;
    std::vector<icy::CohesiveZone2D> czs;

    icy::Node2D* AddNode();
    icy::Element2D* AddElement();
    icy::CohesiveZone2D* AddCZ();
};
#endif
