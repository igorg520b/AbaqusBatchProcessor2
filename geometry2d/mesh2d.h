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
#include "element2d.h"
#include "cohesivezone2d.h"

namespace icy { struct Mesh2D; }

struct icy::Mesh2D
{
    std::vector<icy::Node2D> nodes;
    std::vector<icy::Element2D> elems;
    std::vector<icy::CohesiveZone2D> czs;

    icy::Node2D* AddNode();
    icy::Element2D* AddElement();
    icy::CohesiveZone2D* AddCZ();
};
#endif
