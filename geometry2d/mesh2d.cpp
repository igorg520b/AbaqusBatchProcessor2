#include "mesh2d.h"
#include "element2d.h"
#include "cohesivezone2d.h"
#include "czinsertiontool2d.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>

#include <Eigen/Core>
#include <spdlog/spdlog.h>

icy::Node2D* icy::Mesh2D::AddNode()
{
    nodes.emplace_back();
    nodes.back().globId = (int)nodes.size()-1;
    return &nodes.back();
}

icy::Element2D* icy::Mesh2D::AddElement()
{
    elems.emplace_back();
    elems.back().elemId = (int)elems.size()-1;
    return &elems.back();
}

icy::CohesiveZone2D* icy::Mesh2D::AddCZ()
{
    czs.emplace_back();
    return &czs.back();
}

