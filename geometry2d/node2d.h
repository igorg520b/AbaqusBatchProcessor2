#ifndef NODE2D_H
#define NODE2D_H

#include <functional>
#include <vector>
#include <Eigen/Core>

namespace icy { struct Node2D; struct Element2D; struct CohesiveZone2D; }

struct icy::Node2D
{
    Node2D() { Reset(); incident_faces.reserve(8); }

    void Reset();

    int globId;               // id in fragment; id in mesh; id in freenode list;
    bool surface;
    int group = 0;

    Eigen::Vector2d x0, f;    // pos-initial, pos-current, velocity-current, pos-tentative

    std::vector<uint32_t> incident_faces; // list of faces that connect to this node; to highest bits indicate face # in element

};

#endif
