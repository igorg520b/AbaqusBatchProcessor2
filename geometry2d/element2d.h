#ifndef ELEMENT2d_H
#define ELEMENT2d_H

#include "node2d.h"

namespace icy { struct Element2D; struct Node2D; }

struct icy::Element2D
{
    Element2D() {Reset();}

    int nds[3] {};
    int grainId;
    int elemId;
    std::vector<uint32_t> elem_incident_faces; // faces that are incident to any of elem's nodes

    void Reset();
    void Precompute(std::vector<Node2D> &nodes);

    constexpr static int fi[3][2] {{1,2},{2,0},{0,1}};
};

#endif // ELEMENT123_H
