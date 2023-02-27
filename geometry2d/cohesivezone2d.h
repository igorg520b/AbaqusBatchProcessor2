#ifndef COHESIVEZONE2d_H
#define COHESIVEZONE2d_H

#include "node2d.h"
#include <cstdint>

namespace icy {struct CohesiveZone2D; struct Node2D; struct Element2D;}

struct icy::CohesiveZone2D
{
    CohesiveZone2D() {Reset();}
    int elems[2];                  // each CZ connects two elements
    uint8_t faceIds[2];
    int nds[4];
    void Reset();
};

#endif // COHESIVEZONE_H
