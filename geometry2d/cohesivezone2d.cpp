#include "cohesivezone2D.h"
#include <iomanip>
#include <Eigen/Geometry>

void icy::CohesiveZone2D::Reset()
{
    elems[0] = elems[1] = -1;
    for(int i=0;i<4;i++) nds[i] = -1;
}
