#include "node2d.h"


void icy::Node2D::Reset()
{
    x0.setZero();
    f.setZero();
    globId = -1;
    incident_faces.clear();
    surface = false;
}



