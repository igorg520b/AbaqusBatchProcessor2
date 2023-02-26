#include "node2d.h"


void icy::Node2D::Reset()
{
    x0.setZero();
    globId = -1;
    pinned = false;
    incident_faces.clear();
    surface = false;
}

void icy::Node2D::InitializeFromAnother(icy::Node2D *other)
{
    x0 = other->x0;
    pinned = other->pinned;
}


