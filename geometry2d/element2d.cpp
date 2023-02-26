#include "element2d.h"
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <complex>
#include <cmath>


void icy::Element2D::Reset()
{
    elem_incident_faces.clear();
    elem_incident_faces.reserve(32);
}


void icy::Element2D::Precompute(std::vector<Node2D> &nodes)
{
    Eigen::Matrix2d Dm;
    Dm << nodes[nds[0]].x0 - nodes[nds[2]].x0, nodes[nds[1]].x0 - nodes[nds[2]].x0;

    double area_initial = Dm.determinant()/2;
    if(area_initial==0) throw std::runtime_error("element's initial area is zero");
    else if(area_initial < 0)
    {
        std::swap(nds[0],nds[1]);
    }
}


