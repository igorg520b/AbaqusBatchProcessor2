#ifndef CZINSERTIONTOOL_H
#define CZINSERTIONTOOL_H

#include "node2d.h"
#include "mesh2d.h"
#include <tuple>
#include <algorithm>

namespace icy {class CZInsertionTool2D;}

class icy::CZInsertionTool2D
{
public:
    struct NodeExtension
    {
        std::vector<Element2D*> adj_elems;
        std::vector<int> adj_grains;
    };

    struct Facet
    {
        icy::Element2D *elems[2] {};
        int facet_idx[2] {};
        int orientation; // index of the node in facet 1 that matches node 0 in facet 0
        std::tuple<int,int> key;
        static std::tuple<int,int> make_key(int nd0idx, int nd1idx)
        {
            if(nd0idx > nd1idx) return std::tuple<int,int> (nd1idx,nd0idx);
            else return std::tuple<int,int> (nd0idx,nd1idx);
        }
    };

    void InsertCZs(icy::Mesh2D &mesh);

private:
    void ReplaceNodeInElement(icy::Element2D *elem, int whichNode, int replacement);
};

#endif // CZINSERTIONTOOL_H
