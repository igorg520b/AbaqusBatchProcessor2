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
        static std::tuple<int,int> make_key(icy::Node2D *nd0, icy::Node2D *nd1)
        {
            int nds[2] = {nd0->globId,nd1->globId};
            std::sort(std::begin(nds),std::end(nds));
            return std::tuple<int,int> (nds[0],nds[1]);
        }
    };

    void InsertCZs(icy::Mesh2D &mesh);

private:
    void ReplaceNodeInElement(icy::Element2D *elem, icy::Node2D *whichNode, icy::Node2D *replacement);
};

#endif // CZINSERTIONTOOL_H
