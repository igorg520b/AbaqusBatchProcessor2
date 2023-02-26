#include "czinsertiontool2d.h"
#include "element2d.h"
#include "cohesivezone2d.h"
#include <iostream>
#include <map>
#include <vector>
#include <spdlog/spdlog.h>

void icy::CZInsertionTool2D::InsertCZs(icy::Mesh2D &mesh)
{
    std::size_t nNodes = mesh.nodes.size();

    std::vector<NodeExtension> nExt;
    nExt.resize(nNodes);

    for(icy::Element2D *elem : mesh.elems)
    {
        for(int k=0;k<3;k++)
        {
            nExt[elem->nds[k]->globId].adj_elems.push_back(elem);
            nExt[elem->nds[k]->globId].adj_grains.push_back(elem->grainId);
        }
    }

    int connected_to_grain = 0, disconnected = 0;
    for(NodeExtension &ndex : nExt)
    {
        auto &v = ndex.adj_grains;
        std::sort(v.begin(),v.end());
        v.erase(std::unique(v.begin(),v.end()), v.end());
        if(v.size()==1) disconnected++;
        else if(v.size()>1) connected_to_grain++;
        else throw std::runtime_error("disconnected grain");
    }

    std::cout << "nNodes " << nNodes << std::endl;
    std::cout << "connected_to_grain " << connected_to_grain << std::endl;
    std::cout << "disconnected " << disconnected << std::endl;

    // FIND ALL FACETS
    std::map<std::tuple<int,int>,icy::CZInsertionTool2D::Facet> facets;

    for(icy::Element2D *e : mesh.elems)
    {
        for(int k=0;k<3;k++)
        {
            std::tuple<int,int> key = Facet::make_key(e->nds[Element2D::fi[k][0]], e->nds[Element2D::fi[k][1]]);
            Facet facet;
            facet.key = key;
            facet.elems[0] = e;
            facet.facet_idx[0] = k;
            auto result = facets.insert({key,facet});
            if(result.second == false)
            {
                Facet &f = result.first->second;
                f.elems[1] = e;
                f.facet_idx[1] = k;
            }
        }
    }

    int facets_exterior = 0, facets_interior = 0, facets_unassigned = 0, facets_connecting_grains = 0;

    for(auto &kvp : facets)
    {
        Facet &f = kvp.second;
        if(f.elems[1] != nullptr) facets_interior++;
        else if(f.elems[0] != nullptr) facets_exterior++;
        else facets_unassigned++;

        if(f.elems[1] != nullptr && f.elems[0]->grainId != f.elems[1]->grainId)
        {
            facets_connecting_grains++;
            // synchronize the nodes
            icy::Node2D *ref_node = f.elems[0]->nds[Element2D::fi[f.facet_idx[0]][0]];
            f.orientation = -1;
            for(int k=0;k<2;k++)
                if(ref_node == f.elems[1]->nds[Element2D::fi[f.facet_idx[1]][k]]) f.orientation = k;
            if(f.orientation == -1)
            {
                spdlog::critical("no matching node in facet");
                throw std::runtime_error("no matching node in facet");
            }

            for(int k=0;k<2;k++)
                if(f.elems[0]->nds[Element2D::fi[f.facet_idx[0]][k]] !=
                        f.elems[1]->nds[Element2D::fi[f.facet_idx[1]][(2-k+f.orientation)%2]])
                {
                    spdlog::critical("facet matching error");
                    throw std::runtime_error("facet matching error");
                }
        }
    }
    std::cout << "facets_exterior " << facets_exterior << std::endl;
    std::cout << "facets_interior " << facets_interior << std::endl;
    std::cout << "facets_unassigned " << facets_unassigned << std::endl;
    std::cout << "facets_connecting_grains " << facets_connecting_grains << std::endl;



    // for each node whose adj_grains.size()>1,
    // insert adj_grains.size()-1 new nodes
    // next: in each adj_elems, substitute the node

    for(int i=0;i<nNodes;i++)
    {
        NodeExtension &ndex = nExt[i];
        std::vector<int> &v = ndex.adj_grains;
        if(v.size() == 1) continue;

        std::map<int, icy::Node2D*> insertedNodes;
        icy::Node2D *nd = mesh.nodes[i];
        insertedNodes.insert({v[0],nd});

        for(int i=1;i<v.size();i++)
        {
            icy::Node2D *insertedNode = mesh.AddNode();
            insertedNode->InitializeFromAnother(nd);
            insertedNodes.insert({v[i],insertedNode});
        }

        for(icy::Element2D *elem : ndex.adj_elems)
        {
            ReplaceNodeInElement(elem, nd, insertedNodes.at(elem->grainId));
        }

    }
    spdlog::info("cz nodes inserted; nodes.size() {}",mesh.nodes.size());

    // create cohesive zones (but the nodes are still fused)
    mesh.czs.reserve(facets_connecting_grains);
    for(auto &kvp : facets)
    {
        Facet &f = kvp.second;

        if(f.elems[1] != nullptr && f.elems[0]->grainId != f.elems[1]->grainId)
        {

            icy::CohesiveZone2D *cz = mesh.AddCZ();
            icy::Element2D *e0 = f.elems[0];
            icy::Element2D *e1 = f.elems[1];
            cz->elems[0] = e0;
            cz->elems[1] = e1;
            cz->faceIds[0] = (uint8_t)f.facet_idx[0];
            cz->faceIds[1] = (uint8_t)f.facet_idx[1];
            cz->nds[0] = e0->nds[Element2D::fi[f.facet_idx[0]][0]];
            cz->nds[1] = e0->nds[Element2D::fi[f.facet_idx[0]][1]];

            cz->nds[2] = e1->nds[Element2D::fi[f.facet_idx[1]][(0+f.orientation)%2]];
            cz->nds[3] = e1->nds[Element2D::fi[f.facet_idx[1]][(1+f.orientation)%2]];

            for(int k=0;k<2;k++)
                if(cz->nds[k]->x0 != cz->nds[k+2]->x0)
                {
                    for(int j=0;j<4;j++) spdlog::critical("nd {}; {},{}",j,cz->nds[j]->x0[0],cz->nds[j]->x0[1]);
                    spdlog::critical("f.orientation {}", f.orientation);
                    throw std::runtime_error("cz insertion error  - 2D");
                }
        }
    }
    spdlog::info("czs inserted; czs.size() {}",mesh.czs.size());
}

void icy::CZInsertionTool2D::ReplaceNodeInElement(icy::Element2D *elem, icy::Node2D *whichNode, icy::Node2D *replacement)
{
    if(whichNode == replacement) return;
    for(int i=0;i<3;i++)
    {
        if(elem->nds[i]==whichNode)
        {
            elem->nds[i] = replacement;
            return;
        }
    }
    throw std::runtime_error("ReplaceNodeInElement: node not found");
}
