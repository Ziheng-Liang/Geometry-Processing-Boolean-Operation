#ifndef IGL_BOL_Node_H
#define IGL_BOL_Node_H

#include "polygon.h"
#include <Eigen/Dense>

namespace igl
{
    namespace bol
    {
    	using namespace std;
        struct Node {
            Node* left;
            Node* right;
            int size;
            Eigen::VectorXi index;
            vector<tuple<int,int>> edges;
            vector<Polygon*> polygons;
        };

        inline void remove_edge_from_node(Node* node, int i, int j) {
        	for (int k = 0; k < node->size; k++) {
        		if ((get<0>(node->edges.at(k)) == i && get<1>(node->edges.at(k)) == j) || 
        			(get<1>(node->edges.at(k)) == i && get<0>(node->edges.at(k)) == i)) {
        			node->edges.erase(node->edges.begin()+k);
        			node->size --;
        			break;
        		}
        	}
        	for (int k = 0; k < node->polygons.size(); k++) {
        		if (exist_edges(node->polygons.at(k), i, j) != -1) {
        			node->polygons.erase(node->polygons.begin()+k);
        			break;
        		}
        	}
        }
    }
}

#endif