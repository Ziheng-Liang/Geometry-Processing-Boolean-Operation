#include "polygon.h"
#include <Eigen/Dense>

#ifndef IGL_BOL_Node_H
#define IGL_BOL_Node_H

namespace igl
{
    namespace bol
    {
        struct Node {
            Node* left;
            Node* right;
            int size;
            Eigen::RowVectorXi index;
            Eigen::MatrixXi edges;
            std::vector<Polygon*> polygons;
        };

        void remove_edge_from_node(Node* node, int i, int j) {
        	for (int k = 0; k < node->size; k++) {
        		if ((node->edges(k, 0) == i && node->edges(k, 1) == j) || 
        			(node->edges(k, j) == i && node->edges(k, 1) == i)) {
        			node->edges.erase(node->edges.begin()+k);
        			node->size --;
        			break;
        		}
        	}
        	for (int k = 0; k < node->polygons.size(); k++) {
        		if (exist_edges(node->polygons.at(k), i, j)) {
        			node->polygons.erase(node->polygons.begin()+k);
        			break;
        		}
        	}
        }
    }
}

#endif