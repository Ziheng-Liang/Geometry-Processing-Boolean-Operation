#ifndef IGL_BOL_Node_H
#define IGL_BOL_Node_H

#include "polygon.h"
#include <Eigen/Dense>
#include <iostream>

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
        	std::cout << "remove_edge_from_node1" << std::endl;
        	for (int k = 0; k < node->size; k++) {
        		std::cout << k << std::endl;
        		std::cout << node->edges.size() << std::endl;
        		if ((get<0>(node->edges.at(k)) == i && get<1>(node->edges.at(k)) == j) || 
        			(get<1>(node->edges.at(k)) == i && get<0>(node->edges.at(k)) == i)) {
        			node->edges.erase(node->edges.begin()+k);
        			node->size --;
        			break;
        		}
        	}
        	std::cout << "remove_edge_from_node2" << std::endl;
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