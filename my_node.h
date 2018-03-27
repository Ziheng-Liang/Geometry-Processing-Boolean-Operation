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
    }
}

#endif