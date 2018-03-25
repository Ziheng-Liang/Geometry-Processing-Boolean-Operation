#include "my_node.h"

#ifndef IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H
#define IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H


namespace igl
{
  namespace bol
  {
  	void label(MatrixXr V, Eigen::MatrixXi F, Eigen::MatrixXi ODF);
    void constrained_delaunay_triangulation(MatrixXr V, Eigen::MatrixXi F);
    void subdivide(int start, int end, Node* node);
  }
}

