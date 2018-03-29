#include "my_node.h"
#include "boost_include.h"
#include <Eigen/Dense>

#ifndef IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H
#define IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H


namespace igl
{
  namespace bol
  {
  	void contruct_tree(Eigen::MatrixXd V, 
  					   Node* node);

    void delaunay_triangulation(Eigen::MatrixXd V, 
    										Eigen::MatrixXi C, 
    										Eigen::MatrixXi F);

    void subdivide(Eigen::RowVectorXi index, 
    			   Node* node);

    void build_edges(Eigen::MatrixXd V, 
    				 Node* node);

    double angle(Eigen::RowVectorXd p1, 
    			 Eigen::RowVectorXd p2, 
    			 Eigen::RowVectorXd p3);

    bool intersect(Eigen::RowVector3d p1, 
    			   Eigen::RowVector3d p2, 
    			   Eigen::RowVector3d p3, 
    			   Eigen::RowVector3d p4);

    void get_circle_center(Eigen::RowVector3d p1, 
    					   Eigen::RowVector3d p2, 
    					   Eigen::RowVector3d p3, 
    					   Eigen::RowVector3d center);

    double dmnop(Eigen::RowVector3d m, 
    		  Eigen::RowVector3d n, 
    		  Eigen::RowVector3d o, 
    		  Eigen::RowVector3d p);
  }
}

#endif