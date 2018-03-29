#include "my_node.h"
#include "boost_include.h"
#include <Eigen/Dense>

#ifndef IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H
#define IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H


namespace igl
{
  namespace bol
  {
    using namespace Eigen;
    using namespace std;
    void constrained_delaunay_triangulation(MatrixXd V, 
                                            MatrixXi C, 
                                            MatrixXi F);

  	void contruct_tree(MatrixXd V, 
  					   Node* node,
                       MatrixXd projectV);

    void add_constrained(MatrixXd V, 
                         MatrixXi C, 
                         Node* node);

    void break_polygons(MatrixXd V, 
                        RowVectorXi C, 
                        vector<Polygon*> polygons, 
                        RowVectorXi index);


    void delaunay_triangulation(MatrixXd V, 
    							Node* node);

    void subdivide(RowVectorXi index, 
    			   Node* node);

    double angle(RowVectorXd p1, 
    			 RowVectorXd p2, 
    			 RowVectorXd p3);

    bool intersect(RowVector3d p1, 
    			   RowVector3d p2, 
    			   RowVector3d p3, 
    			   RowVector3d p4);

    void get_circle_center(RowVector3d p1, 
    					   RowVector3d p2, 
    					   RowVector3d p3, 
    					   RowVector3d center);

    double dmnop(RowVector3d m, 
    		  RowVector3d n, 
    		  RowVector3d o, 
    		  RowVector3d p);
  }
}

#endif