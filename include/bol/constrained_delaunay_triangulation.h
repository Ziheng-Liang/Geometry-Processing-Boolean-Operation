#ifndef IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H
#define IGL_BOL_CONSTRAINED_DELAUNAY_TRIANGULATION_H

#include "my_node.h"
#include "polygon.h"
#include <Eigen/Dense>


namespace igl
{
  namespace bol
  {
    using namespace Eigen;
    using namespace std;
    
    inline void constrained_delaunay_triangulation(const MatrixXd &V, 
                                            const MatrixXi &C, 
                                            MatrixXi &F);

  	inline void contruct_tree(const MatrixXd &V, 
  					           Node* node,
                       MatrixXd &projectV);

    inline void add_constrained(MatrixXd V, 
                         MatrixXi C, 
                         Node* node);

    inline void break_polygons(MatrixXd V, 
                        RowVectorXi C, 
                        vector<Polygon*> polygons, 
                        RowVectorXi index);


    inline void delaunay_triangulation(MatrixXd V, 
    							Node* node);

    inline void subdivide(VectorXi index, 
    			   Node* node);

    inline double angle(RowVectorXd p1, 
    			 RowVectorXd p2, 
    			 RowVectorXd p3);

    inline bool intersect(const RowVectorXd &p1, 
                   const RowVectorXd &p2, 
                   const RowVectorXd &p3, 
                   const RowVectorXd &p4);

    inline void get_circle_center(const RowVectorXd &p1, 
                           const RowVectorXd &p2, 
                           const RowVectorXd &p3, 
                           RowVectorXd &center);

    inline double dmnop(RowVector3d m, 
    		  RowVector3d n, 
    		  RowVector3d o, 
    		  RowVector3d p);

    inline bool colinear(const RowVectorXd &p1, 
                  const RowVectorXd &p2, 
                  const RowVectorXd &p3);
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "constrained_delaunay_triangulation.cpp"
#endif

#endif