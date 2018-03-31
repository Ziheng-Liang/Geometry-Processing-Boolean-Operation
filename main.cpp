#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 
#include <boost/multiprecision/cpp_int.hpp> 
#include "boost_include.h"
// #include "triangle.h"
#include "constrained_delaunay_triangulation.h"
// #include "constrained_delaunay_triangulation.cpp"
#include "polygon.h"

// motivation
// current
// future work

int main()
{
  using namespace boost::multiprecision;
  using namespace igl::bol;
  Eigen::MatrixXd V(4,2);
  V<<0,0,
     1,0,
     1,1,
     0,1;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     0,2,3;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  std::cout<<"Hello, mesh: "<<std::endl<<L*V<<std::endl;
  cpp_rational v = 0.5;
  std::cout<< v << std::endl;
  VectorXr vt(3);
  vt(0) = 0.5;
  vt(1) = 0.5;
  vt(2) = 0.5;
  std::cout << vt << std::endl;
  Eigen::MatrixXi temp;
  Eigen::MatrixXi F1;
  constrained_delaunay_triangulation(V, temp, F1);

  // MatrixXr vnew = V.cast<cpp_rational> ();

  // std::cout << vnew << std::endl;

  // Matrix33r A, B;
  // A(0,0) = 0;
  // A(0,1) = 0;
  // A(0,2) = 0;

  // A(1,0) = 0;
  // A(1,1) = 1;
  // A(1,2) = 0;

  // A(2,0) = 1;
  // A(2,1) = 0;
  // A(2,2) = 0;

  // B(0,0) = 0;
  // B(0,1) = 0;
  // B(0,2) = -1;

  // B(1,0) = 2;
  // B(1,1) = 2;
  // B(1,2) = 0;

  // B(2,0) = 0;
  // B(2,1) = 0;
  // B(2,2) = 1;

   
  // std::vector<RowVector3r> r =igl::bol::t2t_intersect(A, B);
  // for (int i = 0; i<r.size(); i++){
  //   std::cout << r[i] << std::endl;
  // }
  return 0; 
}
