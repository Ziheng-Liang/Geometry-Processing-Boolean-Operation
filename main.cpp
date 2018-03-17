#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 
#include <boost/multiprecision/cpp_int.hpp> 
#include "boost_include.h"
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

  return 0;
}
