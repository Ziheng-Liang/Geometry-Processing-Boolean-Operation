#ifndef IGL_BOL_BOOST_INCLUDE_H
#define IGL_BOL_BOOST_INCLUDE_H

#include <Eigen/Core>
#include <boost/multiprecision/cpp_int.hpp> 
#include <Eigen/LU>
#include <Eigen/Dense>
namespace igl {
	namespace bol{
		typedef boost::multiprecision::cpp_rational rat;
		typedef Eigen::Matrix<rat,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
		typedef Eigen::Matrix<rat,Eigen::Dynamic,3> MatrixX3r; //list of triangle
		typedef Eigen::Matrix<rat,3,3> Matrix33r; //This can a triangle
		typedef Eigen::Matrix<rat,Eigen::Dynamic,1> VectorXr;
		typedef Eigen::Matrix<rat, 1, Eigen::Dynamic> RowVectorXr;
		typedef Eigen::Matrix<rat, 1, 3> RowVector3r; // This can be a point
	}


}

namespace Eigen{
	template <> struct NumTraits<igl::bol::rat>: Eigen::GenericNumTraits<igl::bol::rat>{
	typedef igl::bol::rat Real;
    typedef igl::bol::rat NonInteger;
    typedef igl::bol::rat Nested;
    typedef igl::bol::rat Literal;

     enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 150,
      MulCost = 100
    	};
	};

}




#endif
