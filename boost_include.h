#ifndef IGL_BOL_BOOST_INCLUDE_H
#define IGL_BOL_BOOST_INCLUDE_H

#include <Eigen/Core>
#include <boost/multiprecision/cpp_int.hpp> 

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

#endif
