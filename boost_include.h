#ifndef IGL_BOL_BOOST_INCLUDE_H
#define IGL_BOL_BOOST_INCLUDE_H

#include <Eigen/Core>
#include <boost/multiprecision/cpp_int.hpp> 

namespace igl {
	namespace bol{
		typedef Eigen::Matrix<boost::multiprecision::cpp_rational,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
		typedef Eigen::Matrix<boost::multiprecision::cpp_rational,Eigen::Dynamic,3> MatrixX3r;
		typedef Eigen::Matrix<boost::multiprecision::cpp_rational,Eigen::Dynamic,1> VectorXr;
		typedef Eigen::Matrix<boost::multiprecision::cpp_rational, 1, Eigen::Dynamic> RowVectorXr;

	}


}

#endif
