#ifndef IGL_BOL_TRIANGLE_H
#define IGL_BOL_TRIANGLE_H

#include <Eigen/Core>
#include "boost_include.h"
#include <vector>
namespace igl {
	namespace bol{
		//Check if a triagnle is degenerate or not
		inline bool is_degenerate(Matrix33r V);

		//Intersect triangle A and B
		inline std::vector<RowVector3r> t2t_intersect(const Matrix33r & A, const Matrix33r & B);

		inline void t2t_intersect_on_A(const Matrix33r & A, const Matrix33r & B, MatrixXr & AV, Eigen::MatrixXi & AF);

		inline Eigen::RowVector3d rat_to_double(RowVector3r & to_cast);

	}
}

#ifndef IGL_STATIC_LIBRARY
#  include "triangle.cpp"
#endif


#endif