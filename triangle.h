#ifndef IGL_BOL_TRIANGLE_H
#define IGL_BOL_TRIANGLE_H

#include <Eigen/Core>
#include "boost_include.h"
#include <vector>
namespace igl {
	namespace bol{
		//Check if a triagnle is degenerate or not
		bool is_degenerate(Matrix33r V);

		//Intersect triangle A and B
		std::vector<RowVector3r> t2t_intersect(const Matrix33r & A, const Matrix33r & B);
	}
}



#endif