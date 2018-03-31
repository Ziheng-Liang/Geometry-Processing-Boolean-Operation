#ifndef IGL_BOL_CONVEXHULL_H
#define IGL_BOL_CONVEXHULL_H

#include <Eigen/Core>
#include "boost_include.h"
#include <vector>
namespace igl {
	namespace bol{
		//https://en.wikipedia.org/wiki/Graham_scan
		void convexhull(const MatrixXr & V, MatrixXr & H);

	}
}



#endif