#include <Eigen/Core>
#include "boost_include.h"
#include <algorithm>

namespace igl {
	namespace bol {
		void axis_aligned_boudning_box(MatrixXr V, 
								  Eigen::MatrixXi F, 
								  MatrixXr AABB);
		bool aabb_intersect(RowVectorXr First, 
							RowVectorXr Second);
	}
}