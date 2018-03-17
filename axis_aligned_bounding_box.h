#include <Eigen/Core>
#include <boostinclude>

namespace igl {
	namespace bol {
		void axis_aligned_boudning_box(MatrixXr V, 
								  Eigen::MatrixXi F, 
								  MatrixXr AABB);
		bool aabb_intersect(RowVectorXr First, 
							RowVectorXr Second);
	}
}