#include "boost_include.h"
#include <Eigen/Core>

namespace igl {
	namespace bol {
		bool point_on_plane(RowVectorXr x1, 
							RowVectorXr x2, 
							RowVectorXr x3, 
							RowVectorXr x4);
	}
}