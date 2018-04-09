#include "supporting_plane.h"

bool point_on_plane(RowVectorXr x1, RowVectorXr x2, RowVectorXr x3, RowVectorXr x4) {
	using namespace Eigen;
	MatrixXd coplane = MatrixXd::Zero(4, 4);
	coplane << x1(0), x2(0), x3(0), x4(0),
			   x1(1), x2(1), x3(1), x4(1),
			   x1(2), x2(2), x3(2), x4(2),
			   1,     1,     1,     1;
	return coplane.determinant() == 0;
}