#include "axis_aligned_boudning_box.h"

void axis_aligned_boudning_box(MatrixXr V, Eigen::MatrixXi F, MatrixXr AABB) {
	using namespace Eigen;
	using namespace std;
	int face_num = F.rows();
	AABB = MatrixXi(face_num, 6);
	for (int i = 0; i < face_num; i++) {
		for (int j = 0; j < 3; j++) {
			AABB(i, j) = min(V(F(i, 0), j), min(V(F(i, 1), j), V(F(i, 2), j)));
			AABB(i, j + 3) = max(V(F(i, 0), j), max(V(F(i, 1), j), V(F(i, 2), j)));
		}
	}
}

bool aabb_intersect(RowVectorXr First, RowVectorXr Second) {
	for (int i = 0; i < 3; i++) {
		if (First(i) > Second(i + 3)) 
			return false;
		if (Second(i) > First(i + 3))
			return false
	}
	return true;
}

