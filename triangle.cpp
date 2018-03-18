#include "triangle.h"
#include <vector>






//http://web.mst.edu/~chaman/home/pubs/2015WimoTriangleTrianglePublished.pdf
//A : 3 x 3, each row represents the vectics 
std::vector<RowVector3r> intersect(const Matrix33r & A, const Matrix33r & B){
	RowVector3r pA = A.row(0);
	RowVector3r U = A.row(1) - pA;
	RowVector3r V = A.row(2) - pA;

	RowVector3r pB = B.row(0);
	RowVector3r S = B.row(1) - pB;
	RowVector3r T = B.row(2) - pB;

	RowVector3r AP = B.row(0) - A.row(0);

	rat l = (U.cross(V)).dot(U.cross(V));
	RowVector3r alpha = S.cross(U.cross(V)) / l;
	RowVector3r beta =  T.cross(U.cross(V)) / l;
	RowVector3r gamma = AP.cross(U.cross(V)) / l;

	- gamma.dot(U);
	alpha.dot(U);
	beta.dot(U);
	1 - gamma.dot(U);
	
	-1 - gamma.dot(V);
	alpha.dot(V);
	beta.dot(V);
	- gamma.dot(V);

	-gamma.dot(U - V);
	alpha.dot(U-V);
	beta.dot(U-V);
	1 - gamma.dot( U -V);


	


}
