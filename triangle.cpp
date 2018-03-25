#include "triangle.h"
#include <vector>
#include <iostream>
using namespace igl::bol;
bool igl::bol::is_degenerate(Matrix33r V){
	// return (V.row(0).array() == V.row(1).array() ||
	// 		V.row(0).array() == V.row(1).array() ||
	// 		V.row(1).array() == V.row(2).array());
	return false;
}

static int find_different_sign(rat d[3]){
	for (int i = 0; i < 3; i ++){
		if (d[i] * d[(i+1)%3] < 0 && d[i] * d[(i+2)%3] < 0){
			return i;
		}
	}
	assert(false);
}

static int find_zero(rat d[3]){
	for (int i = 0; i < 3; i++){
		if (d[i] == 0) return i;
	}
	assert(false);
}

static int find_non_zero(rat d[3]){
	for (int i = 0; i < 3; i++){
		if (d[i] != 0) return i;
	}
	assert(false);
}


//Find the itersection of two lines
//Precondition: two lines coplanar
static void l2l_intersection(const RowVector3r &d1, const RowVector3r &x1, 
	const RowVector3r &d2, const RowVector3r &x2, RowVector3r & p){
	//Line 1: = x1 + t * d1;
	//Line 2: = x2 + s * d2
	//Solve this linear equation (over determined)
	//We know its on the same plane before hand, therefore we just need two axis
	MatrixXr left;
	left.resize(2,2);
	left(0,0) = d1(0,0);
	left(0,1) = - d2(0,0);
	left(1,0) = d1(0,1);
	left(1,1) = - d2(0,1);
	VectorXr right;
	right.resize(2,1);
	right(0,0) = x2(0,0) - x1(0,0);
	right(1,0) = x2(0,1) - x1(0,1);
	
	auto ts = left.inverse() * right;
	auto t = ts(0,0);
	p = x1 + t * d1;
}

//Find a point from intersection of two plane
static void find_point_from_two_plane(const RowVector3r & n1, const rat & d1, 
	const RowVector3r & n2, const rat & d2, RowVector3r p){
	//
	for (int i = 0; i < 3; i++){
		int a0 = i;
		int a1 = (i+1)%3;
		int a2 = (i+2)%3;
		if ((n1(a0) * n2(a1) - n1(a1) * n2(a0)) != 0){ //if it is invertible
			//Then a2 is good to assume,
			rat fixed = 1; //set to 1; or anything

			//Find a point on the line, suppose z = 1
			MatrixXr left;
			left.resize(2,2);
			left(0,0) = n1(0,a0);
			left(0,1) = n1(0,a1);
			left(1,0) = n2(0,a0);
			left(1,1) = n2(0,a1);
			VectorXr right;
			right.resize(2,1);
			right(0,0) = -d1 - n1(0,a2);
			right(1,0) = -d2 - n2(0,a2);
			MatrixXr result = left.inverse() * right;

			p(0,a0) = result(0,0);
			p(0,a1) = result(1,0);
			p(0,a2) = 1;
			return;
		}
	}
	assert(false);
}

//For solving 2x2 system of linear equaltions
static void solve2x2(rat a0, rat a1, rat a2, rat a3, rat b1, rat b2){

}


//Find the t value that cooresponds a poing on line.
//Precondition: p is on the line, dicretion d is non zero
static rat point_on_line(const RowVector3r & p, const RowVector3r & x, const RowVector3r &d){
	RowVector3r v = p - x;
	for (int i = 0; i < 3; i++){
		if (d(0,i) != 0) return v(0,i)/d(0,i);
	}
}


//https://www.tandfonline.com/doi/pdf/10.1080/10867651.1997.10487472?needAccess=true
//A : 3 x 3, each row represents the vectics 
std::vector<RowVector3r> igl::bol::t2t_intersect(const Matrix33r & A, const Matrix33r & B){
	assert(!is_degenerate(A));
	assert(!is_degenerate(B));

	RowVector3r nA = (A.row(1) - A.row(0)).cross(A.row(2) - A.row(0));
	rat dA = -nA.dot(A.row(0));
	rat sdB2A[3];

	int num_sd_B2A0 = 0;
	for (int i = 0; i < 3; i++){
		sdB2A[i] = nA.dot(B.row(i)) + dA;
		if (sdB2A[i] == 0) num_sd_B2A0++;
	}


	RowVector3r nB = (B.row(1) - B.row(0)).cross(B.row(2) - B.row(0));
	rat dB = -nB.dot(A.row(0));
	rat sdA2B[3];
	int num_sd_A2B0 = 0;
	for (int i = 0; i < 3; i++){
		sdA2B[i] = nB.dot(A.row(i)) + dB;
		if (sdA2B[i] == 0) num_sd_A2B0++;
	}


	std::cout << nA << std::endl;
	std::cout << nB << std::endl;
	std::vector<RowVector3r> return_v;

	if (sdB2A[0] == 0 && sdB2A[1] == 0 && sdB2A[2] == 0){ //coplanar


	} else {
		if ((sdB2A[0] > 0 && sdB2A[1] > 0 && sdB2A[2] > 0) 
			|| (sdB2A[0] < 0 && sdB2A[1] < 0 && sdB2A[2] < 0)){//no intersection
			return return_v;
		} else if ((sdA2B[0] > 0 && sdA2B[1] > 0 && sdA2B[2] > 0) 
			|| (sdA2B[0] < 0 && sdA2B[1] < 0 && sdA2B[2] < 0)){ //no intersection
			return return_v;
		} else { //possible intersection
			//Line direction: 
			//We cant use square root, since its on rational space
			RowVector3r ld = nA.cross(nB);

			RowVector3r O;
			find_point_from_two_plane(nA, dA, nB, dB, O);

			std::cout << "Here" << std::endl;
			//Do interval test on A0
			rat tA0, tA1;
			RowVector3r pA0, pA1;
			{
				assert(num_sd_A2B0 < 3);
				if (num_sd_A2B0 == 2){
					int v = find_non_zero(sdA2B);
					int v1 = (v+1)%3;
					int v2 = (v+2)%3;
					pA0 = A.row(v1);
					pA1 = A.row(v2);
				} else if (num_sd_A2B0 == 1){
 					int v = find_zero(sdA2B);
 					pA0 = pA1 = A.row(v);
				} else {
					int v0, v1, v2;//v1 is the one that stands out
					v1 = find_different_sign(sdA2B);
					v0 = (v1+1)%3;
					v2 = (v1+2)%3;


					l2l_intersection(O, ld, A.row(v1), A.row(v0) - A.row(v1), pA0);
					l2l_intersection(O, ld, A.row(v2), A.row(v0) - A.row(v2), pA1);
					
				}
				tA0 = point_on_line(pA0, O, ld);
				tA1 = point_on_line(pA1, O, ld);
				if (tA0 > tA1){
					rat temp = tA0;
					tA0 = tA1;
					tA1 = temp;

					auto tempp = pA0;
					pA0 = pA1;
					pA1 = tempp;
				}

			}

 			//Do interval test on B
			rat tB0, tB1;
			RowVector3r pB0, pB1;
			{
				assert(num_sd_B2A0 < 3);
				if (num_sd_B2A0 == 2){
					int v = find_non_zero(sdB2A);
					int v1 = (v+1)%3;
					int v2 = (v+2)%3;
					pB0 = B.row(v1);
					pB1 = B.row(v2);
				} else if (num_sd_B2A0 == 1){
 					int v = find_zero(sdB2A);
 					pB0 = pB1 = B.row(v);
				} else {
					int v0, v1, v2;//v1 is the one that stands out
					v1 = find_different_sign(sdB2A);
					v0 = (v1+1)%3;
					v2 = (v1+2)%3;


					l2l_intersection(O, ld, B.row(v1), B.row(v0) - B.row(v1), pB0);
					l2l_intersection(O, ld, B.row(v2), B.row(v0) - B.row(v2), pB1);
					
				}
				tB0 = point_on_line(pB0, O, ld);
				tB1 = point_on_line(pB1, O, ld);
				if (tB0 > tB1){
					rat temp = tB0;
					tB0 = tB1;
					tB1 = temp;

					auto tempp = pB0;
					pB0 = pB1;
					pB1 = tempp;
				}

			}
			std::cout << "zeros:" << num_sd_B2A0 << " " << num_sd_A2B0 << std::endl; 
			std::cout << pB0 << " " << pB1 << " " << pA0 << " " << pA1 << std::endl; 
			if (tA0 > tB1 || tA1 < tB0){ //no intersection
				return return_v;
			} 
			if (tB0 >= tA0){
				return_v.push_back(pB0);
			} else {
				return_v.push_back(pA0);
			}

			if (tB1 <= tA1){
				return_v.push_back(pB1);
			} else {
				return_v.push_back(pA1);
			}



		}

	}


}
