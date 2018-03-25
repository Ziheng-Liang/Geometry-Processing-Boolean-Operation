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


//Find the itersection of two lines
//Precondition: two lines coplanar
static void l2l_intersection(const RowVector3r &x1, const RowVector3r &d1, 
	const RowVector3r &x2, const RowVector3r &d2, RowVector3r & p){
	//Line 1: = x1 + t * d1;
	//Line 2: = x2 + s * d2
	//Solve this linear equation (over determined)
	//We know its on the same plane before hand, therefore we just need two axis
	std::cout << "l1: " << d1 << " x1: " << x1 << std::endl;
	std::cout << "l2: " << d2 << " x2: " << x2 << std::endl; 

	for (int i = 0; i < 3; i++){
		int a0 = i;
		int a1 = (i+1)%3;
		int a2 = (i+2)%3;
		if (d1(a0) * (-d2(a1)) - d1(a0) * (-d2(a1)) != 0) { //if it is invertible
			MatrixXr left;
			left.resize(2,2);
			left(0,0) = d1(0,a0);
			left(0,1) = - d2(0,a0);
			left(1,0) = d1(0,a1);
			left(1,1) = - d2(0,a1);
			VectorXr right;
			right.resize(2,1);
			right(0,0) = x2(0,a0) - x1(0,a0);
			right(1,0) = x2(0,a1) - x1(0,a1);
			std::cout << "KK" << std::endl;
			auto ts = left.inverse() * right;
			auto t = ts(0,0);
			p = x1 + t * d1;
			std::cout << "K" << std::endl;

		}


	}

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

//A helper: Find the intersection of a triangle touches a plane and a line
static void t2plane_intersect_line(const Matrix33r & A, rat sd[3], const RowVector3r & p, const RowVector3r & d, RowVector3r & p0, RowVector3r &p1, rat & t0, rat & t1){
	
	int num_sd_0 = 0;
	for (int i = 0; i< 3; i++){
		if (sd[i] == 0) num_sd_0++;
	}
	assert(num_sd_0 < 3);


	if (num_sd_0 == 2){
		int v = find_non_zero(sd);
		int v1 = (v+1)%3;
		int v2 = (v+2)%3;
		p0 = A.row(v1);
		p1 = A.row(v2);
	} else if (num_sd_0 == 1){
		int v0 = find_zero(sd);
		int v1 = (v0+1)%3;
		int v2 = (v0+2)%3;
		if (sd[v1] * sd[v2] < 0) {//different sign, on different sides
			l2l_intersection(p, d, A.row(v0), A.row(v1) - A.row(v0), p0);
			l2l_intersection(p, d, A.row(v2), A.row(v1) - A.row(v2), p1);
		} else { // All on the same side. no need to test
			p0 = p1 = A.row(v0);	
		}
	} else {
		int v0, v1, v2;//v0 is the one that stands out, v1 v2 on opposite side
		v0 = find_different_sign(sd);
		v1 = (v0+1)%3;
		v2 = (v0+2)%3;
		l2l_intersection(p, d, A.row(v1), A.row(v0) - A.row(v1), p0);
		l2l_intersection(p, d, A.row(v2), A.row(v0) - A.row(v2), p1);
		
	}
	t0 = point_on_line(p0, p, d);
	t1 = point_on_line(p1, p, d);
	if (t0 > t1){
		rat temp = t0;
		t0 = t1;
		t1 = temp;

		auto tempp = p0;
		p0 = p1;
		p1 = tempp;
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

	for (int i = 0; i < 3; i++){
		sdB2A[i] = nA.dot(B.row(i)) + dA;
	}


	RowVector3r nB = (B.row(1) - B.row(0)).cross(B.row(2) - B.row(0));
	rat dB = -nB.dot(A.row(0));
	rat sdA2B[3];
	for (int i = 0; i < 3; i++){
		sdA2B[i] = nB.dot(A.row(i)) + dB;
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

			std::cout << "Here" << ld << std::endl;
			//Do interval test on A0
			rat tA0, tA1;
			RowVector3r pA0, pA1;
			t2plane_intersect_line(A, sdA2B, O, ld, pA0, pA1, tA0, tA1);

 			//Do interval test on B
			rat tB0, tB1;
			RowVector3r pB0, pB1;
			t2plane_intersect_line(B, sdB2A, O, ld, pB0, pB1, tB0, tB1);
			// std::cout << "zeros:" << num_sd_B2A0 << " " << num_sd_A2B0 << std::endl; 
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
