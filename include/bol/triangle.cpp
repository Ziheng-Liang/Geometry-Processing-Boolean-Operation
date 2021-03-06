#include "triangle.h"
#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <unordered_set>
namespace igl{
namespace bol{


inline bool is_degenerate(Matrix33r V){
	RowVector3r a  = V.row(1) - V.row(0);
	RowVector3r b = V.row(2) - V.row(0);

	return (a.cross(b).array() == 0).all();
}


static bool vectors_equal(const RowVector3r & a, const RowVector3r & b){
	return (a.array() == b.array()).all();
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
	const RowVector3r & n2, const rat & d2, RowVector3r & p){
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


//Check whether a point lies between two point.
//Precondition, 3 points are colinear. b != a
//Using a bounding box technique to check, much efficient.
//return if it is inside
static bool p_lies_ls(const RowVector3r &p, const RowVector3r &a, const RowVector3r & b){
	rat min_coord[3];
	min_coord[0] = a(0,0);
	min_coord[1] = a(0,1);
	min_coord[2] = a(0,2);

	rat max_coord[3];
	max_coord[0] = a(0,0);
	max_coord[1] = a(0,1);
	max_coord[2] = a(0,2);

	for (int i = 0; i < 3; i++){
		if (b(0, i) < min_coord[i]) min_coord[i] = b(0, i);
		if (b(0, i) > max_coord[i]) max_coord[i] = b(0, i);
	}

	for (int i = 0; i < 3; i++){
		if (p(0, i) < min_coord[i] || p(0,i) > max_coord[i]) return false; 
	}
	return true;

}

//Find the itersection of two lines
//Precondition: two lines coplanar
static void l2l_intersection(const RowVector3r &x1, const RowVector3r &d1, 
	const RowVector3r &x2, const RowVector3r &d2, RowVector3r & p){
	//Line 1: = x1 + t * d1;
	//Line 2: = x2 + s * d2
	//Solve this linear equation (over determined)
	//We know its on the same plane before hand, therefore we just need two axis

	for (int i = 0; i < 3; i++){
		int a0 = i;
		int a1 = (i+1)%3;
		int a2 = (i+2)%3;
		if (d1(a0) * (-d2(a1)) - d1(a1) * (-d2(a0)) != 0) { //if it is invertible
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

			auto ts = left.inverse() * right;
			auto t = ts(0,0);
			p = x1 + t * d1;

		}


	}

}

static bool points_colinear(const RowVector3r & a, const RowVector3r & b, const RowVector3r & c){
	RowVector3r d1 = a - b;
	RowVector3r d2 = c - b;
	RowVector3r crossprod =  d1.cross(d2);
	return (crossprod.array() == 0).all();
}


//Find the intersection of two line segment a0-a1 and b0-b1
//Precondition: two line segment coplanar, a0, a1 not the same point, b0 b1 not the same point
inline std::vector<RowVector3r> ls2ls_intersection(const RowVector3r &a0, const RowVector3r &a1, 
			const RowVector3r &b0, const RowVector3r &b1){
	assert(!(a0.array() == a1.array()).all());
	assert(!(b0.array() == b1.array()).all());
	
	std::vector<RowVector3r> return_v;
	//two lines can be paralle
	
	RowVector3r da = a1 - a0;
	RowVector3r db = b1 - b0;

	RowVector3r cdadb = da.cross(db);
	if ((cdadb.array() == 0).all() ){ //da cross db = 0, da db are paralell
		RowVector3r dc = a1 - b0;
		RowVector3r cdcdb = dc.cross(db);

		if ((cdcdb.array() == 0).all()){ //if they are collinear
			//We have 4 case
			bool b0_ina = p_lies_ls(b0, a0, a1);
			bool b1_ina = p_lies_ls(b1, a0, a1);
			//case 1: one point in b lies in a0-a1. intersection is b-a
			//case 2: two point in b lies in a0-a1, intersection is b0b1
			//case 3: 0 point in b lies in a0-a1, 
				//subcase 1: 0 point of a lies in b0-b1, no intersection
				//subcase 2: 2 points of a lies in b0-b1, intersection is a0a1

			if (b0_ina && !b1_ina){
				rat dotb1b0 = (a1 - a0).dot(b1 - a0);
				if (dotb1b0 > 0){		//a0-b0-a1-b1//
					return_v.push_back(b0);
					return_v.push_back(a1);
					return return_v;
				} else if (dotb1b0 < 0){ // b1-a0-b0-a1//
					return_v.push_back(a0);
					return_v.push_back(b0);
					return return_v;
				}

			} else if (!b0_ina && b1_ina){
				rat dotb1b0 = (a1-a0).dot(b0 - a0);
				if (dotb1b0 > 0){ //a0-b1-a1-b0
					return_v.push_back(b1);
					return_v.push_back(a1);
					return return_v;
				} else if (dotb1b0 < 0){ //b0-a0-b1-a1
					return_v.push_back(a0);
					return_v.push_back(b1);
					return return_v;
				}

			} else if (b0_ina && b1_ina){ // intersection is b0b1
				return_v.push_back(b0);
				return_v.push_back(b1);
				return return_v;
			} else if (!b0_ina && !b1_ina){
				bool a0_inb = p_lies_ls(a0, b0, b1);
				if (a0_inb){ //intersection is a0a1
					return_v.push_back(a0);
					return_v.push_back(a1);
					return return_v;
				} else { // no intersection
					return return_v;
				}
			} 
		} else { //parallel but not colinear, then they should not intersect
			return return_v;
		}
	} else { // if they are not parallel. find the intersection, and check if the point is inside of two line segment
		RowVector3r intersect_point; 
		l2l_intersection(a0, da, b0, db, intersect_point);
		bool p_ina = p_lies_ls(intersect_point, a0, a1);
		bool p_inb = p_lies_ls(intersect_point, b0, b1);
		if (p_ina && p_inb){
			return_v.push_back(intersect_point);
			return return_v;
		} else {
			return return_v;
		}
	}
	//Here should have capture all the cases.
	assert(false);
}



//Find the t value that cooresponds a poing on line.
//Precondition: p is on the line, dicretion d is non zero
static rat point_on_line(const RowVector3r & p, const RowVector3r & x, const RowVector3r &d){
	RowVector3r v = p - x;
	for (int i = 0; i < 3; i++){
		if (d(0,i) != 0) return v(0,i)/d(0,i);
	}
}

//Check whether point is inside the triangle
static bool point_in_triangle(const RowVector3r & q, const Matrix33r & A){
	RowVector3r a = A.row(0);
	RowVector3r b = A.row(1);
	RowVector3r c = A.row(2);
    auto u = b - a;
    auto v = c - a;
    auto n = u.cross(v);
    auto w = q - a;
    rat a0 = u.cross(w).dot(n)/n.dot(n);
    rat a1 = w.cross(v).dot(n)/n.dot(n);
    rat a2 = 1 - a0 - a1;
    return (a0 >= 0 && a0 <=1) && (a1 >= 0 && a1 <=1)  &&(a2 >= 0 && a2 <=1);
}




//Return the list of segments
static std::vector<RowVector3r> coplanar_t2t_intersection(const Matrix33r & A, const Matrix33r & B, bool subdivide_edges){

	std::vector<RowVector3r> return_v;
	//Array of bool to check if a point is in triangle or not
	bool vint[2][3];
	bool no_point_inside = true;
	for (int t = 0; t < 2; t++){
		const Matrix33r & cur_triangle = (t == 0) ? A : B;
		const Matrix33r & other_triangle = (t == 0) ? B : A;
		for (int i = 0; i < 3; i++){
			RowVector3r cur_vertex = cur_triangle.row(i);
			vint[t][i] = point_in_triangle(cur_vertex, other_triangle);
			if (vint[t][i]) no_point_inside = false;
		}
	}
	std::array<std::array<std::vector<RowVector3r>, 3>, 3> intersect_points;
	bool no_intersection = true;
	for (int i = 0; i < 3; i++){
		RowVector3r a0 = A.row(i);
		RowVector3r a1 = A.row((i+1)%3);
		for (int j = 0; j < 3; j++){
			RowVector3r b0 = B.row(j);
			RowVector3r b1 = B.row((j+1)%3);
			auto intersections = ls2ls_intersection(a0, a1, b0, b1);
			// for (int k = 0; k < intersections.size(); k ++){
			// }
			if (intersections.size() != 0) no_intersection = false;
			intersect_points[i][j] = intersections;
		}

	}

	if (no_intersection && no_point_inside){ //No overlapping area
		return return_v;
	}

	//Add A's intersection loop through [0, 1] [1, 2] [2, 3]
	//Here it could add potentially two duplicate point
	for (int t = 0; t < 2; t++){
		const Matrix33r & cur_triangle = (t == 0) ? A : B;

		for (int i = 0; i < 3; i++){
			if (vint[t][i] && vint[t][(i+1)%3]){ //If both vertex is inside
				return_v.push_back(cur_triangle.row(i));
				return_v.push_back(cur_triangle.row((i+1)%3));
			} else if (vint[t][i]){ //if one vertex is inside, the other is strictly outside
				return_v.push_back(cur_triangle.row(i));
				bool vertex_added = false;
				for (int j = 0; j < 3; j++){
					int m, n;
					m = (t == 0) ? i : j;
					n = (t == 0) ? j : i;
					if (intersect_points[m][n].size()==1 && !vectors_equal(intersect_points[m][n][0], cur_triangle.row(i))){
						return_v.push_back(intersect_points[m][n][0]);
						vertex_added = true;
						break;
					}
				}
				if (!vertex_added){
					return_v.push_back(cur_triangle.row(i));
				}
				if (subdivide_edges){
					return_v.push_back(return_v[return_v.size() - 1]);
					return_v.push_back(cur_triangle.row((i+1)%3));
				}
			} else if (vint[t][(i+1)%3]){ // if other vertex is inside
				return_v.push_back(cur_triangle.row((i+1)%3));
				bool vertex_added = false;
				for (int j = 0; j < 3; j++){
					int m, n;
					m = (t == 0) ? i : j;
					n = (t == 0) ? j : i;
					if (intersect_points[m][n].size()==1 && !vectors_equal(intersect_points[m][n][0], cur_triangle.row((i+1)%3))){
						return_v.push_back(intersect_points[m][n][0]);
						vertex_added = true;
						break;
					}
				}
				if (!vertex_added){
					return_v.push_back(cur_triangle.row((i+1)%3));
				}
				if (subdivide_edges){
					return_v.push_back(return_v[return_v.size() - 1]);
					return_v.push_back(cur_triangle.row(i));
				}
			} else {
				int num_intersection = 0;
				for(int j = 0; j < 3; j++){
					int m, n;
					m = (t == 0) ? i : j;
					n = (t == 0) ? j : i;
					if (intersect_points[m][n].size()==1){ 
						return_v.push_back(intersect_points[m][n][0]);
						num_intersection++;
					}
				}
				if (subdivide_edges){
					if (num_intersection == 0){
						return_v.push_back(cur_triangle.row(i));
						return_v.push_back(cur_triangle.row((i+1)%3));
					} else {
						auto first = return_v[return_v.size()-1]; 
						auto second = return_v[return_v.size()-2];
						if (p_lies_ls(first, second, cur_triangle.row(i))){ // v--first--second--v
							return_v.push_back(cur_triangle.row(i));
							return_v.push_back(first);
							return_v.push_back(second);
							return_v.push_back(cur_triangle.row((i+1)%3));
						} else {											// v--second--first--v
							return_v.push_back(cur_triangle.row(i));
							return_v.push_back(second);
							return_v.push_back(first);
							return_v.push_back(cur_triangle.row((i+1)%3));
						}
					}
				}
			}
		}
	}


	return return_v;
}

static std::vector<RowVector3r> coplanar_t2t_intersection(const Matrix33r & A, const Matrix33r & B){
	return coplanar_t2t_intersection(A, B, false);
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
//subdived edges: try to subdive edges along the way. only apply to A if not coplanar
inline std::vector<RowVector3r> t2t_intersect(const Matrix33r & A, const Matrix33r & B, bool subdivide_edges){
	assert(!is_degenerate(A));
	assert(!is_degenerate(B));

	RowVector3r nA = (A.row(1) - A.row(0)).cross(A.row(2) - A.row(0));
	rat dA = -nA.dot(A.row(0));
	rat sdB2A[3];

	for (int i = 0; i < 3; i++){
		sdB2A[i] = nA.dot(B.row(i)) + dA;
	}


	RowVector3r nB = (B.row(1) - B.row(0)).cross(B.row(2) - B.row(0));
	rat dB = -nB.dot(B.row(0));
	rat sdA2B[3];
	for (int i = 0; i < 3; i++){
		sdA2B[i] = nB.dot(A.row(i)) + dB;
	}


	std::vector<RowVector3r> return_v;

	if (sdB2A[0] == 0 && sdB2A[1] == 0 && sdB2A[2] == 0){ //coplanar
		auto result = coplanar_t2t_intersection(A, B, subdivide_edges);
		return result;
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

			//Do interval test on A0
			rat tA0, tA1;
			RowVector3r pA0, pA1;
			t2plane_intersect_line(A, sdA2B, O, ld, pA0, pA1, tA0, tA1);

 			//Do interval test on B
			rat tB0, tB1;
			RowVector3r pB0, pB1;
			t2plane_intersect_line(B, sdB2A, O, ld, pB0, pB1, tB0, tB1);
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

			if (subdivide_edges){
				auto first = return_v[0];
				auto second = return_v[1];
				for (int i = 0; i < 3; i++){
					bool first_colinear = points_colinear(first, A.row(i), A.row((i+1)%3));
					bool second_colinear = points_colinear(second, A.row(i), A.row((i+1)%3));
					if (!first_colinear && ! second_colinear){
						return_v.push_back(A.row(i));
						return_v.push_back(A.row((i+1)%3));
					} else if (first_colinear && second_colinear){
						if (p_lies_ls(first,  A.row(i), second)){ // v--first--second--v
							return_v.push_back(A.row(i));
							return_v.push_back(first);
							return_v.push_back(second);
							return_v.push_back(A.row((i+1)%3));
						} else {
							return_v.push_back(A.row(i));
							return_v.push_back(second);
							return_v.push_back(first);
							return_v.push_back(A.row((i+1)%3));
						}
					} else if (first_colinear && !second_colinear){
						return_v.push_back(A.row(i));
						return_v.push_back(first);
						return_v.push_back(first);
						return_v.push_back(A.row((i+1)%3));
					} else if (second_colinear && !first_colinear){
						return_v.push_back(A.row(i));
						return_v.push_back(second);
						return_v.push_back(second);
						return_v.push_back(A.row((i+1)%3));
					}
				}

			}
		}

	}
	return return_v;

}


inline std::vector<RowVector3r> t2t_intersect(const Matrix33r & A, const Matrix33r & B){
	return t2t_intersect(A, B, false);
}



inline Eigen::RowVector3d rat_to_double(RowVector3r &to_cast){
	Eigen::RowVector3d r;
	for (int i = 0; i < to_cast.size(); i++){
		r(i) = static_cast<double>(to_cast(i));
	}
	return r;
}

static int find_first_rowvector(const std::vector<RowVector3r> & v, const RowVector3r & q){
	for (int i = 0; i < v.size(); i++){
		if (vectors_equal(v[i], q)){
			return i;
		}
	}
	return -1;
}


inline void t2t_intersect_on_A(const Matrix33r & A, const Matrix33r & B, MatrixXr & AV, Eigen::MatrixXi & AF){
	std::vector<RowVector3r> list_of_ls = t2t_intersect(A, B, true);
	std::vector<RowVector3r> list_of_v;

	//insert all vertices in A
	list_of_v.push_back(A.row(0));
	list_of_v.push_back(A.row(1));
	list_of_v.push_back(A.row(2));

	auto hash=[](const std::pair<int, int>& p){
	    return p.first + p.second;
	};

	auto equal=[](const std::pair<int, int> & p1, const std::pair<int, int> & p2){
	    return (p1.first == p2.first && p1.second == p2.second) || (p1.first == p2.second && p1.second == p2.first);
	};

	std::unordered_set<std::pair<int, int>, decltype(hash), decltype(equal)> set_of_c(0, hash, equal);

	for (int i = 0; i < list_of_ls.size()/2; i++){
		auto first_v = list_of_ls[i*2];
		auto second_v = list_of_ls[i*2+1];

		int first_index = find_first_rowvector(list_of_v, first_v);
		int second_index = find_first_rowvector(list_of_v, second_v);

		if (first_index == -1) {
			list_of_v.push_back(first_v);
			first_index = list_of_v.size() - 1;
		}

		if (second_index == -1) {
			list_of_v.push_back(second_v);
			second_index = list_of_v.size() - 1;
		}
		if (vectors_equal(first_v, second_v)){ // if two points are the same, we dont need constraints
		} else {//append this relation to set_of_c
			set_of_c.insert(std::make_pair(first_index, second_index));
		}
	}



	AV.resize(list_of_v.size(),3);
	AF.resize(set_of_c.size(),2);
	for(int i = 0; i < list_of_v.size(); i++){
		AV.row(i) = list_of_v[i];
	}

	int index = 0;
	for (auto&& c_pair :set_of_c){
		AF(index, 0) = c_pair.first;
		AF(index, 1) = c_pair.second;
		index++;
	}

}

}}

