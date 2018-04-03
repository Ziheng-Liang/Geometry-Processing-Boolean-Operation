#include "constrained_delaunay_triangulation.h"
#include <igl/sort.h>
#include <igl/slice.h>
#include <igl/colon.h>
#include <math.h> 
#include <iostream>

using namespace igl::bol;
using namespace Eigen;
using namespace std;

void igl::bol::constrained_delaunay_triangulation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &C, Eigen::MatrixXi & F) {
	Node* root = new Node();
	MatrixXd projectV;
	cout << "test1" << endl;
	contruct_tree(V, root, projectV);
	cout << "test2" << endl;
	delaunay_triangulation(projectV, root);
	cout << "test3" << endl;
	add_constrained(projectV, C, root);
	cout << "test4" << endl;
	cout << root->polygons.size() << endl;
	F = Eigen::MatrixXi::Zero(root->polygons.size(), 3);
	for (int i = 0; i < F.rows(); i++) {
		F(i, 0) = root->polygons.at(i)->vertex(0);
		F(i, 1) = root->polygons.at(i)->vertex(1);
		F(i, 2) = root->polygons.at(i)->vertex(2);
	}
}

void igl::bol::contruct_tree(const MatrixXd &V, Node* node, MatrixXd &projectV) {
	assert (V.rows() >= 3);
	cout << "test1.0" << endl;
	RowVector3d p1 = V.row(0);
	RowVector3d p2 = V.row(1);
	RowVector3d p3 = V.row(2);
	cout << "test1.1" << endl;
	RowVector3d up = (p2-p1).cross(p3-p1);
	RowVector3d x = (p2-p1).normalized();
	RowVector3d y = up.cross(x).normalized();
	cout << x << endl;
	cout << y << endl;
	cout << "test1.2" << endl;
	MatrixXd Vn = V.rowwise() - p1;
	MatrixXd A = MatrixXd::Zero(2, 3);
	A << x(0), x(1), x(2), y(0), y(1), y(2);
	MatrixXd piA = A.transpose()*(A * A.transpose()).inverse();
	cout << "test1.3" << endl;
	projectV = V * piA;
	MatrixXd temp;
	MatrixXi xindex;

	cout << "test1.4" << endl;
	igl::sort(projectV, 1, true, temp, xindex);
	cout << "test1.5" << endl;
	subdivide(xindex.col(0), node);
	cout << "test1.6" << endl;
}

void igl::bol::subdivide(Eigen::VectorXi index, Node* node) {
	node->index = index;
	int n = index.rows();
	node->size = n;
	if ((index.rows()) <= 3) {
		return;
	}
	Node* left_node = new Node();
	node->left = left_node;
	subdivide(index.head(n/2), left_node);

	Node* right_node = new Node();
	node->right = right_node;
	subdivide(index.tail(n-n/2), right_node);
}

void igl::bol::add_constrained(Eigen::MatrixXd V, Eigen::MatrixXi C, Node* node) {
	for (int i = 0; i < C.rows(); i++) {
		break_polygons(V, C.row(i), node->polygons, node->index);
	}
}

void igl::bol::break_polygons(Eigen::MatrixXd V, Eigen::RowVectorXi C, std::vector<Polygon*> polygons, Eigen::RowVectorXi index) {
	int vidx1 = -1, vidx2 = -1;
	Polygon* start, *next;
	bool intersection = false;
	for (int i = 0; i < polygons.size(); i++) {
		vidx1 = find_vertex(polygons.at(i), C(0));
		if (vidx1 != -1) {
			for (int j = 0; j < polygons.at(i)->size - 1; j++) {
				if (intersect(V.row(C(0)), V.row(C(1)), 
							  V.row(polygons.at(i)->vertex(j)), 
							  V.row(polygons.at(i)->vertex(j+1)))) {
					intersection = true;
					next = polygons.at(i)->adjacent_polygon.at(j);
					break;
				}
			}
			if (intersection) {
				start = polygons.at(i);
				break;
			}
		}		
		vidx2 = find_vertex(polygons.at(i), C(1));
		if (vidx2 != -1) {
			for (int j = 0; j < polygons.at(i)->size - 1; j++) {
				if (intersect(V.row(C(0)), V.row(C(1)), 
							  V.row(polygons.at(i)->vertex(j)), 
							  V.row(polygons.at(i)->vertex(j+1)))) {
					intersection = true;
					next = polygons.at(i)->adjacent_polygon.at(j);
					break;
				}
			}
			if (intersection) {
				start = polygons.at(i);
				break;
			}
		}
	}
	polygons.erase(std::remove(polygons.begin(), polygons.end(), start), polygons.end());
	Polygon* new_polygon;
	while (vidx1 == -1 || vidx2 == -1) {
		merge(start, next, new_polygon);
		polygons.erase(std::remove(polygons.begin(), polygons.end(), next), polygons.end());
		for (int j = 0; j < new_polygon->size - 1; j++) {
			if (intersect(V.row(C(0)), V.row(C(1)), 
						  V.row(new_polygon->vertex(j)), 
						  V.row(new_polygon->vertex(j+1)))) {
				next = new_polygon->adjacent_polygon.at(j);
				start = new_polygon;
				break;
			}
		}
		if (vidx1 == -1) {
			vidx1 = find_vertex(start, C(0));
		}
		if (vidx2 == -1) {
			vidx2 = find_vertex(start, C(1));
		}
	}

	MatrixXi xindex;
	MatrixXd sub_V;
	VectorXi cols;
	igl::colon(0,V.cols()-1, cols);
	igl::slice(V,new_polygon->vertex,cols,sub_V);
	igl::slice(index,new_polygon->vertex,cols,xindex);
	Node* node = new Node();
	subdivide(xindex, node);
	delaunay_triangulation(sub_V, node);
	polygons.insert(polygons.end(), node->polygons.begin(), node->polygons.end());
}


void igl::bol::delaunay_triangulation(MatrixXd V, Node* node) {
	if (node->left) {
		assert(node->right);
	}
	if (node->right) {
		assert(node->left);
	}
	if (!node->left && !node->right) {
		cout << "node->size" << endl;
		cout << node->size << endl;
		if (node->size == 3) {
			Polygon* polygon = new Polygon();
			polygon->vertex = VectorXi::Zero(node->size);
			for (int i = 0; i < node->size; i++) {
				polygon->vertex(i) = node->index(i);
				node->edges.push_back(make_tuple(node->index(i), node->index((i+1) % node->size)));
			}
			polygon->size = node->size;
			polygon->adjacent_polygon = {NULL, NULL, NULL};
			polygon->adjacent_index = {NULL, NULL, NULL};
			node->polygons.push_back(polygon);
		}
		else {
			for (int i = 0; i < node->size; i++) {
				node->edges.push_back(make_tuple(node->index(i), node->index((i+1) % node->size)));
			}
		}
	}
	else {
		cout << "test2.1" << endl;
		cout << "node->right->edges.size()" << endl;
		cout << node->right->edges.size() << endl;
		delaunay_triangulation(V, node->left);
		delaunay_triangulation(V, node->right);
		cout << "node->right->edges.size()" << endl;
		cout << node->right->edges.size() << endl;

		// get partial Vs
		MatrixXd left_sub_V;
		MatrixXd right_sub_V;
		VectorXi cols;
		igl::colon(0,V.cols()-1, cols);
		igl::slice(V,node->left->index,cols,left_sub_V);
		igl::slice(V,node->right->index,cols,right_sub_V);
		cout << "test2.2" << endl;
		cout << left_sub_V << endl;

		// get vertical orders of each sub tree
		MatrixXi left_vert_index;
		MatrixXi right_vert_index;
		MatrixXd temp;
		igl::sort(left_sub_V, 1, true, temp, left_vert_index);
		igl::sort(right_sub_V, 1, true, temp, right_vert_index);
		cout << "test2.25" << endl;

		// find the base LR-edge
		int l_vidx = 0, r_vidx = 0;
		while (l_vidx < left_vert_index.rows() && r_vidx < right_vert_index.rows()) {
			RowVectorXd p1 = left_sub_V.row(left_vert_index(l_vidx,1));
			RowVectorXd p2 = right_sub_V.row(right_vert_index(r_vidx,1));
			bool has_intersection = false;
			cout << "test2.27" << endl;
			for (int i = 0; i < node->left->index.rows(); i++) {
				// need to skip self
				RowVectorXd p3 = V.row(get<0>(node->left->edges.at(i)));
				RowVectorXd p4 = V.row(get<1>(node->left->edges.at(i)));
				cout << "test2.275" << endl;
				if (intersect(p1, p2, p3, p4)) {
					l_vidx ++;
					has_intersection = true;
					break;
				}
			}
			cout << "test2.28" << endl;
			for (int i = 0; i < node->right->index.rows(); i++) {
				RowVectorXd p3 = V.row(get<0>(node->right->edges.at(i)));
				RowVectorXd p4 = V.row(get<1>(node->right->edges.at(i)));
				if (intersect(p1, p2, p3, p4)) {
					r_vidx ++;
					has_intersection = true;
					break;
				}
			}
			if (!has_intersection) {
				break;
			}
		}
		cout << "test2.3" << endl;
		cout << node->index << endl;
		cout << node->left->index << endl;
		cout << node->right->index << endl;

		vector<tuple<int,int>> new_edges;
		new_edges.push_back(make_tuple(node->left->index(l_vidx), node->right->index(r_vidx)));

		int l_ridx = node->left->index(l_vidx);
		int r_ridx = node->right->index(r_vidx);
		cout << "test2.4" << endl;
		cout << node->left->index(l_vidx) << endl;
		cout << node->right->index(r_vidx) << endl;
		// if no more edges are found, triangulation is complete
		int l_candidate_final;
		int r_candidate_final;
		do {
			cout << "loooooooooop" << endl;
			// find candidates from both side
			vector<tuple<int,int>> r_candidate;
			vector<tuple<int,int>> l_candidate;
			for (int i = 0; i < node->right->edges.size(); i++) {
				int other;

				// find edges contains right node of base LR edge
				if (get<0>(node->right->edges.at(i)) == r_ridx){
					other = get<1>(node->right->edges.at(i));
				}
				else if (get<1>(node->right->edges.at(i)) == r_ridx) {
					other = get<0>(node->right->edges.at(i));
				}
				else {
					continue;
				}
				// to do: need to check angle clock vs counter-clock
				double a = angle(V.row(r_ridx), V.row(l_ridx), V.row(other));
				cout << "angle: " << endl;
				cout << r_ridx << endl;
				cout << other << endl;
				cout << l_ridx << endl;
				cout << a << endl;
				if (a > 3.1415926535) {
					break;
				}
				if (a < 0) {
					break;
				}
				cout << "angle:" << endl;
				cout << a << endl;
				// add it to a specific location based on angle
				bool added = false;
				for (int j = 0; j < r_candidate.size(); j++) {
					if (other == get<0>(r_candidate.at(j))) {
						added = true;
					}
					else if (a < get<1>(r_candidate.at(j))){
						r_candidate.insert(r_candidate.begin() + j, make_tuple(other, a));
						cout << "added to r_candidate1" << endl;
						added = true;
					}
					// add it to the end if not added before
				}
				if (!added) {
					cout << "added to r_candidate2" << endl;
					r_candidate.push_back(make_tuple(other, a));
				}
			}

			cout << "test2.5" << endl;

			for (int i = 0; i < node->left->edges.size(); i++) {
				int other;

				// find edges contains left node of base LR edge
				if (get<0>(node->left->edges.at(i)) == l_ridx){
					other = get<1>(node->left->edges.at(i));
				}
				else if (get<1>(node->left->edges.at(i)) == l_ridx) {
					other = get<0>(node->left->edges.at(i));
				}
				else {
					continue;
				}
				// to do: need to check angle clock vs counter-clock
				double a = angle(V.row(l_ridx), V.row(other), V.row(r_ridx));
				cout << "angle: " << endl;
				cout << r_ridx << endl;
				cout << other << endl;
				cout << l_ridx << endl;
				cout << a << endl;
				if (a > 3.1415926535) {
					break;
				}
				if (a < 0) {
					break;
				}
				bool added = false;
				// add it to a specific location based on angle
				for (int j = 0; j < l_candidate.size(); j++) {
					if (other == get<0>(l_candidate.at(j))) {
						added = true;
					}
					else if (a < get<1>(l_candidate.at(j))){
						l_candidate.insert(l_candidate.begin() + j, make_tuple(other, a));
						added = true;
					}
				}
				if (!added) {
					l_candidate.push_back(make_tuple(other, a));
				}
			}
			// find final candidate from both side
			cout << "test2.51" << endl;
			RowVectorXd p_left = V.row(l_ridx);
			RowVectorXd p_right = V.row(r_ridx);
			RowVectorXd center;
			cout << r_candidate.size() << endl;
			l_candidate_final = -1;
			r_candidate_final = -1;
			for (int i = 0; i < r_candidate.size(); i++) {
				cout << "test2.515" << endl;
				if ((i + 1) == r_candidate.size()) {
					r_candidate_final = get<0>(r_candidate.at(i));
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(r_candidate.at(i))), center);
				// if next candidate is outside, this is the final candidate
				if ((p_left-center).dot(p_left-center) < (V.row(get<0>(r_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(r_candidate.at(i+1))))) {
					r_candidate_final = get<0>(r_candidate.at(i));
					break;
				}
				// if not, remove this edge
				else {
					cout << "test2.52" << endl;
					remove_edge_from_node(node->right, r_ridx, get<0>(r_candidate.at(i)));
					cout << "test2.53" << endl;
				}
			}

			cout << "test2.6" << endl;
			cout << l_candidate.size() << endl;
			if (l_candidate.size() == 1) {
				l_candidate_final = get<0>(l_candidate.at(0));
			}
			for (int i = 0; i < l_candidate.size(); i++) {
				cout << "test2.61" << endl;
				if ((i + 1) == l_candidate.size()) {
					l_candidate_final = get<0>(l_candidate.at(i));
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(l_candidate.at(i))), center);
				if ((p_left-center).dot(p_left-center) < (V.row(get<0>(l_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(l_candidate.at(i+1))))) {
					cout << "test2.611" << endl;
					l_candidate_final = get<0>(l_candidate.at(i));
					break;
				}
				else {
					cout << "test2.62" << endl;
					remove_edge_from_node(node->left, l_ridx, get<0>(l_candidate.at(i)));
					cout << "test2.63" << endl;
				}
			}

			// final step
			cout << "test2.7" << endl;
			cout << l_candidate_final << endl;
			cout << r_candidate_final << endl;
			if (l_candidate_final == -1 && r_candidate_final != -1) {
				r_ridx = r_candidate_final;
			}
			else if (l_candidate_final != -1 && r_candidate_final == -1) {
				l_ridx = l_candidate_final;
			}
			else if (l_candidate_final != -1 && r_candidate_final != -1) {
				get_circle_center(p_left, p_right, V.row(l_candidate_final), center);
				if ((p_left-center).dot(p_left-center) < (V.row(r_candidate_final) - center)
					.dot(V.row(r_candidate_final) - center)) {
					l_ridx = l_candidate_final;
				}
				else {
					r_ridx = r_candidate_final;
				}
			}
			else {
				break;
			}
			Polygon* new_polygon =  new Polygon();

			cout << "test2.75" << endl;
			new_polygon->size = 3;
			cout << "test2.76" << endl;
			new_polygon->vertex = VectorXi::Zero(3);
			new_polygon->vertex(0) = r_ridx;
			cout << "test2.77" << endl;
			new_polygon->adjacent_polygon.push_back(NULL);
			cout << "test2.78" << endl;
			new_polygon->vertex(1) = l_ridx;
			cout << "test2.8" << endl;
			if (l_ridx == get<0>(new_edges.at(new_edges.size()-1))) {
				new_polygon->vertex(2) = get<1>(new_edges.at(new_edges.size()-1));
				for (int i = 0; i < node->right->polygons.size();i++) {
					int oppo_idx = exist_edges(node->right->polygons.at(i), 
											   new_polygon->vertex(2), 
											   new_polygon->vertex(0));
					if (oppo_idx != -1) {
						new_polygon->adjacent_polygon.push_back(node->polygons.at(node->polygons.size()-1));
						new_polygon->adjacent_polygon.push_back(node->right->polygons.at(i));
						new_polygon->adjacent_index.push_back(NULL);
						new_polygon->adjacent_index.push_back(exist_edges(new_polygon->adjacent_polygon.at(1), 
																		  new_polygon->vertex(1), 
																		  new_polygon->vertex(2)));
						new_polygon->adjacent_index.push_back(oppo_idx);	
						new_polygon->adjacent_polygon.at(1)->adjacent_index.at(new_polygon->adjacent_index.at(1)) = 1;
						new_polygon->adjacent_polygon.at(2)->adjacent_index.at(new_polygon->adjacent_index.at(2)) = 1;
						break;
					}
				}
			}
			else {
				cout << "test2.9" << endl;
				new_polygon->vertex(2) = get<0>(new_edges.at(new_edges.size()-1));
				for (int i = 0; i < node->right->polygons.size();i++) {
					int oppo_idx = exist_edges(node->left->polygons.at(i), 
											   new_polygon->vertex(1), 
											   new_polygon->vertex(2));
					if (oppo_idx != -1) {
						new_polygon->adjacent_polygon.push_back(node->left->polygons.at(i));
						new_polygon->adjacent_polygon.push_back(node->polygons.at(node->polygons.size()-1));
						new_polygon->adjacent_index.push_back(NULL);
						new_polygon->adjacent_index.push_back(oppo_idx);
						new_polygon->adjacent_index.push_back(exist_edges(new_polygon->adjacent_polygon.at(2), 
																		  new_polygon->vertex(2), 
																		  new_polygon->vertex(0)));
						new_polygon->adjacent_polygon.at(1)->adjacent_index.at(new_polygon->adjacent_index.at(1)) = 1;
						new_polygon->adjacent_polygon.at(2)->adjacent_index.at(new_polygon->adjacent_index.at(2)) = 1;
						break;
					}
				}
			}
			new_edges.push_back(make_tuple(l_ridx, r_ridx));
			node->polygons.push_back(new_polygon);
			cout << new_edges.size() << endl;
			cout << l_ridx << endl;
			cout << r_ridx << endl;
			cout << "end loooooooooop" << endl;
		} while ((l_candidate_final != -1) || (r_candidate_final != -1));
		// add all new edges into the list
		node->edges.insert(node->edges.end(), node->left->edges.begin(), node->left->edges.end());
		node->edges.insert(node->edges.end(), node->right->edges.begin(), node->right->edges.end());
		node->edges.insert(node->edges.end(), new_edges.begin(), new_edges.end());
	}

}


// need to determine clockwise angle or counter-clockwise angle
// double igl::bol::angle(RowVectorXd p1, RowVectorXd p2, RowVectorXd p3) {
// 	// cout << "inside angle" << endl;
// 	// cout << p1 << endl;
// 	// cout << p2 << endl;
// 	// cout << p3 << endl;
// 	// cout << (p2 - p1).dot(p3 - p1) / ((p2 - p1).norm() * (p3 - p1).norm()) << endl;
// 	return acos((p2 - p1).dot(p3 - p1) / ((p2 - p1).norm() * (p3 - p1).norm()));

// }

double igl::bol::angle(RowVectorXd p1, RowVectorXd p2, RowVectorXd p3) {
	double Dir_C_to_A = atan2(p2(1) - p1(1), p2(0)- p1(0));
	double Dir_C_to_B = atan2(p3(1) - p1(1), p3(0)- p1(0));
	double Angle_ACB = Dir_C_to_A - Dir_C_to_B;

	// Handle wrap around
	const double Pi = acos(-1);  // or use some π constant
	if (Angle_ACB > Pi) Angle_ACB -= 2*Pi;
	else if (Angle_ACB < -Pi) Angle_ACB += 2*Pi;

	// Answer is in the range of [-pi...pi]
	return Angle_ACB;
}
// bool igl::bol::intersect(RowVector3d p1, RowVector3d p2, RowVector3d p3, RowVector3d p4) {
// 	RowVector3d v1 = p2 - p1;
// 	RowVector3d v2 = p4 - p3;
// 	return (v1.cross(v2)).dot(p3 - p1) == 0;
// }

// bool igl::bol::intersect(RowVector3d p1, RowVector3d p2, RowVector3d p3, RowVector3d p4) {
// 	RowVector3d v1 = p2 - p1;
// 	RowVector3d v2 = p4 - p3;
// 	return (v1.cross(v2)).dot(p3 - p1) == 0;
// }

bool igl::bol::intersect(const RowVectorXd &p1, const RowVectorXd &p2, const RowVectorXd &p3, const RowVectorXd &p4) {
	assert(p1.cols() == 2);
	if (p1==p3 || p1==p4 || p2==p3 || p2==p4) {
		return false;
	}
	RowVectorXd v1 = p2 - p1;
	RowVectorXd v2 = p4 - p3;
	RowVectorXd d = p3 - p1;
	// colinear
	double det = v1(0) * v2(1) - v1(1) * v2(0);
	if (det < 0.0000001) {
		return false;
	}

	double r = (d(0) * v2(1) - d(1) * v2(0)) / det;
    double s = (v1(0) * d(1) - v1(0) * d(0)) / det;

	return !(r < 0 || r > 1 || s < 0 || s > 1);
}

void igl::bol::get_circle_center(const RowVectorXd &p1, const RowVectorXd &p2, const RowVectorXd &p3, RowVectorXd &center) {
	assert(p1.cols() == 2);
	double ma = (p2(1) - p1(1)) / (p2(0) - p1(0));
	double mb = (p3(1) - p2(1)) / (p3(0) - p2(0));
	double x = (ma*mb*(p1(1)-p3(1)) + mb*(p1(0)+p2(0) - ma*(p2(0)+p3(0))))/(2*(mb-ma));
	double y = -1/ma * (x - (p1(0) + p2(0))/2) + (p1(1)+p2(1))/2;
	center = VectorXd::Zero(2);
	center << x, y;
}

// void igl::bol::get_circle_center(RowVector3d p1, RowVector3d p2, RowVector3d p3, RowVector3d center) {
// 	RowVector3d a = (p1 + p2)/2;
// 	RowVector3d b = (p2 + p3)/2;
// 	RowVector3d up = (p1-p2).cross(p1-p3);
// 	RowVector3d c = up.cross(p1-p2);
// 	RowVector3d d = up.cross(p2-p3);
// 	double m = (dmnop(a,c,d,c)*dmnop(d,c,b,a)-dmnop(a,c,b,c)*dmnop(d,c,d,c)) / 
// 			   (dmnop(b,a,b,a)*dmnop(d,c,d,c)-dmnop(d,c,b,a)*dmnop(d,c,b,a));
// 	center = p1 + m*(p2-p1);
// }

// double igl::bol::dmnop(RowVector3d m, RowVector3d n, RowVector3d o, RowVector3d p){
// 	//replace with rational later
// 	return (m - n).dot(o - p);
// }