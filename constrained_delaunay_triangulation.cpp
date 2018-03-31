#include "constrained_delaunay_triangulation.h"
#include <igl/sort.h>
#include <igl/slice.h>
#include <igl/colon.h>
#include <math.h> 
#include <iostream>

using namespace igl::bol;
using namespace Eigen;
using namespace std;

void igl::bol::constrained_delaunay_triangulation(Eigen::MatrixXd V, Eigen::MatrixXi C, Eigen::MatrixXi F) {
	Node* root = new Node();
	MatrixXd projectV;
	cout << "test1" << endl;
	contruct_tree(V, root, projectV);
	cout << "test2" << endl;
	delaunay_triangulation(projectV, root);
	cout << "test3" << endl;
	add_constrained(projectV, C, root);
	cout << "test4" << endl;
	F = Eigen::MatrixXi::Zero(root->polygons.size(), 3);
	for (int i = 0; i < F.rows(); i++) {
		F(i, 0) = root->polygons.at(i)->vertex(0);
		F(i, 1) = root->polygons.at(i)->vertex(1);
		F(i, 2) = root->polygons.at(i)->vertex(2);
	}
}

void igl::bol::contruct_tree(MatrixXd V, Node* node, MatrixXd projectV) {
	using namespace std;
	using namespace Eigen;
	assert (V.rows() >= 3);
	RowVector3d p1 = V.row(0);
	RowVector3d p2 = V.row(1);
	RowVector3d p3 = V.row(2);
	RowVector3d up = (p2-p1).cross(p3-p1);
	RowVector3d x = (p2-p2).normalized();
	RowVector3d y = x.cross(up).normalized();
	MatrixXd Vn = V.rowwise() - p1;
	MatrixXd A = MatrixXd::Zero(2, 3);
	A << x(0), x(1), x(2), y(0), y(1), y(2);
	MatrixXd piA = A.transpose()*(A * A.transpose()).inverse();
	projectV = V * piA;
	MatrixXd temp;
	MatrixXi xindex;
	igl::sort(projectV, 1, true, temp, xindex);
	subdivide(xindex, node);
}

void igl::bol::subdivide(Eigen::RowVectorXi index, Node* node) {
	if ((index.rows()) <= 3) {
		return;
	}
	node->index = index;
	int n = index.rows();
	node->size = n;
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
		delaunay_triangulation(V, node->left);
		delaunay_triangulation(V, node->right);

		// get partial Vs
		MatrixXd left_sub_V;
		MatrixXd right_sub_V;
		VectorXi cols;
		igl::colon(0,V.cols()-1, cols);
		igl::slice(V,node->left->index,cols,left_sub_V);
		igl::slice(V,node->right->index,cols,right_sub_V);

		// get vertical orders of each sub tree
		MatrixXi left_vert_index;
		MatrixXi right_vert_index;
		MatrixXd temp;
		igl::sort(left_sub_V, 1, true, temp, left_vert_index);
		igl::sort(right_sub_V, 1, true, temp, right_vert_index);

		// find the base LR-edge
		int l_vidx = 0, r_vidx = 0;
		while (l_vidx < left_vert_index.rows() && r_vidx < right_vert_index.rows()) {
			RowVectorXd p1 = left_sub_V.row(left_vert_index(l_vidx,0));
			RowVectorXd p2 = right_sub_V.row(right_vert_index(r_vidx,0));
			bool has_intersection = false;
			for (int i = 0; i < node->left->index.rows(); i++) {
				// need to skip self
				RowVectorXd p3 = V.row(get<0>(node->left->edges.at(i)));
				RowVectorXd p4 = V.row(get<1>(node->left->edges.at(i)));
				if (intersect(p1, p2, p3, p4)) {
					l_vidx ++;
					has_intersection = true;
					break;
				}
			}
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

		vector<tuple<int,int>> new_edges;
		new_edges.push_back(make_tuple(node->left->index(l_vidx), node->right->index(r_vidx)));

		int l_candidate_final = -1;
		int r_candidate_final = -1;
		int l_ridx = node->left->index(l_vidx);
		int r_ridx = node->right->index(r_vidx);

		// if no more edges are found, triangulation is complete
		while (l_candidate_final != -1 && r_candidate_final != -1) {
			l_candidate_final = -1;
			r_candidate_final = -1;
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
				double a = angle(V.row(r_ridx), V.row(other), V.row(l_ridx));

				// add it to a specific location based on angle
				for (int j = 0; j < r_candidate.size(); j++) {
					if (a < get<1>(r_candidate.at(j))){
						r_candidate.insert(r_candidate.begin() + j, make_tuple(other, a));
						break;
					}
					// add it to the end if not added before
					if (j == r_candidate.size()) {
						r_candidate.push_back(make_tuple(other, a));
					}
				}
			}

			for (int i = 0; i < node->left->edges.size(); i++) {
				int other;

				// find edges contains left node of base LR edge
				if (get<0>(node->left->edges.at(i)) == r_ridx){
					other = get<1>(node->left->edges.at(i));
				}
				else if (get<1>(node->left->edges.at(i)) == r_ridx) {
					other = get<0>(node->left->edges.at(i));
				}
				else {
					continue;
				}
				// to do: need to check angle clock vs counter-clock
				double a = angle(V.row(l_ridx), V.row(other), V.row(r_ridx));

				// add it to a specific location based on angle
				for (int j = 0; j < l_candidate.size(); j++) {
					if (a < get<1>(l_candidate.at(j))){
						l_candidate.insert(l_candidate.begin() + j, make_tuple(other, a));
						break;
					}
					// add it to the end if not added before
					if (j == l_candidate.size()) {
						l_candidate.push_back(make_tuple(other, a));
					}
				}
			}

			// find final candidate from both side
			RowVectorXd p_left = V.row(l_ridx);
			RowVectorXd p_right = V.row(r_ridx);
			RowVectorXd center;
			int l_candidate_final = -1;
			int r_candidate_final = -1;
			for (int i = 0; i < r_candidate.size() - 1; i++) {
				// terminate if angle is greater than 180
				if (get<1>(r_candidate.at(i)) > 3.1415926535) {
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
					remove_edge_from_node(node, r_ridx, get<0>(r_candidate.at(i)));
				}
			}

			for (int i = 0; i < l_candidate.size() - 1; i++) {
				if (get<1>(l_candidate.at(i)) > 3.1415926535) {
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(l_candidate.at(i))), center);
				if ((p_left-center).dot(p_left-center) < (V.row(get<0>(l_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(l_candidate.at(i+1))))) {
					l_candidate_final = get<0>(l_candidate.at(i));
					break;
				}
				else {
					remove_edge_from_node(node, l_ridx, get<0>(l_candidate.at(i)));
				}
			}

			// final step
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
			Polygon* new_polygon;
			new_polygon->size = 3;
			new_polygon->vertex = VectorXi::Zero(3);
			new_polygon->vertex(0) = r_ridx;
			new_polygon->adjacent_polygon.push_back(NULL);
			new_polygon->vertex(1) = l_ridx;
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
		}
		// add all new edges into the list
		node->edges.insert(node->edges.end(), node->left->edges.begin(), node->left->edges.end());
		node->edges.insert(node->edges.end(), node->right->edges.begin(), node->right->edges.end());
		node->edges.insert(node->edges.end(), new_edges.begin(), new_edges.end());
	}

}


// need to determine clockwise angle or counter-clockwise angle
double igl::bol::angle(RowVectorXd p1, RowVectorXd p2, RowVectorXd p3) {
	Eigen::RowVectorXd p1d = p1;
	Eigen::RowVectorXd p2d = p3;
	Eigen::RowVectorXd p3d = p1;
	return acos((p2d - p1d).dot(p3d - p1d) / ((p2d - p1d).norm() * (p3d - p1d).norm()));
}



bool igl::bol::intersect(RowVector3d p1, RowVector3d p2, RowVector3d p3, RowVector3d p4) {
	RowVector3d v1 = p2 - p1;
	RowVector3d v2 = p4 - p3;
	return (v1.cross(v2)).dot(p3 - p1) == 0;
}

void igl::bol::get_circle_center(RowVector3d p1, RowVector3d p2, RowVector3d p3, RowVector3d center) {
	RowVector3d a = (p1 + p2)/2;
	RowVector3d b = (p2 + p3)/2;
	RowVector3d up = (p1-p2).cross(p1-p3);
	RowVector3d c = up.cross(p1-p2);
	RowVector3d d = up.cross(p2-p3);
	double m = (dmnop(a,c,d,c)*dmnop(d,c,b,a)-dmnop(a,c,b,c)*dmnop(d,c,d,c)) / 
			   (dmnop(b,a,b,a)*dmnop(d,c,d,c)-dmnop(d,c,b,a)*dmnop(d,c,b,a));
	center = p1 + m*(p2-p1);
}

double igl::bol::dmnop(RowVector3d m, RowVector3d n, RowVector3d o, RowVector3d p){
	//replace with rational later
	return (m - n).dot(o - p);
}