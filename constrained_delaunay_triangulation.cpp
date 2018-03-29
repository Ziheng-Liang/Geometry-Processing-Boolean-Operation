#include "constrained_delaunay_triangulation.h"
#include <igl/sort.h>
#include <igl/slice.h>
#include <igl/colon.h>
#include <math.h> 

using namespace igl::bol;
using namespace Eigen;
using namespace std;



void igl::bol::contruct_tree(MatrixXd V, Node* node) {
	using namespace std;
	using namespace Eigen;
	assert (V.rows() >= 3);
	RowVector3d p1 = V.row(0);
	RowVector3d p2 = V.row(1);
	RowVector3d p3 = V.row(2);
	RowVector3d up = (p2-p1).cross(p3-p1);
	RowVector3d x = (p2-p2);
	RowVector3d y = x.cross(up);
	MatrixXd Vn = V.rowwise() - p1;
	MatrixXd A = MatrixXd::Zero(2, 3);
	A << x(0), x(1), x(2), y(0), y(1), y(2);
	MatrixXd piA = A.transpose()*(A * A.transpose()).inverse();
	MatrixXd projectV = V * piA;
	MatrixXd temp;
	MatrixXi xindex;
	igl::sort(projectV, 1, true, temp, xindex);
	subdivide(xindex, node);
}

void igl::bol::constrained_delaunay_triangulation(Eigen::MatrixXd V, Eigen::MatrixXi C, Eigen::MatrixXi F) {
	Node* root = new Node();
	contruct_tree(V, root);
	delaunay_triangulation(V, root);
	add_constrained(V, C, root);
	F = Eigen::MatrixXi::Zero(root->polygons.size, 3);
	for (int i = 0; i < F.rows(); i++) {
		F(i, 0) = root->polygons.at(i)->vertex.at(0);
		F(i, 1) = root->polygons.at(i)->vertex.at(1);
		F(i, 2) = root->polygons.at(i)->vertex.at(2);
	}
}

void igl::bol::add_constrained(Eigen::MatrixXd V, Eigen::MatrixXi C, Node* node) {
	for (int i = 0; i < C.rows(); i++) {
		break_polygons(V, C.row(i), node->polygons)
	}
}

void igl::bol::break_polygons(Eigen::MatrixXd V, Eigen::RowVectorXi C, std::vector<Polygon*> polygons, Eigen::RowVectorXi index) {
	int vidx1 = -1, vidx2 = -1;
	polygon* start, next;
	bool intersection = false;
	for (int i = 0; i < polygons.size(); i++) {
		vidx1 = find(polygons.at(i), C(0));
		if (vidx1 != -1) {
			for (int j = 0; j < polygons.at(i)->size - 1; j++) {
				if (intersect(V.row(C(0)), V.row(C(1)), 
							  V.row(polygons.at(i)->vertex(j), 
							  V.row(polygons.at(i)->vertex(j+1))))) {
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
		vidx2 = find(polygons.at(i), C(1));
		if (vidx2 != -1) {
			for (int j = 0; j < polygons.at(i)->size - 1; j++) {
				if (intersect(V.row(C(0)), V.row(C(1)), 
							  V.row(polygons.at(i)->vertex(j), 
							  V.row(polygons.at(i)->vertex(j+1))))) {
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
	polygons.erase(std::remove(polygons.begin(), polygons.end(), start), vec.end());
	while (vidx1 == -1 || vidx2 == -1) {
		Polygon* new_polygon;
		merge(start, next, new_polygon);
		polygons.erase(std::remove(polygons.begin(), polygons.end(), next), vec.end());
		for (int j = 0; j < new_polygon->size - 1; j++) {
			if (intersect(V.row(C(0)), V.row(C(1)), 
						  V.row(new_polygon->vertex(j), 
						  V.row(new_polygon->vertex(j+1))))) {
				next = new_polygon->adjacent_polygon.at(j);
				start = new_polygon
				break;
			}
		}
		if (vidx1 == -1) {
			vidx1 = find(polygons.at(i), C(0));
		}
		if (vidx2 == -1) {
			vidx2 = find(polygons.at(i), C(1));
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


void igl::bol::delaunay_triangulation(MatrixXd V, Node* node) {
	if (node->left) {
		assert(node->right);
		build_edges(V, node->left);
	}
	if (node->right) {
		assert(node->left);
		build_edges(V, node->right);
	}
	if (!node->left && !node->right) {
		Polygon* polygon = new Polygon();
		polygon->size = node->size;
		VectorXi edges = VectorXi::Zero(node->size, 2);
		for (int i = 0; i < node->size; i++) {
			edges(i, 0) = node->index(i);
			edges(i, 1) = node->index((i+1) % node->size);
		}
		polygon->edges = edges;
		node->edges = edges;
		node->polygons.push_back(polygon);
	}
	else {
		MatrixXd left_sub_V;
		MatrixXd right_sub_V;
		VectorXi cols;
		igl::colon(0,V.cols()-1, cols);
		igl::slice(V,node->left->index,cols,left_sub_V);
		igl::slice(V,node->right->index,cols,right_sub_V);
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
				RowVectorXd p3 = V.row(node->left->edges(i, 0));
				RowVectorXd p4 = V.row(node->left->edges(i, 1));
				if (intersect(p1, p2, p3, p4)) {
					l_vidx ++;
					has_intersection = true;
					break;
				}
			}
			for (int i = 0; i < node->right->index.rows(); i++) {
				RowVectorXd p3 = V.row(node->right->edges(i, 0));
				RowVectorXd p4 = V.row(node->right->edges(i, 1));
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
		MatrixXi new_edges;
		new_edges.resize(new_edges.rows()+1, 2);
		new_edges.row(new_edges.rows() - 1) << node->left->index(l_vidx), node->right->index(r_vidx);

		int l_candidate_final = -1;
		int r_candidate_final = -1;
		int l_ridx = node->left->index(l_vidx);
		int r_ridx = node->right->index(r_vidx);

		while (l_candidate_final != -1 || r_candidate_final != -1) {
			l_candidate_final = -1;
			r_candidate_final = -1;
			// find candidates from both side
			vector<tuple<int,int>> r_candidate;
			vector<tuple<int,int>> l_candidate;
			for (int i = 0; i < node->right->edges.rows(); i++) {
				int other;
				if (node->right->edges(i, 0) == r_ridx){
					other = node->right->edges(i, 1);
				}
				else if (node->right->edges(i, 1) == r_ridx) {
					other = node->right->edges(i, 0);
				}
				double a = angle(V.row(r_ridx), V.row(other), V.row(l_ridx));
				bool added = false;
				for (int j = 0; j < r_candidate.size(); j++) {
					if (a < get<1>(r_candidate.at(j))){
						r_candidate.insert(r_candidate.begin() + j, make_tuple(other, a));
						added = true;
						break;
					}
				}
				if (!added) {
					r_candidate.push_back(make_tuple(other, a));
				}
			}

			for (int i = 0; i < node->left->edges.rows(); i++) {
				int other;
				if (node->left->edges(i, 0) == l_ridx){
					other = node->left->edges(i, 1);
				}
				else if (node->left->edges(i, 1) == l_ridx) {
					other = node->left->edges(i, 0);
				}
				double a = angle(V.row(l_ridx), V.row(other), V.row(r_ridx));
				bool added = false;
				for (int j = 0; j < l_candidate.size(); j++) {
					if (a < get<1>(l_candidate.at(j))){
						l_candidate.insert(l_candidate.begin() + j, make_tuple(other, a));
						added = true;
						break;
					}
				}
				if (!added) {
					l_candidate.push_back(make_tuple(other, a));
				}
			}

			// find final candidate from both side
			RowVectorXd p_left = V.row(l_ridx);
			RowVectorXd p_right = V.row(r_ridx);
			RowVectorXd center;
			int l_candidate_final = -1;
			int r_candidate_final = -1;
			for (int i = 0; i < r_candidate.size() - 1; i++) {
				if (get<1>(r_candidate.at(i)) > 3.1415926535) {
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(r_candidate.at(i))), center);
				if ((p_left-center).dot(p_left-center) < (V.row(get<0>(r_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(r_candidate.at(i+1))))) {
					r_candidate_final = get<0>(r_candidate.at(i));
					break;
				}
				else {
					remove_edge_from_node(node->polygons, p_right, get<0>(r_candidate.at(i)));
				}
			}

			for (int i = 0; i < l_candidate.size() - 1; i++) {
				if (get<1>(l_candidate.at(i)) > 3.1415926535) {
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(l_candidate.at(i))), center);
				if ((p_left-center).dot(p_left-center) < (V.row(get<0>(r_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(r_candidate.at(i+1))))) {
					r_candidate_final = get<0>(l_candidate.at(i));
					break;
				}
				else {
					remove_edge_from_node(node->polygons, p_left, get<0>(l_candidate.at(i)));
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
			new_edges.resize(new_edges.rows()+1, 2);
			new_edges.row(new_edges.rows() - 1) << l_ridx, r_ridx;
		}
		
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