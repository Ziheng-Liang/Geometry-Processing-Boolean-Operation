#include "constrained_delaunay_triangulation.h"
#include <igl/sort.h>
#include <igl/slice.h>
#include <igl/colon.h>
#include <math.h> 
#include <iostream>
#include <igl/unique_simplices.h>
using namespace igl::bol;
using namespace Eigen;
using namespace std;

void igl::bol::constrained_delaunay_triangulation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &C, Eigen::MatrixXi & F) {
	Node* root = new Node();
	MatrixXd projectV;
	contruct_tree(V, root, projectV);
	delaunay_triangulation(projectV, root);
	add_constrained(projectV, C, root);
	MatrixXi FF = Eigen::MatrixXi::Zero(root->polygons.size(), 3);
	for (int i = 0; i < FF.rows(); i++) {
		FF(i, 0) = root->polygons.at(i)->vertex(0);
		FF(i, 1) = root->polygons.at(i)->vertex(1);
		FF(i, 2) = root->polygons.at(i)->vertex(2);
	}
	igl::unique_simplices(FF, F);
}

void igl::bol::contruct_tree(const MatrixXd &V, Node* node, MatrixXd &projectV) {
	assert (V.rows() >= 3);
	RowVector3d p1 = V.row(0);
	RowVector3d p2 = V.row(1);
	RowVector3d p3 = V.row(2);
	RowVector3d up = (p2-p1).cross(p3-p1);
	RowVector3d x = (p2-p1).normalized();
	RowVector3d y = up.cross(x).normalized();
	MatrixXd Vn = V.rowwise() - p1;
	MatrixXd A = MatrixXd::Zero(2, 3);
	A << x(0), x(1), x(2), y(0), y(1), y(2);
	MatrixXd piA = A.transpose()*(A * A.transpose()).inverse();
	projectV = V * piA;
	MatrixXd temp;
	MatrixXi xindex;
	projectV = MatrixXd::Zero(V.rows(), 2);
	projectV.col(0) = V.col(0);
	projectV.col(1) = V.col(1);
	igl::sort(projectV, 1, true, temp, xindex);
	subdivide(xindex.col(0), node);
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
	Polygon *start;
	Polygon *next;
	bool intersection = false;

	// loop through all polygons
	for (int i = 0; i < polygons.size(); i++) {
		// find polygon that contains one point of the constrain
		vidx1 = find_vertex(polygons.at(i), C(0));
		if (vidx1 != -1) {
			// check if any edge intersect the constrain edge
			for (int j = 1; j < polygons.at(i)->size; j++) {
				if (intersect(V.row(C(0)), V.row(C(1)), 
							  V.row(polygons.at(i)->vertex(j-1)), 
							  V.row(polygons.at(i)->vertex(j)))) {
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
			for (int j = 1; j < polygons.at(i)->size; j++) {
				if (intersect(V.row(C(0)), V.row(C(1)), 
							  V.row(polygons.at(i)->vertex(j-1)), 
							  V.row(polygons.at(i)->vertex(j)))) {
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
	if (!start && !next) {
		return;
	}
	polygons.erase(std::remove(polygons.begin(), polygons.end(), start), polygons.end());
	Polygon* new_polygon = new Polygon();
	do {
		merge(start, next, new_polygon);
		for (int j = 0; j < new_polygon->size; j ++) {
			if (new_polygon->adjacent_polygon.at(j)){
			}
			else {
			}
		}

		polygons.erase(std::remove(polygons.begin(), polygons.end(), next), polygons.end());
		for (int j = 0; j < new_polygon->size; j++) {
			if (intersect(V.row(C(0)), V.row(C(1)), 
						  V.row(new_polygon->vertex(j)), 
						  V.row(new_polygon->vertex((j+1)%new_polygon->size)))) {
				next = new_polygon->adjacent_polygon.at(j);
				start = new_polygon;
				break;
			}
		}

		// if both edge are inside vertex, done
		vidx1 = find_vertex(new_polygon, C(0));
		vidx2 = find_vertex(new_polygon, C(1));
		if (vidx1==-1 || vidx2==-1) {
			new_polygon = new Polygon();
		}
	}while ((vidx1 == -1 || vidx2 == -1));
	Polygon *a = new Polygon();
	Polygon *b = new Polygon();
	split(new_polygon, a, b, C(0), C(1));

	MatrixXi xindex;
	MatrixXi temp;
	igl::sort(a->vertex, 1, true, temp, xindex);
	// igl::slice(V,node->left->index,cols,left_sub_V);
	// igl::slice(V,node->right->index,cols,right_sub_V);

	Node* node = new Node();
	subdivide(temp, node);
	delaunay_triangulation(V, node);
	polygons.insert(polygons.end(), node->polygons.begin(), node->polygons.end());

	igl::sort(b->vertex, 1, true, temp, xindex);
	// igl::slice(V,node->left->index,cols,left_sub_V);
	// igl::slice(V,node->right->index,cols,right_sub_V);

	node = new Node();
	subdivide(temp, node);
	delaunay_triangulation(V, node);
	polygons.insert(polygons.end(), node->polygons.begin(), node->polygons.end());
	MatrixXi F = Eigen::MatrixXi::Zero(polygons.size(), 3);
	for (int i = 0; i < F.rows(); i++) {
		F(i, 0) = polygons.at(i)->vertex(0);
		F(i, 1) = polygons.at(i)->vertex(1);
		F(i, 2) = polygons.at(i)->vertex(2);
	}
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
			if (colinear(V.row(node->index(0)), V.row(node->index(1)), V.row(node->index(2)))){
				for (int i = 0; i < node->size; i++) {
					if (((V.row(node->index(i)) - V.row(node->index((i+1)%3)))
						.dot(V.row(node->index(i)) - V.row(node->index((i+2)%3)))) < 0) {
						node->edges.push_back(make_tuple(node->index((i+2) % node->size), node->index(i)));
						node->edges.push_back(make_tuple(node->index(i), node->index((i+1) % node->size)));
						node->size = 2;
					}
				}
			}
			else {
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
			RowVectorXd p1 = left_sub_V.row(left_vert_index(l_vidx,1));
			RowVectorXd p2 = right_sub_V.row(right_vert_index(r_vidx,1));
			bool has_intersection = false;
	
			for (int i = 0; i < node->left->edges.size(); i++) {
				// need to skip self
				RowVectorXd p3 = V.row(get<0>(node->left->edges.at(i)));
				RowVectorXd p4 = V.row(get<1>(node->left->edges.at(i)));
		
				if (intersect(p1, p2, p3, p4)) {
					l_vidx ++;
					has_intersection = true;
					break;
				}
			}
	
			for (int i = 0; i < node->right->edges.size(); i++) {
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
		new_edges.push_back(make_tuple(node->left->index(left_vert_index(l_vidx,1)), node->right->index(right_vert_index(r_vidx,1))));

		int l_ridx = node->left->index(left_vert_index(l_vidx,1));
		int r_ridx = node->right->index(right_vert_index(r_vidx,1));



		// if no more edges are found, triangulation is complete
		int l_candidate_final;
		int r_candidate_final;
		do {
			// find candidates from both side
			vector<tuple<int,double>> r_candidate;
			vector<tuple<int,double>> l_candidate;
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
	
	
	
	
	
				if (a > 3.1415926535) {
					continue;
				}
				if (a < 0) {
					continue;
				}
		
		
				// add it to a specific location based on angle
				bool added = false;
				for (int j = 0; j < r_candidate.size(); j++) {
					if (other == get<0>(r_candidate.at(j))) {
						added = true;
						break;
					}
					else if (a < get<1>(r_candidate.at(j))){
						r_candidate.insert(r_candidate.begin() + j, make_tuple(other, a));
						added = true;
						break;
					}
					// add it to the end if not added before
				}
				if (!added) {
					r_candidate.push_back(make_tuple(other, a));
				}
			}

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
		
				if (a > 3.1415926535) {
					continue;
				}
				if (a < 0) {
					continue;
				}
				bool added = false;
				// add it to a specific location based on angle
				for (int j = 0; j < l_candidate.size(); j++) {
					if (other == get<0>(l_candidate.at(j))) {
						added = true;
						break;
					}
					else if (a < get<1>(l_candidate.at(j))){
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

			l_candidate_final = -1;
			r_candidate_final = -1;
			for (int i = 0; i < r_candidate.size(); i++) {
				if ((i + 1) == r_candidate.size()) {
					r_candidate_final = get<0>(r_candidate.at(i));
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(r_candidate.at(i))), center);
				// if next candidate is outside, this is the final candidate
				if ((p_left-center).dot(p_left-center) <= (V.row(get<0>(r_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(r_candidate.at(i+1))) - center)) {
					r_candidate_final = get<0>(r_candidate.at(i));
					break;
				}
				// if not, remove this edge
				else {
					remove_edge_from_node(node->right, r_ridx, get<0>(r_candidate.at(i)));
				}
			}

			if (l_candidate.size() == 1) {
				l_candidate_final = get<0>(l_candidate.at(0));
			}
			for (int i = 0; i < l_candidate.size(); i++) {
		
				if ((i + 1) == l_candidate.size()) {
					l_candidate_final = get<0>(l_candidate.at(i));
					break;
				}
				get_circle_center(p_left, p_right, V.row(get<0>(l_candidate.at(i))), center);
				if ((p_left-center).dot(p_left-center) <= (V.row(get<0>(l_candidate.at(i+1))) - center)
					.dot(V.row(get<0>(l_candidate.at(i+1))) - center)) {
					l_candidate_final = get<0>(l_candidate.at(i));
					break;
				}
				else {
			
					remove_edge_from_node(node->left, l_ridx, get<0>(l_candidate.at(i)));
			
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

			Polygon* new_polygon = new Polygon();
			new_polygon->size = 3;
			new_polygon->vertex = VectorXi::Zero(3);
			new_polygon->vertex(0) = r_ridx;
			new_polygon->vertex(1) = l_ridx;
			// up is always NULL, will be updated in the future

			new_polygon->adjacent_polygon.push_back(NULL);
			new_polygon->adjacent_index.push_back(NULL);
			if (l_ridx == get<0>(new_edges.at(new_edges.size()-1))) {
				new_polygon->vertex(2) = get<1>(new_edges.at(new_edges.size()-1));
				// if it is not bottom edege
				if (node->polygons.size() != 0) {
					new_polygon->adjacent_polygon.push_back(node->polygons.at(node->polygons.size()-1));
					new_polygon->adjacent_index.push_back(exist_edges(new_polygon->adjacent_polygon.at(1), 
																  new_polygon->vertex(1), 
																  new_polygon->vertex(2)));
				}
				// if it is bottom edge
				else {
					new_polygon->adjacent_polygon.push_back(NULL);
					new_polygon->adjacent_index.push_back(NULL);
				}
				bool added = false;
				for (int i = 0; i < node->right->polygons.size();i++) {
					int oppo_idx = exist_edges(node->right->polygons.at(i), 
											   new_polygon->vertex(2), 
											   new_polygon->vertex(0));
					if (oppo_idx != -1) {
						new_polygon->adjacent_polygon.push_back(node->right->polygons.at(i));
						new_polygon->adjacent_index.push_back(oppo_idx);
						added = true;
						break;
					}
				}
				if (!added) {
					new_polygon->adjacent_polygon.push_back(NULL);
					new_polygon->adjacent_index.push_back(NULL);
				}
			}
			else {
				new_polygon->vertex(2) = get<0>(new_edges.at(new_edges.size()-1));
	
				bool added = false;
				for (int i = 0; i < node->left->polygons.size();i++) {
					int oppo_idx = exist_edges(node->left->polygons.at(i), 
											   new_polygon->vertex(1), 
											   new_polygon->vertex(2));
					if (oppo_idx != -1) {
						new_polygon->adjacent_polygon.push_back(node->left->polygons.at(i));
						new_polygon->adjacent_index.push_back(oppo_idx);
						added = true;
						break;
					}
				}
				if (!added) {
					new_polygon->adjacent_polygon.push_back(NULL);
					new_polygon->adjacent_index.push_back(NULL);
				}
	
				if (node->polygons.size() != 0) {
					new_polygon->adjacent_polygon.push_back(node->polygons.at(node->polygons.size()-1));
					new_polygon->adjacent_index.push_back(exist_edges(new_polygon->adjacent_polygon.at(2), 
																			  new_polygon->vertex(2), 
																			  new_polygon->vertex(0)));
				}
				else {
					new_polygon->adjacent_polygon.push_back(NULL);
					new_polygon->adjacent_index.push_back(NULL);
				}
			}


			if (new_polygon->adjacent_polygon.size() >= 2 && new_polygon->adjacent_polygon.at(1)) {
				new_polygon->adjacent_polygon.at(1)->adjacent_index.at(new_polygon->adjacent_index.at(1)) = 1;
				new_polygon->adjacent_polygon.at(1)->adjacent_polygon.at(new_polygon->adjacent_index.at(1)) = new_polygon;
			}
			if (new_polygon->adjacent_polygon.size() >= 3 && new_polygon->adjacent_polygon.at(2)) {
				new_polygon->adjacent_polygon.at(2)->adjacent_index.at(new_polygon->adjacent_index.at(2)) = 2;
				new_polygon->adjacent_polygon.at(2)->adjacent_polygon.at(new_polygon->adjacent_index.at(2)) = new_polygon;
			}

			new_edges.push_back(make_tuple(l_ridx, r_ridx));
			node->polygons.push_back(new_polygon);
	
	
		} while ((l_candidate_final != -1) || (r_candidate_final != -1));
		// add all new edges into the list
		node->edges.insert(node->edges.end(), node->left->edges.begin(), node->left->edges.end());
		node->edges.insert(node->edges.end(), node->right->edges.begin(), node->right->edges.end());
		node->edges.insert(node->edges.end(), new_edges.begin(), new_edges.end());
		
		
		
		node->polygons.insert(node->polygons.end(), node->left->polygons.begin(), node->left->polygons.end());
		node->polygons.insert(node->polygons.end(), node->right->polygons.begin(), node->right->polygons.end());
		
		
		
		MatrixXi F = MatrixXi::Zero(node->polygons.size(), 3);
		for (int i = 0; i < F.rows(); i++) {
			F(i, 0) = node->polygons.at(i)->vertex(0);
			F(i, 1) = node->polygons.at(i)->vertex(1);
			F(i, 2) = node->polygons.at(i)->vertex(2);
		}
	}

}


// need to determine clockwise angle or counter-clockwise angle
// double igl::bol::angle(RowVectorXd p1, RowVectorXd p2, RowVectorXd p3) {
// 	/
// 	/
// 	/
// 	/
// 	/
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
	// RowVectorXd v1 = p2 - p1;
	// RowVectorXd v2 = p4 - p3;
	// RowVectorXd d = p3 - p1;
	// // colinear
	// double det = v1(0) * v2(1) - v1(1) * v2(0);
	// if (abs(det) < 1e-10) {
	// 	return false;
	// }

	// double r = (d(0) * v2(1) - d(1) * v2(0)) / det;
 //    double s = (v1(0) * d(1) - v1(0) * d(0)) / det;
	// return !(r < 0 || r > 1 || s < 0 || s > 1);
	return (((p1(0)-p3(0))*(p4(1)-p3(1)) - (p1(1)-p3(1))*(p4(0)-p3(0)))
            * ((p2(0)-p3(0))*(p4(1)-p3(1)) - (p2(1)-p3(1))*(p4(0)-p3(0))) < 0)
            &&
           (((p3(0)-p1(0))*(p2(1)-p1(1)) - (p3(1)-p1(1))*(p2(0)-p1(0)))
            * ((p4(0)-p1(0))*(p2(1)-p1(1)) - (p4(1)-p1(1))*(p2(0)-p1(0))) < 0);
}

// void igl::bol::get_circle_center(const RowVectorXd &p1, const RowVectorXd &p2, const RowVectorXd &p3, RowVectorXd &center) {
// 	assert(p1.cols() == 2);
// 	double ma = (p2(1) - p1(1)) / (p2(0) - p1(0));
// 	double mb = (p3(1) - p2(1)) / (p3(0) - p2(0));
// 	double x = (ma*mb*(p1(1)-p3(1)) + mb*(p1(0)+p2(0)) - ma*(p2(0)+p3(0)))/(2*(mb-ma));
// 	double y = -1/ma * (x - (p1(0) + p2(0))/2) + (p1(1)+p2(1))/2;
// 	center = VectorXd::Zero(2);
// 	center << x, y;
// }

bool igl::bol::colinear(const RowVectorXd &p1, const RowVectorXd &p2, const RowVectorXd &p3)  {
	return abs(p1(0) * (p2(1) - p3(1)) + p2(0) * (p3(1) - p1(1)) + p3(0) * (p1(1) - p2(1))) < 1e-10;
}

void igl::bol::get_circle_center(const RowVectorXd &p1, const RowVectorXd &p2, const RowVectorXd &p3, RowVectorXd &center) {
	RowVectorXd a = (p1 + p2)/2;
	RowVectorXd b = (p2 + p3)/2;
	RowVectorXd c(2);
	c << (-(p2 - p1))(1), (p2 - p1)(0);
	RowVectorXd d(2);
	d << (-(p3 - p2))(1), (p3 - p2)(0);
	double t = c(0) * d(1) - c(1) * d(0);
	double n = (b(0) - a(0)) * d(1)
                 - (b(1) - a(1)) * d(0);
    center = a + c*(n/t);
}