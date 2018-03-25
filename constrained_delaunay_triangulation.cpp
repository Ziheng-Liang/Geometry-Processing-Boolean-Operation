#include "constrained_delaunay_triangulation.h"
#include <igl/sort.h>

void contruct_tree(MatrixXr V, Node* node) {
	using namespace std;
	using namespace eigen;
	assert (V.rows() >= 3);
	RowVectorXr p1 = V.rows(0);
	RowVectorXr p2 = V.rows(1);
	RowVectorXr p3 = V.rows(2);
	RowVectorXr up = (p2-p1).corss(p3-p1);
	RowVectorXr x = (p2-p2);
	RowVectorXr y = x.corss(up);
	MatrixXr Vn = V.rowwise() - p1;
	MatrixXr A = MatrixXr::Zero(2, 3);
	A << x(0), x(1), x(2), y(0), y(1), y(2);
	MatrixXr piA = A.transpose()*(A * A.transpose()).inverse();
	MatrixXr projectV = V * piA;
	MatrixXr temp;
	MatrixXi xindex;
	igl::sort(projectV, 1, temp, xindex);
	subdivide(xindex, node);
}

void constrained_delaunay_triangulation(MatrixXr V, Eigen::MatrixXi C) {
	Node root = new Node();
	contruct_tree(V, root);

}

void subdivide(RowVectorXi index, Node* node) {
	if (index.rows()) <= 3 {
		return;
	}
	node.index = index;
	int n = index.rows();
	Node left_node = new Node();
	node->left = left_node;
	subdivide(index.head(n/2), left_node);

	Node right_node = new Node();
	node->right = right_node;
	subdivide(index.tail(n-n/2), end, right_node);
}

void build_edges(MatrixXr V, Node* node, std::Vector<std::tuple<int,int>> edge_list) {
	if (node->left) {
		build_edges(V, node->left, edge_list);
	}
	if (node->right) {
		build_edges(V, node->right, edge_list);
	}

void get_circle_center(RowVectorXr p1, RowVectorXr p2, RowVectorXr p3, RowVectorXr center) {
	RowVectorXr pa = (p1 + p2)/2;
	RowVectorXr pb = (p2 + p3)/2;
	RowVectorXr up = (p1-p2).cross(p1-p3);
	RowVectorXr pc = up.corss(p1-p2);
	RowVectorXr pd = up.corss(p2-p3);
	double m = (dmnop(a,c,d,c)*dmnop(d,c,b,a)-dmnop(a,c,b,c)*dmnop(d,c,d,c)) / 
			   (dmnop(b,a,b,a)*dmnop(d,c,d,c)-dmnop(d,c,b,a)*dmnop(d,c,b,a));
	center = p1 + m*(p2-p1);
}

double dmnop(RowVectorXr m, RowVectorXr n, RowVectorXr o, RowVectorXr p){
	//replace with rational later
	return (m - n).dot(o - p);
}