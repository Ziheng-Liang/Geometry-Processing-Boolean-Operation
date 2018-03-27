#ifndef IGL_BOL_Polygon_H
#define IGL_BOL_Polygon_H
namespace igl
{
    namespace bol
    {
        struct Polygon {
            int size; // number of vertex
            Eigen::VectorXi vertex; // ordered vertex list
            std::vector<Polygon*> adjacent_polygon;
        };

        void split (Polygon* ab, Polygon* a, Polygon* b, std::tuple<int,int> new_edge) {
            int idx1=-1, idx2=-1;
        	for (int i = 0; i < ab->size; i++) {
                int current = ab->vertex(i);
                if (current == std::get<0>(new_edge)) {
                    idx1 = i;
                }
                else if (current == std::get<1>(new_edge)) {
                    idx2 = i;
                }
            }
            if (idx2 < idx1) {
                std::swap(idx1,idx2);
            }
            assert(idx1!=-1);
            assert(idx2!=-1);
            assert(idx1!=-idx2);
            b->size = idx2 - idx1 + 1;
            a->size = ab->size - b->size + 2;
            a->vertex = Eigen::VectorXi::Zero(a->size);
            b->vertex = Eigen::VectorXi::Zero(b->size); 
            a->vertex.head(idx1+1) = ab->vertex.head(idx1+1);
            a->vertex.tail(a->size-idx1-1) = ab->vertex.tail(a->size-idx1-1);
            b->vertex = ab->vertex.segment(idx1, idx2+1);

            for (int i = 0; i < ab->size; i++) {
                if (i < idx1) {
                    a->adjacent_polygon.push_back(ab->adjacent_polygon.at(i));
                }
                if (i == idx1) {
                    a->adjacent_polygon.push_back(b);
                }
                if (i >= idx1 && i < idx2) {
                    b->adjacent_polygon.push_back(ab->adjacent_polygon.at(i));
                }
                if (i == idx2) {
                    b->adjacent_polygon.push_back(a);
                }
                if (i >= idx2) {
                    a->adjacent_polygon.push_back(ab->adjacent_polygon.at(i));
                }
            }
        }

        void merge (Polygon* a, Polygon* b, Polygon* ab) {
            int idx1, idx2;
            for (int i = 0; i < a->size; i++) {
                if (a->adjacent_polygon.at(i) == b) {
                    idx1 = i;
                    break;
                }
            }
            for (int i = 0; i < b->size; i++) {
                if (b->adjacent_polygon.at(i) == a) {
                    idx2 = i;
                    break;
                }
            }
            assert(idx1!=-1);
            assert(idx2!=-1);
            ab->size = a->size + b->size - 2;
            ab->vertex = Eigen::VectorXi::Zero(ab->size);
            for (int i = 0; i < ab->size; i++) {
                if (i < idx1) {
                    ab->vertex(i) = a->vertex(i); 
                    ab->adjacent_polygon.push_back(a->adjacent_polygon.at(i));
                }
                else if (i <= idx1 + b->size - 1) {
                    ab->vertex(i) = b->vertex((i - idx1 + idx2 + 1)%b->size);
                    ab->adjacent_polygon.push_back(a->adjacent_polygon.at((i - idx1 + idx2 + 1)%b->size));
                }
                else {
                    ab->vertex(i) = a->vertex(i + 1 - b->size);
                    ab->adjacent_polygon.push_back(a->adjacent_polygon.at(i + 1 - b->size));
                }
            }
        } 
    }
}


#endif