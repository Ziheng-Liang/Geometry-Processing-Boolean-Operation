#ifndef IGL_BOL_Polygon_H
#define IGL_BOL_Polygon_H

#include <Eigen/Dense>
#include <vector>


namespace igl
{
    namespace bol
    {
        struct Polygon {
            int size; // number of vertex
            Eigen::VectorXi vertex; // ordered vertex list
            std::vector<Polygon*> adjacent_polygon;
            // self = adjacent_polygon.at(i)->adjacent_polygon.at(adjacent_index.at(i))
            std::vector<int> adjacent_index; 
        };

        inline void split (Polygon* ab, Polygon* a, Polygon* b, std::tuple<int,int> new_edge) {
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
                    a->adjacent_index.push_back(i);
                    ab->adjacent_polygon.at(i)->adjacent_polygon.at(ab->adjacent_index.at(i)) = a;
                }
                if (i == idx1) {
                    a->adjacent_polygon.push_back(b);
                    a->adjacent_index.push_back(i);
                    ab->adjacent_polygon.at(i)->adjacent_polygon.at(ab->adjacent_index.at(i)) = a;
                }
                if (i >= idx1 && i < idx2) {
                    b->adjacent_polygon.push_back(ab->adjacent_polygon.at(i));
                    b->adjacent_index.push_back(i);
                    ab->adjacent_polygon.at(i)->adjacent_polygon.at(ab->adjacent_index.at(i)) = b;
                }
                if (i == idx2) {
                    b->adjacent_polygon.push_back(a);
                    b->adjacent_index.push_back(i);
                    ab->adjacent_polygon.at(i)->adjacent_polygon.at(ab->adjacent_index.at(i)) = b;
                }
                if (i >= idx2) {
                    a->adjacent_polygon.push_back(ab->adjacent_polygon.at(i));
                    a->adjacent_index.push_back(i);
                    ab->adjacent_polygon.at(i)->adjacent_polygon.at(ab->adjacent_index.at(i)) = a;
                }
            }
        }

        inline void merge (Polygon* a, Polygon* b, Polygon* ab) {
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
                    ab->adjacent_index.push_back(a->adjacent_index.at(i));
                }
                else if (i <= idx1 + b->size - 1) {
                    ab->vertex(i) = b->vertex((i - idx1 + idx2 + 1)%b->size);
                    ab->adjacent_polygon.push_back(a->adjacent_polygon.at((i - idx1 + idx2 + 1)%b->size));
                    ab->adjacent_index.push_back(a->adjacent_index.at((i - idx1 + idx2 + 1)%b->size));
                }
                else {
                    ab->vertex(i) = a->vertex(i + 1 - b->size);
                    ab->adjacent_polygon.push_back(a->adjacent_polygon.at(i + 1 - b->size));
                    ab->adjacent_index.push_back(a->adjacent_index.at(i));
                }
                ab->adjacent_polygon.at(i)->adjacent_index.at(ab->adjacent_index.at(i))= i;
            }
        }


        inline int find_vertex(Polygon* p, int index) {
            for (int i = 0; i < p->size; i++) {
                if (p->vertex(i) == index) {
                    return i;
                }
            }
            return -1;
        }
        inline int exist_edges(Polygon* p, int i, int j) {
            for (int k = 0; k < p->size - 1; k++) {
                if ((p->vertex(k) == i && p->vertex(k+1) == j) ||
                    (p->vertex(k) == j && p->vertex(k+1) == i)) {
                    return i;
                }
            }
            return -1;
        }
    }
}


#endif