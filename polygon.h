#ifndef IGL_BOL_Polygon_H
#define IGL_BOL_Polygon_H
namespace igl
{
    namespace bol
    {
        struct Polygon {
            int size;
            Eigen::MatrixXi edges;
            std::vector<Polygon*> adjacent_polygon;
        };

        // void (Polygon ab, Polygon a, Polygon b, std::tuple<int,int>> new_edge) {
        	
        // }

    }
}


#endif