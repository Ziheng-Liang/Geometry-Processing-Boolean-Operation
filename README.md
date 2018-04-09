# Geometry-Processing-Boolean-Operation
A toolkit for some basics operations used in mesh boolean operation

## Dependency
- [libigl](https://github.com/libigl/libigl)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [boost/multiprecision](https://www.boost.org/)

## Example
Triangle triangle intersection:
```
#include <Eigen/Dense>
#include <triangle.h>
#include <iostream>
#include <vector>
int main()
{
  using namespace igl::bol;
  Matrix33r A, B;
  A(0,0) = 0;
  A(0,1) = 0;
  A(0,2) = 0;

  A(1,0) = 0;
  A(1,1) = 2;
  A(1,2) = 0;

  A(2,0) = 2;
  A(2,1) = 0;
  A(2,2) = 0;

  B(0,0) = 0;
  B(0,1) = 1;
  B(0,2) = 1;

  B(1,0) = 0;
  B(1,1) = 1;
  B(1,2) = -1;

  B(2,0) = 1;
  B(2,1) = 1;
  B(2,2) = -1;
  std::vector<RowVector3r> r = igl::bol::t2t_intersect(A, B);
  for (int i = 0; i<r.size(); i++){
    std::cout << r[i] << std::endl;
  }

}
```
