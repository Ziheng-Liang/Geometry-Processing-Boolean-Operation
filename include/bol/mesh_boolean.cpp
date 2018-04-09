#include "mesh_boolean.h"
#include "remesh_self_intersections.h"



template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC>
bool mesh_boolean(
  const Eigen::MatrixBase<DerivedVA > & VA,
  const Eigen::MatrixBase<DerivedFA > & FA,
  const Eigen::MatrixBase<DerivedVB > & VB,
  const Eigen::MatrixBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC)
{
  using namespace igl::bol;
  //First combine VA FA, VB, FB into VV, FF
  Eigen::Matrix<typename DerivedVA::Scalar,Eigen::Dynamic,3> VV(VA.rows() + VB.rows(), 3);
  DerivedFC FF(FA.rows() + FB.rows(), 3);
  // Can't use comma initializer
  for(int a = 0;a<VA.rows();a++)
  {
    for(int d = 0;d<3;d++) VV(a,d) = VA(a,d);
  }
  for(int b = 0;b<VB.rows();b++)
  {
    for(int d = 0;d<3;d++) VV(VA.rows()+b,d) = VB(b,d);
  }
  FF.block(0, 0, FA.rows(), 3) = FA;
  // Eigen struggles to assign nothing to nothing and will assert if FB is empty
  if(FB.rows() > 0)
  {
    FF.block(FA.rows(), 0, FB.rows(), 3) = FB.array() + VA.rows();
  }
  
  //stage 1: remesh_self intersection on VV FF
  //call remesh self intersection
  //remesh_self_intersection
  
  //stage 2: extract cell
  
  
  //stage 3: winding number propagation


}
