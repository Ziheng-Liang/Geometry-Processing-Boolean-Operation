#include "mesh_boolean.h"




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
  //First combine VA FA, VB, FB into VV, FF
  
  
  //stage 1: remesh_self intersection
  
  
  //stage 2: extract cell
  
  
  //stage 3: winding number propagation


}
