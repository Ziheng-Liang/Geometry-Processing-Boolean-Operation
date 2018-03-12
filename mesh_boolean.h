
#ifndef IGL_BOL_MESH_BOOLEAN_H
#define IGL_BOL_MESH_BOOLEAN_H

#include <igl/MeshBooleanType.h>
namespace igl
{
	namespace bol
	{
      //  Inputs:
      //    VA  #VA by 3 list of vertex positions of first mesh
      //    FA  #FA by 3 list of triangle indices into VA
      //    VB  #VB by 3 list of vertex positions of second mesh
      //    FB  #FB by 3 list of triangle indices into VB
      //    type  type of boolean operation
      //  Outputs:
      //    VC  #VC by 3 list of vertex positions of boolean result mesh
      //    FC  #FC by 3 list of triangle indices into VC
      //  Returns true ff inputs induce a piecewise constant winding number
      //    field and type is valid
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
          Eigen::PlainObjectBase<DerivedFC > & FC);
	}
}
#endif
