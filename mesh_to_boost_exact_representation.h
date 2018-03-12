
#ifndef IGL_BOL_MESH_TO_BOOST_EXACT_REPRESENTATION_H
#define IGL_BOL_MESH_TO_BOOST_EXACT_REPRESENTATION_H


namespace igl
{
      // Convert a mesh (V,F) to a list of boost exact arithmetics
      //
      // Inputs:
      //   V  #V by 3 list of vertex positions
      //   F  #F by 3 list of triangle indices
      // Outputs:
      //   TO-DO: decide what should we put here
      template <
        typename DerivedV,
        typename DerivedF>
      void mesh_to_boost_exact_representation(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedF> & F);

}

#endif
