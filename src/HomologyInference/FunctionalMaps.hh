/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

template<typename T>
MatX<T>
map_using_fmap(
        const MatXd& _basis_A,
        const MatXd& _basis_B,
        const MatXd& _C,
        const MatX<T>& _data_A);

std::vector<ExternalProperty<VH, Complex>>
map_fields_fmap(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const MatXd& _basis_A,
        const MatXd& _basis_B,
        const MatXd& _C,
        const std::vector<ExternalProperty<VH, Complex>>& _fields_A);

}
