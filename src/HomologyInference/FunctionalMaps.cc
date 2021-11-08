/*
 * Author: Janis Born
 */
#include "FunctionalMaps.hh"

#include <HomologyInference/Utils/Out.hh>
#include <HomologyInference/Utils/Timer.hh>

namespace HomologyInference
{

template<typename T>
MatX<T>
map_using_fmap(
        const MatXd& _basis_A,
        const MatXd& _basis_B,
        const MatXd& _C,
        const MatX<T>& _data_A)
{
    ISM_ASSERT_EQ(_data_A.rows(), _basis_A.rows());
    ISM_ASSERT_EQ(_basis_A.cols(), _basis_B.cols());
    const int k = _basis_A.cols();
    ISM_ASSERT_EQ(_C.rows(), k);
    ISM_ASSERT_EQ(_C.cols(), k);

    // Use pseudo-inverse of basis_A to convert _data_A to coefficient vector a.
    auto basis_A_inv = (_basis_A.transpose() * _basis_A).inverse() * _basis_A.transpose();
    MatX<T> a = basis_A_inv * _data_A;

    MatX<T> result = _basis_B * _C * a;
    ISM_ASSERT_EQ(result.rows(), _basis_B.rows());
    ISM_ASSERT_EQ(result.cols(), _data_A.cols());
    return result;
}
template MatX<double> map_using_fmap<double>(const MatXd&, const MatXd&, const MatXd&, const MatX<double>&);
template MatX<Complex> map_using_fmap<Complex>(const MatXd&, const MatXd&, const MatXd&, const MatX<Complex>&);

std::vector<ExternalProperty<VH, Complex>>
map_fields_fmap(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const MatXd& _basis_A,
        const MatXd& _basis_B,
        const MatXd& _C,
        const std::vector<ExternalProperty<VH, Complex>>& _fields_A)
{
    Timer timer(__FUNCTION__);

    ISM_ASSERT_EQ(_basis_A.rows(), _mesh_A.n_vertices());
    ISM_ASSERT_EQ(_basis_B.rows(), _mesh_B.n_vertices());
    ISM_ASSERT_EQ(_basis_A.cols(), _basis_B.cols());
    const int k = _basis_A.cols();
    ISM_ASSERT_EQ(_C.rows(), k);
    ISM_ASSERT_EQ(_C.cols(), k);

    // Convert fields properties to matrix. Each column is one field.
    MatX<Complex> f_A(_mesh_A.n_vertices(), _fields_A.size());
    for (int row = 0; row < f_A.rows(); ++row)
        for (int col = 0; col < f_A.cols(); ++col)
            f_A(row, col) = _fields_A[col][VH(row)];

    // Apply functional map.
    MatX<Complex> f_B = map_using_fmap(_basis_A, _basis_B, _C, f_A);
    ISM_ASSERT_EQ(f_B.rows(), _mesh_B.n_vertices());
    ISM_ASSERT_EQ(f_B.cols(), f_A.cols());

    // Convert back from matrix to properties.
    std::vector<ExternalProperty<VH, Complex>> fields_B;
    for (int col = 0; col < f_B.cols(); ++col)
    {
        ExternalProperty<VH, Complex> field_B(_mesh_B);
        for (const auto& vh_B : _mesh_B.vertices())
        {
            const int row = vh_B.idx();
            const Complex& v = f_B(row, col);
            const Complex v_normalized = v / std::abs(v);
            field_B[vh_B] = v_normalized;
        }
        fields_B.push_back(field_B);
    }
    return fields_B;
}

}
