/*
 * Author: Janis Born
 */
#include "IntegerMatrix.hh"

namespace HomologyInference
{

int determinant(const MatXi& _M)
{
    ISM_ASSERT_EQ(_M.rows(), _M.cols());
    const int n = _M.rows();
    if (n == 1)
    {
        return _M(0, 0);
    }
    else if (n == 2)
    {
        return _M(0, 0) * _M(1, 1) - _M(0, 1) * _M(1, 0);
    }
    else if (n == 3)
    {
        return   _M(0,0) * _M(1,1) * _M(2,2)
               + _M(0,1) * _M(1,2) * _M(2,0)
               + _M(0,2) * _M(1,0) * _M(2,1)
               - _M(2,0) * _M(1,1) * _M(0,2)
               - _M(2,1) * _M(1,2) * _M(0,0)
               - _M(2,2) * _M(1,0) * _M(0,1);
    }
    else
    {
        // Minor expansion
        int det = 0;
        int sign = 1;
        for (int col = 0; col < n; ++col)
        {
            const int a_ij = _M(0, col);
            if (a_ij != 0)
            {
                MatXi minor_M(n - 1, n - 1);
                minor_M.leftCols(col) = _M.bottomLeftCorner(n-1, col);
                minor_M.rightCols(n-1-col) = _M.bottomRightCorner(n-1, n-1-col);
                det += sign * a_ij * determinant(minor_M);
            }
            sign *= -1;
        }
        return det;
    }
}

bool is_invertible(const MatXi& _X)
{
    const int det = determinant(_X);
    return (det != 0);
}

bool is_unimodular(const MatXi& _X)
{
    const int det = determinant(_X);
    return (det == 1 || det == -1);
}

bool is_volume_preserving(const MatXi& _X)
{
    const int det = determinant(_X);
    return (det == 1);
}

MatXi closest_integer_matrix(const MatXd& _Y)
{
    MatXi result(_Y.rows(), _Y.cols());
    for (int row = 0; row < result.rows(); ++row)
        for (int col = 0; col < result.rows(); ++col)
            result(row, col) = std::round(_Y(row, col));
    return result;
}

MatXi integer_unimodular_matrix_inverse(const MatXi& _M)
{
    ISM_ASSERT_EQ(_M.rows(), _M.cols());
    ISM_ASSERT(is_unimodular(_M));
    const MatXd M_real = _M.cast<double>();
    const MatXd M_real_inv = M_real.inverse();

    // The result should be close to an integer matrix
    const MatXd M_real_inv_rounded = M_real_inv.array().round();
    const MatXi M_inv = M_real_inv_rounded.cast<int>();
    ISM_ASSERT(M_inv * _M == MatXi::Identity(_M.rows(), _M.cols()));
    return M_inv;
}

}
