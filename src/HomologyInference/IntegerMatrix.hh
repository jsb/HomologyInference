/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

using IntegerMatrixPredicate = std::function<bool (const MatXi&)>;

bool is_invertible(const MatXi& _X);
bool is_unimodular(const MatXi& _X);
bool is_volume_preserving(const MatXi& _X);

/// Closest integer matrix via rounding.
MatXi closest_integer_matrix(const MatXd& _Y);

/// Compute inverse of an integer matrix _M.
/// Asserts that _M is indeed unimodular (integer invertible).
MatXi integer_unimodular_matrix_inverse(const MatXi& _M);

/// Computes the determinant of an integer matrix using recursive minor expansion.
int determinant(const MatXi& _M);

}
