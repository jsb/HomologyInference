/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

/// True if double represents an integer up to eps
bool is_integer(
        const double _d,
        const double _eps = 0.0);

/// True if double represents an integer up to eps
bool is_integer(
        const VecXd& _vd,
        const double _eps = 0.0);

/// Convert double to integer and throw
/// if double representation was not up to eps.
int to_integer(
        const double _d,
        const double _eps = 0.0);

/// Convert double to integer and throw
/// if double representation was not up to eps.
VecXi to_integer(
        const VecXd& _vd,
        const double _eps = 0.0);

}
