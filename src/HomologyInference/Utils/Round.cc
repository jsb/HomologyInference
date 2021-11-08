/*
 * Author: Patrick Schmidt
 */

#include <HomologyInference/Utils/CotanWeights.hh>

namespace HomologyInference
{

/// True if double exactly represents an integer
bool is_integer(
        const double _d,
        const double _eps)
{
    double intpart;
    return modf(_d, &intpart) <= _eps;
}

/// True if double exactly represents an integer
bool is_integer(
        const VecXd& _vd,
        const double _eps)
{
    for (int i = 0; i < _vd.size(); ++i)
    {
        if (!is_integer(_vd[i], _eps))
            return false;
    }

    return true;
}

int to_integer(
        const double _d,
        const double _eps)
{
    ISM_ASSERT(is_integer(_d, _eps));
    return round(_d);
}

VecXi to_integer(
        const VecXd& _vd,
        const double _eps)
{
    VecXi vi(_vd.size());

    for (int i = 0; i < _vd.size(); ++i)
        vi[i] = to_integer(_vd[i], _eps);

    return vi;
}

}
