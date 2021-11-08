/*
 * Author: Patrick Schmidt
 */

#include "Helpers.hh"

#include <Eigen/SparseQR>

namespace HomologyInference
{

int rank(const SparseMatrix& _M)
{
    ISM_INFO("Starting QR factorization");

    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> QR;
    QR.compute(_M);

    ISM_INFO("Finished QR factorization");

    const int r = QR.rank();

    ISM_ASSERT_EQ(_M.rows(), _M.cols());
    ISM_DEBUG_OUT("Matrix has rank n - " << _M.rows() - r);

    return r;
}

double total_area(const TriMesh &_mesh)
{
    double result = 0.0;
    for (auto fh : _mesh.faces())
        result += _mesh.calc_face_area(fh);

    return result;
}

std::pair<Vec3d, Vec3d> bounding_box(
        const TriMesh& _mesh)
{
    Vec3d pmin(+INF_DOUBLE, +INF_DOUBLE, +INF_DOUBLE);
    Vec3d pmax(-INF_DOUBLE, -INF_DOUBLE, -INF_DOUBLE);
    for (const auto& vh : _mesh.vertices())
    {
        const Vec3d& p = _mesh.point(vh);
        pmin = pmin.cwiseMin(p);
        pmax = pmax.cwiseMax(p);
    }
    return {pmin, pmax};
}

double bounding_box_diagonal(
        const TriMesh& _mesh)
{
    const auto [pmin, pmax] = bounding_box(_mesh);
    return (pmax - pmin).norm();
}

Color log_color(
        const double _val,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color)
{
    ISM_ASSERT_G(_val, 0.0);

    double lambda = (std::log(_val) - std::log(_min)) / (std::log(_max) - std::log(_min));
    lambda = std::min(std::max(lambda, 0.0), 1.0);

    return (1.0 - lambda) * _min_color + lambda * _max_color;
}

bool incident(const TriMesh& _mesh, const VH _vh, const FH _fh)
{
    ISM_ASSERT(_mesh.is_valid_handle(_fh));
    for (const auto& vh : _mesh.fv_range(_fh))
        if (vh == _vh)
            return true;
    return false;
}

VH closest_vertex(
        const TriMesh& _mesh,
        const Vec3d& _p)
{
    double closest_dist = INF_DOUBLE;
    VH closest_vh;
    for (auto v : _mesh.vertices())
    {
        const double dist = (_mesh.point(v) - _p).norm();
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_vh = v;
        }
    }

    return closest_vh;
}

}
