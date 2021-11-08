/*
 * Author: Janis Born
 */

#include "MeshNormalization.hh"
#include <HomologyInference/Utils/Helpers.hh>

namespace HomologyInference
{

void center_mesh(
        TriMesh &_mesh,
        const TriMesh::Point &_new_origin)
{
    TriMesh::Point cog = {0,0,0};
    for (const VH vh : _mesh.vertices())
        cog += _mesh.point(vh);
    cog /= _mesh.n_vertices();
    const TriMesh::Point shift = _new_origin - cog;
    for (const VH vh : _mesh.vertices())
    {
        _mesh.point(vh) += shift;
    }

//    ISM_DEBUG_OUT(__FUNCTION__ << " moved mesh by " << shift);
}

void normalize_surface_area(
        TriMesh &_mesh,
        double _new_area)
{
    const double area = total_area(_mesh);
    const double scale = std::sqrt(_new_area) / std::sqrt(area);
    for (const VH vh : _mesh.vertices())
        _mesh.point(vh) *= scale;
    ISM_ASSERT_EPS(total_area(_mesh), _new_area, 1e-6);

//    ISM_DEBUG_OUT(__FUNCTION__ << " scaled mesh by " << scale);
}

void normalize_mesh(
        TriMesh& _mesh)
{
    center_mesh(_mesh, {0, 0, 0});
    normalize_surface_area(_mesh, 1.0);
}

Eigen::Affine3d compute_rigid_alignment(
        const MatXd &V_A, // points of A as rows
        const MatXd &V_B, // points of B as rows
        bool _allow_scaling)
{
    ISM_ASSERT_EQ(V_A.cols(), 3);
    ISM_ASSERT_EQ(V_B.cols(), 3);
    ISM_ASSERT_EQ(V_A.rows(), V_B.rows());

    Eigen::Vector3d cog_A = V_A.colwise().mean();
    Eigen::MatrixXd V_A_centered = V_A;
    V_A_centered.rowwise() -= cog_A.transpose();

    Eigen::Vector3d cog_B = V_B.colwise().mean();
    Eigen::MatrixXd V_B_centered = V_B;
    V_B_centered.rowwise() -= cog_B.transpose();

    Eigen::Matrix3d cov = V_A_centered.transpose() * V_B_centered;
    Eigen::JacobiSVD<decltype(cov)> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d rot = svd.matrixU() * svd.matrixV().transpose();
    if (rot.determinant() < 0.0)
    {
        Eigen::Vector3d diag(1.0, 1.0, -1.0);
        rot = svd.matrixU() * diag.asDiagonal() * svd.matrixV().transpose();
    }

    double scale = 1.0;
    if (_allow_scaling)
        scale = std::sqrt(V_A_centered.squaredNorm()) / std::sqrt(V_B_centered.squaredNorm());

    return Eigen::Translation3d(cog_A) * (scale * rot * Eigen::Translation3d(-cog_B));
}

Eigen::Affine3d compute_rigid_alignment(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        const TriMesh &_mesh_B,
        bool _allow_scaling)
{
    ISM_ASSERT_GEQ(_vertices_A.size(), 4);
    ISM_ASSERT_GEQ(_vertices_B.size(), 4);
    ISM_ASSERT_EQ(_vertices_A.size(), _vertices_B.size());

    Eigen::MatrixXd V_A(_vertices_A.size(), 3);
    Eigen::MatrixXd V_B(_vertices_B.size(), 3);
    for (int i = 0; i < (int)_vertices_A.size(); ++i)
    {
        V_A.row(i) = _mesh_A.point(_vertices_A[i]);
        V_B.row(i) = _mesh_B.point(_vertices_B[i]);
    }

    return compute_rigid_alignment(V_A, V_B, _allow_scaling);
}

void align_rigid(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        TriMesh &_mesh_B,
        bool _allow_scaling)
{
    Eigen::Affine3d trans = compute_rigid_alignment(_vertices_A, _vertices_B, _mesh_A, _mesh_B, _allow_scaling);

    for (const auto& vh : _mesh_B.vertices())
        _mesh_B.point(vh) = trans * _mesh_B.point(vh);

    _mesh_B.update_normals();
}

Eigen::Affine3d compute_rigid_alignment(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        bool _allow_scaling)
{
    ISM_ASSERT(_vtpm.size_okay(_mesh_A));

    std::vector<VH> samples_A;
    for (const auto& vh_A : _mesh_A.vertices())
        if (_vtpm[vh_A].is_valid())
            samples_A.push_back(vh_A);

    Eigen::MatrixXd V_A(samples_A.size(), 3);
    Eigen::MatrixXd V_B(samples_A.size(), 3);
    for (size_t i = 0; i < samples_A.size(); ++i)
    {
        const VH vh_A = samples_A[i];
        const Vec3d p_A = _mesh_A.point(vh_A);

        const BarycentricPoint& bary_B = _vtpm[vh_A];
        ISM_ASSERT(bary_B.is_valid());
        const Vec3d p_B = bary_B.point(_mesh_B);

        V_A.row(i) = p_A;
        V_B.row(i) = p_B;
    }

    return compute_rigid_alignment(V_A, V_B, _allow_scaling);
}

void align_rigid(
        const TriMesh& _mesh_A,
        TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        bool _allow_scaling)
{
    Eigen::Affine3d trans = compute_rigid_alignment(_mesh_A, _mesh_B, _vtpm, _allow_scaling);

    for (const auto& vh : _mesh_B.vertices())
        _mesh_B.point(vh) = trans * _mesh_B.point(vh);

    _mesh_B.update_normals();
}

}
