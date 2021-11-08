/*
 * Author: Janis Born
 */

#include "VertexToPointMap.hh"

#include <HomologyInference/Extern/ACG_BSP/Geometry/Algorithms.hh>
#include <HomologyInference/Extern/ACG_BSP/Geometry/bsp/TriangleBSPT.hh>
#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/Utils/CotanWeights.hh>

#include <Eigen/SparseCholesky>
#include <memory>

namespace HomologyInference
{

using BSP = OpenMeshTriangleBSPT<TriMesh>;
using BSPptr = std::unique_ptr<BSP>;

BSPptr build_bsp(
        const TriMesh& _mesh)
{
    BSPptr bsp(new BSP(_mesh));
    bsp->reserve(_mesh.n_faces());
    for (const FH fh : _mesh.faces())
        bsp->push_back(fh);
    bsp->build(10, 100); // max vertices per leaf 10, max depth 100

    return bsp;
}

bool is_dense(const VertexToPointMap& _vtpm)
{
    for (const auto& bary : _vtpm.container_)
        if (!bary.is_valid())
            return false;
    return true;
}

VertexToPointMap generate_vtpm_by_projection(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B)
{
    BSPptr bspt_B = build_bsp(_mesh_B);

    VertexToPointMap vtpm(_mesh_A);
    for (const auto& vh_A : _mesh_A.vertices())
    {
        // Query point
        const auto p_A = to_acg3(_mesh_A.point(vh_A));
        // Closest triangle
        const FH fh_B = bspt_B->nearest(p_A).handle;
        ISM_ASSERT(fh_B.is_valid());
        const auto c_B = _mesh_B.fv_range(fh_B).to_array<3>([&](auto _vh_B){ return to_acg3(_mesh_B.point(_vh_B)); });
        // Point on triangle
        const Vec3d p_B = to_eigen(ACG::Geometry::closestPointTri(p_A, c_B[0], c_B[1], c_B[2]));
        // Barycentric point
        vtpm[vh_A] = BarycentricPoint(p_B, fh_B, _mesh_B);
        ISM_ASSERT_FINITE(vtpm[vh_A].alpha());
        ISM_ASSERT_FINITE(vtpm[vh_A].beta());
    }
    return vtpm;
}

VertexToPointMap
generate_vtpm_from_landmarks(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const std::vector<VH>& _landmarks_A,
        const std::vector<VH>& _landmarks_B)
{
    ISM_ASSERT_EQ(_landmarks_A.size(), _landmarks_B.size());
    VertexToPointMap vtpm(_mesh_A, BarycentricPoint());
    const size_t n = _landmarks_A.size();
    for (size_t i = 0; i < n; ++i)
    {
        const VH& vh_A = _landmarks_A[i];
        const VH& vh_B = _landmarks_B[i];
        ISM_ASSERT(_mesh_A.is_valid_handle(vh_A));
        ISM_ASSERT(_mesh_B.is_valid_handle(vh_B));

        BarycentricPoint bary_B = BarycentricPoint(vh_B, _mesh_B);
        vtpm[vh_A] = bary_B;
    }
    return vtpm;
}

VertexToPointMap generate_vtpm_from_overlay_meshes(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const TriMesh& _overlay_on_A,
        const TriMesh& _overlay_on_B)
{
    ISM_ASSERT_EQ(_overlay_on_A.n_vertices(), _overlay_on_B.n_vertices());

    // Build bsp trees for overlay on A
    // and mesh B
    BSPptr bsp_overlay_A = build_bsp(_overlay_on_A);
    BSPptr bsp_B = build_bsp(_mesh_B);

    VertexToPointMap vtpm(_mesh_A);
    for (const auto v : _mesh_A.vertices())
    {
        // Express vertex via barycentric coordiantes on overlay A
        const auto p_source = _mesh_A.point(v);
        const auto nn_overlay = bsp_overlay_A->nearest(to_acg3(p_source));
        ISM_ASSERT(nn_overlay.handle.is_valid());
        ISM_ASSERT(_overlay_on_A.is_valid_handle(nn_overlay.handle));
        ISM_ASSERT_L(nn_overlay.dist, 1e-4);
        const auto bary = BarycentricPoint(p_source, nn_overlay.handle, _overlay_on_A);

        // Evaluate barycentric coordinates on overlay B
        auto p_target = bary.point(_overlay_on_B);
        ISM_ASSERT_FINITE_MAT(p_target);

        // Express as barycentric point on B
        const auto nn_B = bsp_B->nearest(to_acg3(p_target));
        ISM_ASSERT(nn_B.handle.is_valid());
        ISM_ASSERT(_overlay_on_A.is_valid_handle(nn_B.handle));
        ISM_ASSERT_L(nn_B.dist, 1e-4);
        vtpm[v] = BarycentricPoint(p_target, nn_B.handle, _mesh_B);
    }

    return vtpm;
}

ExternalProperty<VH, Complex> interpolate_field(const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        const ExternalProperty<VH, Complex>& _field_A,
        const double _w_interp)
{
    ISM_ASSERT(_vtpm.size_okay(_mesh_A));
    ISM_ASSERT(_field_A.size_okay(_mesh_A));

    // Setup constraint system:
    // diag(w) * A * x = diag(w) * b
    const int constr_rows = _mesh_A.n_vertices() + _mesh_B.n_vertices();
    const int constr_cols = _mesh_B.n_vertices();
    SparseMatrix A(constr_rows, constr_cols);
    VecXd w = VecXd::Zero(constr_rows);
    VecX<Complex> b = VecX<Complex>::Zero(constr_rows);
    {
        std::vector<Triplet> A_triplets;
        int row = 0;

        for (const auto& vh_A : _mesh_A.vertices())
        {
            const BarycentricPoint& bary = _vtpm[vh_A];
            // In a sparse VTPM, vertices of A may have no corresponding point on B.
            if (!bary.is_valid())
                continue;

            const Complex& value = _field_A[vh_A];

            A_triplets.emplace_back(row, bary.vh_a(_mesh_B).idx(), bary.alpha());
            A_triplets.emplace_back(row, bary.vh_b(_mesh_B).idx(), bary.beta());
            A_triplets.emplace_back(row, bary.vh_c(_mesh_B).idx(), bary.gamma());
            b[row] = value;
            w[row] = _w_interp;

            ++row;
        }

        for (const auto& vh_B : _mesh_B.vertices())
        {
            double area = 0.0;
            double weight_sum = 0.0;
            for (const auto& heh_B : _mesh_B.voh_range(vh_B)) {
                const VH& vh_to_B = _mesh_B.to_vertex_handle(heh_B);
                const double weight = cotan_weight(_mesh_B, heh_B);

                A_triplets.emplace_back(row, vh_to_B.idx(), weight);

                const HEH& heh_prev_B = _mesh_B.prev_halfedge_handle(heh_B);
                area += _mesh_B.calc_sector_area(heh_prev_B);
                weight_sum += weight;
            }

            A_triplets.emplace_back(row, vh_B.idx(), -weight_sum);
            b[row] = 0;
            w[row] = area / 3;

            ++row;
        }

        A.setFromTriplets(A_triplets.cbegin(), A_triplets.cend());
    }

    SparseMatrix AtWA = A.transpose() * w.asDiagonal() * A;
    VecX<Complex> AtWb = A.transpose() * w.asDiagonal() * b;

    Eigen::SimplicialLDLT<SparseMatrix> solver(AtWA);
    VecX<Complex> x = solver.solve(AtWb);
    ISM_ASSERT_EQ(x.rows(), _mesh_B.n_vertices());

    ExternalProperty<VH, Complex> field_B(_mesh_B);
    for (const auto& vh_B : _mesh_B.vertices())
    {
        field_B[vh_B] = x[vh_B.idx()] / std::abs(x[vh_B.idx()]);
    }
    return field_B;
}

std::vector<ExternalProperty<VH, Complex>>
interpolate_fields(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        const std::vector<ExternalProperty<VH, Complex>>& _fields_A,
        const double _w_interp)
{
    std::vector<ExternalProperty<VH, Complex>> result;
    for (const auto& f_A : _fields_A)
        result.push_back(interpolate_field(_mesh_A, _mesh_B, _vtpm, f_A, _w_interp));
    return result;
}

void sanitize_vtpm(
        VertexToPointMap& _vtpm,
        const double _epsilon)
{
    const double scale = 1.0 - _epsilon;
    const double o_alpha = 1.0 / 3.0;
    const double o_beta = 1.0 / 3.0;
    for (auto& bary : _vtpm.container_)
    {
        if (!std::isfinite(bary.alpha()) || !std::isfinite(bary.beta()))
            bary.set_alpha_beta(o_alpha, o_beta);
        int tries = 0;
        while (!bary.is_inside_exclusive() && tries < 10000)
        {
            bary.set_alpha_beta(
                (bary.alpha() - o_alpha) * scale + o_alpha,
                (bary.beta()  - o_beta)  * scale + o_beta
            );
            ++tries;
        }
        if (!bary.is_inside_exclusive())
        {
            bary.set_alpha_beta(o_alpha, o_beta);
        }
    }
}

void randomly_permute_vtpm(
        VertexToPointMap& _vtpm,
        int _n_transpositions)
{
    srand(0);

    ISM_ASSERT_GEQ(_vtpm.container_.size(), 2);
    int n_transposed = 0;
    while (n_transposed < _n_transpositions) {
        const size_t i0 = rand() % _vtpm.container_.size();
        const size_t i1 = rand() % _vtpm.container_.size();
        if (i0 != i1)
        {
            std::swap(_vtpm.container_[i0], _vtpm.container_[i1]);
            ++n_transposed;
        }
    }
}

void randomly_permute_vtpm(
        VertexToPointMap& _vtpm,
        double _ratio)
{
    _vtpm = randomly_permuted_vtpm(_vtpm, _ratio);
}

VertexToPointMap randomly_permuted_vtpm(
        const VertexToPointMap& _vtpm,
        int _n_transpositions)
{
    VertexToPointMap result = _vtpm;
    randomly_permute_vtpm(result, _n_transpositions);
    return result;
}

VertexToPointMap randomly_permuted_vtpm(
        const VertexToPointMap& _vtpm,
        double _ratio)
{
    srand(0);

    std::vector<size_t> orig_indices(_vtpm.size());
    std::iota(orig_indices.begin(), orig_indices.end(), 0);

    // Only keep a (random) n-element subset of the original indices.
    std::random_shuffle(orig_indices.begin(), orig_indices.end());
    const size_t n = std::round(_vtpm.size() * _ratio);
    orig_indices.resize(n);

    // Randomly permute whithin this subset.
    std::vector<size_t> new_indices = orig_indices;
    std::random_shuffle(new_indices.begin(), new_indices.end());

    VertexToPointMap result = _vtpm;
    for (size_t i = 0; i < n; ++i)
    {
        const size_t orig_index = orig_indices[i];
        const size_t new_index = new_indices[i];
        result.container_[new_index] = _vtpm.container_[orig_index];
    }
    return result;
}

void subsample_vtpm(
        VertexToPointMap& _vtpm,
        double _ratio)
{
    _vtpm = subsampled_vtpm(_vtpm, _ratio);
}

void subsample_vtpm(
        VertexToPointMap& _vtpm,
        int _n_keep)
{
    _vtpm = subsampled_vtpm(_vtpm, _n_keep);
}

VertexToPointMap subsampled_vtpm(
        const VertexToPointMap& _vtpm,
        double _ratio)
{
    const int n = std::round(_vtpm.size() * _ratio);
    return subsampled_vtpm(_vtpm, n);
}

VertexToPointMap subsampled_vtpm(
        const VertexToPointMap& _vtpm,
        int _n_keep)
{
    srand(684);

    std::vector<size_t> orig_indices(_vtpm.size());
    std::iota(orig_indices.begin(), orig_indices.end(), 0);

    // Only keep a (random) n-element subset of the original indices.
    std::random_shuffle(orig_indices.begin(), orig_indices.end());
    orig_indices.resize(_n_keep);

    // Delete all others samples
    VertexToPointMap result(_vtpm.size());
    for (int i = 0; i < _n_keep; ++i)
    {
        const size_t orig_index = orig_indices[i];
        result.container_[orig_index] = _vtpm.container_[orig_index];
    }

    return result;
}

}
