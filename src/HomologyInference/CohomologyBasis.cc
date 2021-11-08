/*
 * Author: Patrick Schmidt
 */

#include "CohomologyBasis.hh"

#include <HomologyInference/Genus.hh>
#include <HomologyInference/Utils/BarycentricPoint.hh>
#include <HomologyInference/Utils/CotanWeights.hh>
#include <HomologyInference/Utils/Timer.hh>
#include <Eigen/SparseCholesky>
#include <queue>

namespace HomologyInference
{

std::vector<ExternalProperty<HEH, double>>
cohomology_basis(
        const PrimalLoops& _loops,
        const TriMesh& _mesh)
{
    Timer timer(__FUNCTION__);

    // Set up linear system
    const int n_fields = _loops.size();
    const int n_constr = _mesh.n_faces() + _mesh.n_vertices() + n_fields;
    const int n_edges = _mesh.n_edges();
    std::vector<Triplet> triplets;
    MatXd rhs = MatXd::Zero(n_constr, n_fields);

    // 1 if canonical direction, -1 otherwise
    auto he_sign = [&] (const HEH& heh)
    {
        const auto& eh = _mesh.edge_handle(heh);
        const auto& heh0 = _mesh.halfedge_handle(eh, 0);
        return (heh == heh0) ? 1.0 : -1.0;
    };

    // Closedness condition:
    // Integral along every contractible loop is zero.
    for (auto f : _mesh.faces())
    {
        for (auto h : f.halfedges())
            triplets.push_back(Triplet(f.idx(), h.edge().idx(), he_sign(h)));
    }

    // Harmonicity condition:
    // Laplace of field is zero at every vertex
    int offset = _mesh.n_faces();
    for (auto v : _mesh.vertices())
    {
        if (v.is_boundary())
        {
            for (auto h : v.outgoing_halfedges())
                if (h.edge().is_boundary())
                    triplets.push_back(Triplet(offset + v.idx(), h.edge().idx(), he_sign(h) * 1.0 / _mesh.calc_edge_length(h.edge())));
        }
        else
        {
            for (auto h : v.outgoing_halfedges())
                triplets.push_back(Triplet(offset + v.idx(), h.edge().idx(), he_sign(h) * cotan_weight(_mesh, h)));
        }
    }

    // Duality condition:
    // Integral of field i along loop i is 1. All others are 0.
    offset += _mesh.n_vertices();
    for (int i = 0; i < _loops.size(); ++i)
    {
        for (auto heh : _loops[i].hehs)
        {
            const auto& eh = _mesh.edge_handle(heh);
            triplets.push_back(Triplet(offset + i, eh.idx(), he_sign(heh)));
        }
    }
    rhs.block(offset, 0, n_fields, n_fields) = MatXd::Identity(n_fields, n_fields);

    SparseMatrix A(n_constr, n_edges);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Solve A^T*A * x = A^T*b via sparse Cholesky
    const SparseMatrix AtA = A.transpose() * A;
    const MatXd Atrhs = A.transpose() * rhs;
    Eigen::SimplicialLDLT<SparseMatrix> solver;
    solver.compute(AtA);
    ISM_ASSERT(solver.info() == Eigen::Success);
    MatXd x = solver.solve(Atrhs);
    ISM_ASSERT(solver.info() == Eigen::Success);
    ISM_ASSERT_EQ(x.rows(), n_edges);
    ISM_ASSERT_EQ(x.cols(), n_fields);

    double residual = (A * x - rhs).norm();
    ISM_DEBUG_VAR(residual);

    // Convert to halfedge properties
    auto gradients = std::vector<ExternalProperty<HEH, double>>(n_fields, ExternalProperty<HEH, double>(_mesh));
    for (int i = 0; i < n_fields; ++i)
        for (auto h : _mesh.halfedges())
            gradients[i][h] = he_sign(h) * x(h.edge().idx(), i);

    return gradients;
}

ExternalProperty<VH, double>
integrate_field_real(
        const ExternalProperty<HEH, double>& _gradient,
        const TriMesh& _mesh,
        const VH _seed_vh,
        const double _seed_value)
{
    // Flood fill starting from first vertex
    ExternalProperty<VH, double> u(_mesh, NAN_DOUBLE); // integrated function
    ExternalProperty<VH, bool> visited(_mesh, false);

    std::queue<SHEH> queue;
    for (auto h : _mesh.voh_range(_seed_vh))
        queue.push(h);
    u[_seed_vh] = _seed_value;
    visited[_seed_vh] = true;

    while (!queue.empty())
    {
        const auto h = queue.front();
        queue.pop();

        if (visited[h.to()])
            continue;

        u[h.to()] = u[h.from()] + _gradient[h];
        visited[h.to()] = true;

        for (auto h_enq : h.to().outgoing_halfedges())
        {
            if (!visited[h_enq.to()])
                queue.push(h_enq);
        }
    }

    return u;
}

std::vector<ExternalProperty<VH, double>>
integrated_fields_real(
        const std::vector<ExternalProperty<HEH, double>>& _gradients,
        const TriMesh& _mesh,
        const VH _seed_vh,
        const double _seed_value)
{
    std::vector<ExternalProperty<VH, double>> result;
    for (const auto& g : _gradients)
        result.push_back(integrate_field_real(g, _mesh, _seed_vh, _seed_value));
    return result;
}

ExternalProperty<VH, Complex>
integrate_field_complex(
        const ExternalProperty<HEH, double>& _gradient,
        const TriMesh& _mesh,
        const VH _seed_vh,
        const Complex _seed_value)
{
    ExternalProperty<VH, double> u_real = integrate_field_real(_gradient, _mesh, _seed_vh, 0.0);
    ExternalProperty<VH, Complex> u_complex(_mesh, NAN_DOUBLE);
    for (const auto& v : _mesh.vertices())
    {
        const double angle = u_real[v] * 2 * M_PI;
        const Complex rot = std::polar(1.0, angle);
        u_complex[v] = rot * _seed_value;
    }
    return u_complex;
}

std::vector<ExternalProperty<VH, Complex>>
integrated_fields_complex(
        const std::vector<ExternalProperty<HEH, double>>& _gradients,
        const TriMesh& _mesh,
        const VH _seed_vh,
        const Complex _seed_value)
{
    Timer timer(__FUNCTION__);

    std::vector<ExternalProperty<VH, Complex>> result;
    for (const auto& g : _gradients)
        result.push_back(integrate_field_complex(g, _mesh, _seed_vh, _seed_value));
    return result;
}

ExternalProperty<HEH, double>
differentiate_field(
        const ExternalProperty<VH, Complex>& _u,
        const TriMesh& _mesh)
{
    ExternalProperty<HEH, double> gradient(_mesh);
    for (const auto& heh : _mesh.halfedges())
    {
        const VH vh0 = _mesh.from_vertex_handle(heh);
        const VH vh1 = _mesh.to_vertex_handle(heh);
        const Complex u0 = _u[vh0];
        const Complex u1 = _u[vh1];
        const Complex rot = u1 / u0;
        const double angle_diff = std::arg(rot) / (2 * M_PI);
        gradient[heh] = angle_diff;
    }
    return gradient;
}

std::vector<ExternalProperty<HEH, double>>
differentiate_fields(
        const std::vector<ExternalProperty<VH, Complex>>& _us,
        const TriMesh& _mesh)
{
    Timer timer(__FUNCTION__);

    std::vector<ExternalProperty<HEH, double>> result;
    for (const auto& u : _us)
        result.push_back(differentiate_field(u, _mesh));
    return result;
}

double
integrate_loop(
        const ExternalProperty<HEH, double>& _gradient,
        const PrimalLoop& _loop)
{
    double res = 0.0;
    for (auto h : _loop.hehs)
        res += _gradient[h];
    return res;
}

MatXd
integrate_loops(
        const std::vector<ExternalProperty<HEH, double>>& _fields,
        const PrimalLoops& _loops)
{
    ISM_ASSERT(_loops.size() == _fields.size());
    const int n = _loops.size();
    MatXd M_integrated = MatXd::Zero(n, n);
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < n; ++col)
            M_integrated(row, col) = integrate_loop(_fields[col], _loops[row]);
    return M_integrated;
}

ExternalProperty<HEH, double>
combine_fields(
        const TriMesh& _mesh,
        const VecXi& _coeffs,
        const std::vector<ExternalProperty<HEH, double>>& _fields)
{
    const VecXd coeffs_d = _coeffs.cast<double>();
    return combine_fields(_mesh, coeffs_d, _fields);
}

ExternalProperty<HEH, double>
combine_fields(
        const TriMesh& _mesh,
        const VecXd& _coeffs,
        const std::vector<ExternalProperty<HEH, double>>& _fields)
{
    ISM_ASSERT_EQ(_coeffs.size(), _fields.size());
    ExternalProperty<HEH, double> field(_mesh, 0.0);
    for (int i = 0; i < _coeffs.size(); ++i)
    {
        const auto& input_field = _fields[i];
        ISM_ASSERT(input_field.size_okay(_mesh));
        for (const auto& heh : _mesh.halfedges())
            field[heh] += _coeffs[i] * input_field[heh];
    }
    return field;
}

std::vector<ExternalProperty<HEH, double>>
transform_fields(
        const TriMesh& _mesh,
        const MatXi& _M,
        const std::vector<ExternalProperty<HEH, double>>& _fields)
{
    ISM_ASSERT_EQ(_M.cols(), _fields.size());
    ISM_ASSERT(!_fields.empty());
    std::vector<ExternalProperty<HEH, double>> result;
    for (int row = 0; row < _M.rows(); ++row)
    {
        const VecXi M_row = _M.row(row);
        ExternalProperty<HEH, double> field = combine_fields(_mesh, M_row, _fields);
        result.push_back(field);
    }
    return result;
}

}
