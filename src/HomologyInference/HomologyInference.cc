/*
 * Author: Janis Born
 */
#include "HomologyInference.hh"

#include <HomologyInference/CohomologyBasis.hh>
#include <HomologyInference/FunctionalMaps.hh>
#include <HomologyInference/Genus.hh>
#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/IntegerMatrix.hh>
#include <HomologyInference/Utils/Timer.hh>

#include <gurobi_c++.h>

#include <queue>

namespace HomologyInference
{

HomologyInferenceResult
infer_homology_map(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const PrimalLoops& _loops_A,
        const PrimalLoops& _loops_B,
        const VertexToPointMap& _vtpm,
        const HomologyInferenceSettings& _settings,
        const std::vector<ExternalProperty<HEH, double>>& _gradient_field_A_maybe_empty,
        const std::vector<ExternalProperty<HEH, double>>& _gradient_field_B_maybe_empty)
{
    Timer timer(__FUNCTION__);

    HomologyInferenceResult result;
    // Copies of the input data.
    result.mesh_A = _mesh_A;
    result.mesh_B = _mesh_B;
    result.loops_A = _loops_A;
    result.loops_B = _loops_B;
    result.vtpm = _vtpm;

    ISM_ASSERT(_vtpm.size_okay(_mesh_A));
    ISM_ASSERT_EQ(genus(_mesh_A), genus(_mesh_B));
    const int g = genus(_mesh_A);

    // Basis loops
    ISM_ASSERT_EQ((int)_loops_A.size(), 2*g);
    ISM_ASSERT_EQ((int)_loops_B.size(), 2*g);

    // Conformal gradient fields
    result.gf_A = _gradient_field_A_maybe_empty.empty() ?
            cohomology_basis(_loops_A, _mesh_A)
          : _gradient_field_A_maybe_empty;
    result.gf_B = _gradient_field_B_maybe_empty.empty() ?
            cohomology_basis(_loops_B, _mesh_B)
          : _gradient_field_B_maybe_empty;
    ISM_ASSERT_EQ((int)result.gf_A.size(), 2*g);
    ISM_ASSERT_EQ((int)result.gf_B.size(), 2*g);

    // Integrated fields on A
    result.igf_A = integrated_fields_complex(result.gf_A, _mesh_A);

    // Map and interpolate integrated fields from A onto B
    result.igf_A_on_B = interpolate_fields(_mesh_A, _mesh_B, _vtpm, result.igf_A, _settings.w_interpolation);

    // Differentiate on B to get gradient fields again
    result.gf_A_on_B = differentiate_fields(result.igf_A_on_B, _mesh_B);

    // Compute an unconstrained solution.
    if (_settings.objective == HomologyInferenceSettings::FieldAlignment)
        result.M_cohomology_unconstrained = least_squares_field_alignment(_mesh_B, result.gf_A_on_B, result.gf_B);
    if (_settings.objective == HomologyInferenceSettings::IntegralMatching)
        result.M_cohomology_unconstrained = integrate_loops(result.gf_A_on_B, _loops_B);

    // Compute matching of gradient fields on B
    result.omega_A = intersection_form(_mesh_A, _loops_A);
    result.omega_B = intersection_form(_mesh_B, _loops_B);

    MatXd M_gurobi_result = infer_homology_map_gurobi(
                _mesh_B,
                _loops_B,
                result.gf_A_on_B,
                result.gf_B,
                result.omega_A,
                result.omega_B,
                _settings,
                _settings.use_unconstrained_as_init ? result.M_cohomology_unconstrained : MatXd());
    result.M_homology = closest_integer_matrix(M_gurobi_result);

    result.M_cohomology = integer_unimodular_matrix_inverse(result.M_homology).transpose();

    if (_settings.constraints == HomologyInferenceSettings::Symplectic)
    {
        ISM_ASSERT(is_symplectic_homology_map(result.M_homology, result.omega_A, result.omega_B));
        ISM_ASSERT(is_symplectic_cohomology_map(result.M_cohomology, result.omega_A, result.omega_B));
    }

    return result;
}

HomologyInferenceResult
infer_homology_map(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const PrimalLoops& _loops_A,
        const PrimalLoops& _loops_B,
        const MatXd& _basis_A,
        const MatXd& _basis_B,
        const MatXd& _C,
        const HomologyInferenceSettings& _settings)
{
    Timer timer(__FUNCTION__);

    HomologyInferenceResult result;
    // Copies of the input data.
    result.mesh_A = _mesh_A;
    result.mesh_B = _mesh_B;
    result.loops_A = _loops_A;
    result.loops_B = _loops_B;
    result.vtpm.init(_mesh_A, BarycentricPoint()); // "Empty"

    ISM_ASSERT_EQ(genus(_mesh_A), genus(_mesh_B));
    const int g = genus(_mesh_A);

    // Basis loops
    ISM_ASSERT_EQ(_loops_A.size(), 2*g);
    ISM_ASSERT_EQ(_loops_B.size(), 2*g);

    // Conformal gradient fields
    result.gf_A = cohomology_basis(_loops_A, _mesh_A);
    result.gf_B = cohomology_basis(_loops_B, _mesh_B);
    ISM_ASSERT_EQ(result.gf_A.size(), 2*g);
    ISM_ASSERT_EQ(result.gf_B.size(), 2*g);

    // Integrated fields on A
    result.igf_A = integrated_fields_complex(result.gf_A, _mesh_A);
    result.igf_B = integrated_fields_complex(result.gf_B, _mesh_B);

    // Map and integrated fields from A onto B using functional map
    result.igf_A_on_B = map_fields_fmap(_mesh_A, _mesh_B, _basis_A, _basis_B, _C, result.igf_A);

    // Differentiate on B to get gradient fields again
    result.gf_A_on_B = differentiate_fields(result.igf_A_on_B, _mesh_B);

    // Compute an unconstrained least-squares matching.
    if (_settings.objective == HomologyInferenceSettings::FieldAlignment)
        result.M_cohomology_unconstrained = least_squares_field_alignment(_mesh_B, result.gf_A_on_B, result.gf_B);
    if (_settings.objective == HomologyInferenceSettings::IntegralMatching)
        result.M_cohomology_unconstrained = integrate_loops(result.gf_A_on_B, _loops_B);

    // Compute matching of gradient fields on B
    result.omega_A = intersection_form(_mesh_A, _loops_A);
    result.omega_B = intersection_form(_mesh_B, _loops_B);

    MatXd M_gurobi_result = infer_homology_map_gurobi(
                _mesh_B,
                _loops_B,
                result.gf_A_on_B,
                result.gf_B,
                result.omega_A,
                result.omega_B,
                _settings,
                _settings.use_unconstrained_as_init ? result.M_cohomology_unconstrained : MatXd());
    result.M_homology = closest_integer_matrix(M_gurobi_result);

    result.M_cohomology = integer_unimodular_matrix_inverse(result.M_homology).transpose();

    if (_settings.constraints == HomologyInferenceSettings::Symplectic)
    {
        ISM_ASSERT(is_symplectic_homology_map(result.M_homology, result.omega_A, result.omega_B));
        ISM_ASSERT(is_symplectic_cohomology_map(result.M_cohomology, result.omega_A, result.omega_B));
    }

    return result;
}

MatXd
least_squares_field_alignment(
        const TriMesh& _mesh_B,
        const std::vector<ExternalProperty<HEH, double>>& _mapped_fields_B,
        const std::vector<ExternalProperty<HEH, double>>& _basis_fields_B)
{
    Timer timer(__FUNCTION__);

    const int g = genus(_mesh_B);
    MatXd M_real(2*g, 2*g);

    // Solve via Moore-Penrose pseudoinverse
    MatXd Phi_X_A(_mesh_B.n_edges(), 2*g);
    MatXd X_B(_mesh_B.n_edges(), 2*g);
    for (const auto& eh : _mesh_B.edges())
    {
        const int row = eh.idx();
        const auto& heh = _mesh_B.halfedge_handle(eh, 0);
        for (int generator = 0; generator < 2*g; ++generator)
        {
            const int col = generator;
            Phi_X_A(row, col) = _mapped_fields_B[generator][heh];
            X_B(row, col) = _basis_fields_B[generator][heh];
        }
    }
    auto X_B_pseudoinv = (X_B.transpose() * X_B).inverse() * X_B.transpose();
    M_real = X_B_pseudoinv * Phi_X_A;
    ISM_ASSERT_EQ(M_real.rows(), 2*g);
    ISM_ASSERT_EQ(M_real.cols(), 2*g);

    return M_real;
}

MatXd
infer_homology_map_gurobi(
        const TriMesh& _mesh_B,
        const PrimalLoops& _loops_B,
        const std::vector<ExternalProperty<HEH, double>>& _mapped_fields_B,
        const std::vector<ExternalProperty<HEH, double>>& _basis_fields_B,
        const MatXi& _omega_A,
        const MatXi& _omega_B,
        const HomologyInferenceSettings& _settings,
        const MatXd& _initial_M)
{
    Timer timer(__FUNCTION__);

    const int g = genus(_mesh_B);
    ISM_ASSERT_EQ((int)_loops_B.size(), 2*g);
    ISM_ASSERT_EQ((int)_mapped_fields_B.size(), 2*g);
    ISM_ASSERT_EQ((int)_basis_fields_B.size(), 2*g);
    ISM_ASSERT_EQ(_omega_A.rows(), 2*g);
    ISM_ASSERT_EQ(_omega_A.cols(), 2*g);
    ISM_ASSERT_EQ(_omega_B.rows(), 2*g);
    ISM_ASSERT_EQ(_omega_B.cols(), 2*g);
    ISM_ASSERT(_initial_M.size() == 0 || (_initial_M.rows() == 2*g && _initial_M.cols() == 2*g));

    MatXd cohomology_matching(2*g, 2*g);

    try
    {
        GRBEnv env;
        GRBModel m(env);
        if (_settings.timeout > 0.0)
            m.set(GRB_DoubleParam_TimeLimit, _settings.timeout);

        // Create variables
        Eigen::Array<GRBVar, Eigen::Dynamic, Eigen::Dynamic> var_M(2*g, 2*g);
        for (int row = 0; row < 2*g; ++row)
        {
            for (int col = 0; col < 2*g; ++col)
            {
                const std::string varname = "M(" + std::to_string(row) + "," + std::to_string(col) + ")";
                var_M(row, col) = m.addVar(_settings.var_lower_bound, _settings.var_upper_bound, 0.0, GRB_INTEGER, varname);

                // Set initial value if provied
                if (_initial_M.size() != 0)
                    var_M(row, col).set(GRB_DoubleAttr_Start, _initial_M(row, col));
            }
        }

        if (_settings.objective == HomologyInferenceSettings::FieldAlignment)
        {
            // Objective
            //     || Phi_X_A - X_B M ||_F
            MatXd Phi_X_A(_mesh_B.n_edges(), 2*g);
            MatXd X_B(_mesh_B.n_edges(), 2*g);
            for (const auto& eh : _mesh_B.edges())
            {
                const int row = eh.idx();
                const auto& heh = _mesh_B.halfedge_handle(eh, 0);
                for (int generator = 0; generator < 2*g; ++generator)
                {
                    const int col = generator;
                    Phi_X_A(row, col) = _mapped_fields_B[generator][heh];
                    X_B(row, col) = _basis_fields_B[generator][heh];
                }
            }
            GRBQuadExpr objective = 0;
            for (int row = 0; row < _mesh_B.n_edges(); ++row)
            {
                for (int col = 0; col < 2*g; ++col)
                {
                    GRBLinExpr entry_ij = Phi_X_A(row, col);
                    for (int k = 0; k < 2*g; ++k)
                    {
                        entry_ij -= X_B(row, k) * var_M(k, col);
                    }
                    objective += entry_ij * entry_ij;
                }
            }
            m.setObjective(objective);
        }
        if (_settings.objective == HomologyInferenceSettings::IntegralMatching)
        {
            MatXd M_integrated = integrate_loops(_mapped_fields_B, _loops_B);
            GRBQuadExpr objective = 0;
            for (int row = 0; row < 2*g; ++row)
            {
                for (int col = 0; col < 2*g; ++col)
                {
                    GRBLinExpr entry_ij = M_integrated(row, col) - var_M(row, col);
                    objective += entry_ij * entry_ij;
                }
            }
            m.setObjective(objective);
        }

        // Constraints:
        //     M * Omega_A * M^T == Omega_B
        // The matrix M we compute here is a cohomology matching,
        // so we have the "cohomology symplectic constraint" as above.
        m.set(GRB_IntParam_NonConvex, 2);
        for (int row = 0; row < 2*g; ++row)
        {
            for (int col = row + 1; col < 2*g; ++col)
            {
                GRBQuadExpr M_omega_A_Mt_row_col = 0.0;
                for (int i = 0; i < 2*g; ++i)
                    for (int j = 0; j < 2*g; ++j)
                        if (_omega_A(i, j) != 0)
                            M_omega_A_Mt_row_col += var_M(row, i) * _omega_A(i, j) * var_M(col, j);
                const int omega_B_row_col = _omega_B(row, col);
                m.addQConstr(M_omega_A_Mt_row_col == omega_B_row_col);
            }
        }

        // Solve
        m.optimize();
        int status = m.get(GRB_IntAttr_Status);
        ISM_DEBUG_VAR(status);

        // Extract solution
        for (int row = 0; row < 2*g; ++row)
            for (int col = 0; col < 2*g; ++col)
                cohomology_matching(row, col) = var_M(row, col).get(GRB_DoubleAttr_X);
    }
    catch(GRBException e)
    {
        ISM_ERROR("Gurobi error code: " << e.getErrorCode());
        ISM_ERROR_throw("Gurobi error: " << e.getMessage());
    }

    // Gurobi's solution is a cohomology matching.
    // Convert to a homology matching (dual basis) via inverse transpose.
    MatXd homology_matching = cohomology_matching.inverse().transpose();

    MatXi cohomology_matching_i = closest_integer_matrix(cohomology_matching);
    MatXi homology_matching_i = closest_integer_matrix(homology_matching);
    ISM_ASSERT(is_symplectic_cohomology_map(cohomology_matching_i, _omega_A, _omega_B));
    ISM_ASSERT(is_symplectic_homology_map(homology_matching_i, _omega_A, _omega_B));

    return homology_matching;
}

bool
is_symplectic_homology_map(
        const MatXi& _M,
        const MatXi& _omega_A,
        const MatXi& _omega_B)
{
    return (_omega_A == _M.transpose() * _omega_B * _M);
}

bool
is_symplectic_cohomology_map(
        const MatXi& _M,
        const MatXi& _omega_A,
        const MatXi& _omega_B)
{
    return (_M * _omega_A * _M.transpose() == _omega_B);
}

}
