/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/DualPath.hh>
#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/Types.hh>
#include <HomologyInference/VertexToPointMap.hh>

namespace HomologyInference
{

struct HomologyInferenceSettings
{
    enum Objective
    {
        FieldAlignment,
        IntegralMatching,
    };

    enum Constraints
    {
        Integer,
        Symplectic,
    };

    double w_interpolation = 1e-3;
    Objective objective = FieldAlignment;
    Constraints constraints = Symplectic;
    double var_lower_bound = -4.0;
    double var_upper_bound = 4.0;
    bool use_unconstrained_as_init = false;
    double timeout = 60.0;
};

struct HomologyInferenceResult
{
    std::string title; // An optional name to describe this result. Can be displayed in the UI.

    MatXi M_homology;
    MatXi M_cohomology;
    MatXi omega_A;
    MatXi omega_B;

    MatXd M_cohomology_unconstrained;

    std::vector<ExternalProperty<HEH, double>> gf_A;
    std::vector<ExternalProperty<HEH, double>> gf_B;
    std::vector<ExternalProperty<HEH, double>> gf_A_on_B;

    std::vector<ExternalProperty<VH, Complex>> igf_A;
    std::vector<ExternalProperty<VH, Complex>> igf_B;
    std::vector<ExternalProperty<VH, Complex>> igf_A_on_B;

    // Copies of the input data to allow visualization of this struct.
    TriMesh mesh_A;
    TriMesh mesh_B;
    PrimalLoops loops_A;
    PrimalLoops loops_B;
    VertexToPointMap vtpm;
};

// Infer a homology map from a VertexToPointMap
HomologyInferenceResult
infer_homology_map(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const PrimalLoops& _loops_A,
        const PrimalLoops& _loops_B,
        const VertexToPointMap& _vtpm,
        const HomologyInferenceSettings& _settings = HomologyInferenceSettings(),
        const std::vector<ExternalProperty<HEH, double>>& _gradient_field_A_maybe_empty = std::vector<ExternalProperty<HEH, double>>(),
        const std::vector<ExternalProperty<HEH, double>>& _gradient_field_B_maybe_empty = std::vector<ExternalProperty<HEH, double>>());

// Infer a homology map from a Functional Map
HomologyInferenceResult
infer_homology_map(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const PrimalLoops& _loops_A,
        const PrimalLoops& _loops_B,
        const MatXd& _basis_A,
        const MatXd& _basis_B,
        const MatXd& _C,
        const HomologyInferenceSettings& _settings = HomologyInferenceSettings());

MatXd
infer_homology_map_gurobi(
        const TriMesh& _mesh_B,
        const PrimalLoops& _loops_B,
        const std::vector<ExternalProperty<HEH, double>>& _mapped_fields_B,
        const std::vector<ExternalProperty<HEH, double>>& _basis_fields_B,
        const MatXi& _omega_A,
        const MatXi& _omega_B,
        const HomologyInferenceSettings& _settings = HomologyInferenceSettings(),
        const MatXd& _initial_M = MatXd());

MatXd
least_squares_field_alignment(
        const TriMesh& _mesh_B,
        const std::vector<ExternalProperty<HEH, double>>& _mapped_fields_B,
        const std::vector<ExternalProperty<HEH, double>>& _basis_fields_B);

bool
is_symplectic_homology_map(
        const MatXi& _M,
        const MatXi& _omega_A,
        const MatXi& _omega_B);

bool
is_symplectic_cohomology_map(
        const MatXi& _M,
        const MatXi& _omega_A,
        const MatXi& _omega_B);

}
