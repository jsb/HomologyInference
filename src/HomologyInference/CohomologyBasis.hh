/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/Types.hh>

namespace HomologyInference
{

/**
 * Compute 2g closed 1-forms (cocycles), represented as scalar values per directed edge.
 * Integrating form i along a loop in homotopy class i yields 1. (0 for all others.)
 * Gu [2003] Global Conformal Surface Parametrization.
 */
std::vector<ExternalProperty<HEH, double>>
cohomology_basis(
        const PrimalLoops& _loops,
        const TriMesh& _mesh);

/**
 * Integrate a closed 1-form over surface via flood-fill.
 * Creates a (discontinuous) real-valued scalar field with period 1.
 */
ExternalProperty<VH, double>
integrate_field_real(
        const ExternalProperty<HEH, double>& _gradient,
        const TriMesh& _mesh,
        const VH _seed_vh = VH(0),
        const double _seed_value = 0.0);

std::vector<ExternalProperty<VH, double>>
integrated_fields_real(
        const std::vector<ExternalProperty<HEH, double>>& _gradients,
        const TriMesh& _mesh,
        const VH _seed_vh = VH(0),
        const double _seed_value = 0.0);

/**
 * Integrate a closed 1-form over surface via flood-fill.
 * Creates a (continuous) complex-valued scalar field representing angles.
 */
ExternalProperty<VH, Complex>
integrate_field_complex(
        const ExternalProperty<HEH, double>& _gradient,
        const TriMesh& _mesh,
        const VH _seed_vh = VH(0),
        const Complex _seed_value = Complex(1.0));

std::vector<ExternalProperty<VH, Complex>>
integrated_fields_complex(
        const std::vector<ExternalProperty<HEH, double>>& _gradients,
        const TriMesh& _mesh,
        const VH _seed_vh = VH(0),
        const Complex _seed_value = Complex(1.0));

/**
 * Differentiate a complex scalar field.
 */
ExternalProperty<HEH, double>
differentiate_field(
        const ExternalProperty<VH, Complex>& _u,
        const TriMesh& _mesh);

std::vector<ExternalProperty<HEH, double>>
differentiate_fields(
        const std::vector<ExternalProperty<VH, Complex>>& _us,
        const TriMesh& _mesh);

/**
 * Integrate a 1-form along loop.
 */
double
integrate_loop(
        const ExternalProperty<HEH, double>& _gradient,
        const PrimalLoop& _loop);

MatXd
integrate_loops(
        const std::vector<ExternalProperty<HEH, double>>& _fields,
        const PrimalLoops& _loops);

ExternalProperty<HEH, double>
combine_fields(
        const TriMesh& _mesh,
        const VecXi& _coeffs,
        const std::vector<ExternalProperty<HEH, double>>& _fields);

ExternalProperty<HEH, double>
combine_fields(
        const TriMesh& _mesh,
        const VecXd& _coeffs,
        const std::vector<ExternalProperty<HEH, double>>& _fields);

std::vector<ExternalProperty<HEH, double>>
transform_fields(
        const TriMesh& _mesh,
        const MatXi& _M,
        const std::vector<ExternalProperty<HEH, double>>& _fields);

}
