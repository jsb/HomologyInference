/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/Utils/BarycentricPoint.hh>
#include <HomologyInference/Utils/ExternalProperty.hh>

namespace HomologyInference
{

using VertexToPointMap = ExternalProperty<VH, BarycentricPoint>;

/// A VTPM is dense if every vertex of A is mapped to a point on B.
/// In a sparse VTPM, there are some vertices that map to invalid BarycentricPoints.
bool is_dense(const VertexToPointMap& _vtpm);

/// Generates a VTPM from A to B by projecting each vertex of A to the closest point on the surface of B.
VertexToPointMap generate_vtpm_by_projection(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B);

/// Generates a VTPM from A to B from a set of given landmark correspondences.
/// Warning: In general, this is a sparse VTPM.
VertexToPointMap generate_vtpm_from_landmarks(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const std::vector<VH>& _landmarks_A,
        const std::vector<VH>& _landmarks_B);

VertexToPointMap generate_vtpm_from_overlay_meshes(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const TriMesh& _overlay_on_A,
        const TriMesh& _overlay_on_B);

/// Maps the complex scalar field _field_A (given at vertices of A) onto the surface of B using _vtpm.
/// Interpolates and resamples the mapped field onto vertices of B.
/// _w_interp controls the weight of the interpolation constraints (0 --> only harmonic, inf --> only interpolation)
ExternalProperty<VH, Complex> interpolate_field(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        const ExternalProperty<VH, Complex>& _field_A,
        const double _w_interp);

std::vector<ExternalProperty<VH, Complex>> interpolate_fields(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        const std::vector<ExternalProperty<VH, Complex>>& _fields_A,
        const double _w_interp);

/// Sanitizes the VTPM by shrinking all barycentric coordinates a bit towards the face interiors.
void sanitize_vtpm(
        VertexToPointMap& _vtpm,
        const double _epsilon = 1e-6);

void randomly_permute_vtpm(
        VertexToPointMap& _vtpm,
        int _n_transpositions);

void randomly_permute_vtpm(
        VertexToPointMap& _vtpm,
        double _ratio);

VertexToPointMap randomly_permuted_vtpm(
        const VertexToPointMap& _vtpm,
        int _n_transpositions);

VertexToPointMap randomly_permuted_vtpm(
        const VertexToPointMap& _vtpm, double _ratio);

void subsample_vtpm(
        VertexToPointMap& _vtpm,
        double _ratio);

void subsample_vtpm(
        VertexToPointMap& _vtpm,
        int _n_keep);

VertexToPointMap subsampled_vtpm(
        const VertexToPointMap& _vtpm,
        double _ratio);

VertexToPointMap subsampled_vtpm(
        const VertexToPointMap& _vtpm,
        int _n_keep);

}
