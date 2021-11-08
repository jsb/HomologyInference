/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/Types.hh>
#include <HomologyInference/VertexToPointMap.hh>

namespace HomologyInference
{

void center_mesh(
        TriMesh& _mesh,
        const TriMesh::Point& _new_origin = {0,0,0});

void normalize_surface_area(
        TriMesh& _mesh,
        double _new_area = 1.0);

void normalize_mesh(
        TriMesh& _mesh);

/**
 * Compute rigid transformation (+scaling)
 * that aligns points V_B to points V_A
 */
Eigen::Affine3d compute_rigid_alignment(
        const MatXd& V_A, // points of A as rows
        const MatXd& V_B, // points of B as rows
        bool _allow_scaling = false);

/**
 * Compute rigid transformation (+scaling)
 * that aligns mesh_B to mesh_A.
 */
Eigen::Affine3d compute_rigid_alignment(
        const std::vector<VH>& _vertices_A,
        const std::vector<VH>& _vertices_B,
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        bool _allow_scaling = false);

Eigen::Affine3d compute_rigid_alignment(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const VertexToPointMap& _vtpm,
        bool _allow_scaling = false);

void align_rigid(
        const std::vector<VH>& _vertices_A,
        const std::vector<VH>& _vertices_B,
        const TriMesh& _mesh_A, // reference
        TriMesh& _mesh_B, // output
        bool _allow_scaling = false);

void align_rigid(
        const TriMesh& _mesh_A, // reference
        TriMesh& _mesh_B, // output
        const VertexToPointMap& _vtpm,
        bool _allow_scaling = false);

}
