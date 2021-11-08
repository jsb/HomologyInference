/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/Types.hh>
#include <HomologyInference/Viewer/GlowDraw.hh>
#include <HomologyInference/Viewer/RWTHColors.hh>

namespace HomologyInference
{

const DrawStyle homology_landmark_style = WidthWorld(0.01);
const DrawStyle homology_path_style = WidthWorld(0.005, false);

void view_integrable_gradient_field_isolines(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _offset = 0.0,
        const Color _color = RWTH_BLACK,
        const DrawStyle& _style = homology_path_style);

/// Finds an _offset parameter for view_integrable_gradient_field_isolines
/// that yields the shortest isoline.
double guess_best_isoline_offset(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const int _n_samples = 16);

double guess_best_isoline_offset_for_homology_class(
        const TriMesh& _mesh,
        const std::vector<ExternalProperty<HEH, double>>& _gfs,
        const MatXi& _omega,
        const VecXi& _h,
        const int _n_samples = 16);

void view_homology_class_as_isoline(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const std::vector<ExternalProperty<HEH, double>>& _gfs,
        const MatXi& _omega,
        const VecXi& _h,
        const double _offset,
        const Color& _c,
        const DrawStyle& _style = homology_path_style);

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, double>& _field,
        const Color& _color_from,
        const Color& _color_to);

void view_complex_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Complex>& _field);

void view_integrable_gradient_field_on_edges(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _scale = 0.005,
        const double _line_width = 0.001,
        const Color _color = RWTH_BLACK);

void view_integrable_gradient_field_interpolated(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _scale = 0.005,
        const double _line_width = 0.001,
        const Color _color = RWTH_BLACK);

void view_integrable_gradient_field_interpolated(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _scale = 0.005,
        const double _line_width = 0.001,
        const Color _color = RWTH_BLACK);

void view_path(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const PrimalPath& _path,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& _style = homology_path_style);

void view_path(
        const TriMesh& _mesh,
        const PrimalPath& _path,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& _style = homology_path_style);

void view_path(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const DualPath& _path,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& _style = homology_path_style);

void view_path(
        const TriMesh& _mesh,
        const DualPath& _path,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& _style = homology_path_style);

void view_path_with_arrows(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const PrimalPath& _path,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& _style = homology_path_style);

void view_path_with_arrows(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const DualPath& _path,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& _style = homology_path_style);

}
