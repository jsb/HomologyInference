/*
 * Author: Janis Born
 */
#include "View.hh"

#include <HomologyInference/CohomologyBasis.hh>
#include <HomologyInference/Utils/IO.hh>
#include <HomologyInference/Viewer/ComplexColors.hh>
#include <HomologyInference/Viewer/GlowDraw.hh>
#include <HomologyInference/Viewer/MeshView.hh>
#include <HomologyInference/Viewer/RWTHColorGenerator.hh>

namespace HomologyInference
{

void view_integrable_gradient_field_isolines(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _offset,
        const Color _color,
        const DrawStyle& _style)
{
    const auto integrated = integrate_field_real(_field, _mesh);

    for (const auto& fh : _mesh.faces())
    {
        const auto& vh0 = fh.halfedge().from();

        std::map<double, Vec3d> opposite;
        double current_f = integrated[vh0] + _offset;
        for (const auto& heh : fh.halfedges())
        {
            const double next_f = current_f + _field[heh];

            Vec3d p0 = _mesh.point(heh.from());
            Vec3d p1 = _mesh.point(heh.to());
            double min_f = current_f;
            double max_f = next_f;
            if (max_f < min_f)
            {
                std::swap(min_f, max_f);
                std::swap(p0, p1);
            }

            const int min_i = std::ceil(min_f);
            const int max_i = std::floor(max_f);
            for (int i = min_i; i <= max_i; ++i)
            {
                const double sample_f = i;
                const double t = (sample_f - min_f) / (max_f - min_f);
                const Vec3d p = (1 - t) * p0 + t * p1;

                const auto it = opposite.find(sample_f);
                if (it == opposite.end())
                {
                    opposite[sample_f] = p;
                }
                else
                {
                    const Vec3d& p_opp = it->second;
                    if (p.allFinite() && p_opp.allFinite())
                        _draw.line(p, p_opp, _color, _style);
                }
            }
            current_f = next_f;
        }
    }
}

double guess_best_isoline_offset(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const int _n_samples)
{
    double best_offset = 0.0;
    double best_total_length = INF_DOUBLE;
    for (int sample_i = 0; sample_i < _n_samples; ++sample_i)
    {
        const double offset = static_cast<double>(sample_i) / _n_samples;
        GlowDraw draw;
        view_integrable_gradient_field_isolines(draw, _mesh, _field, offset);
        double total_length = 0.0;
        for (const auto& e : draw.m_lines.edges())
            total_length += pm::edge_length(e, draw.pos_lines);

        if (total_length < best_total_length)
        {
            best_total_length = total_length;
            best_offset = offset;
        }
    }
    return best_offset;
}

double guess_best_isoline_offset_for_homology_class(
        const TriMesh& _mesh,
        const std::vector<ExternalProperty<HEH, double>>& _gfs,
        const MatXi& _omega,
        const VecXi& _h,
        const int _n_samples)
{
    const VecXi h_ortho = _omega * _h;
    auto gf = combine_fields(_mesh, h_ortho, _gfs);
    return guess_best_isoline_offset(_mesh, gf, _n_samples);
}

void view_homology_class_as_isoline(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const std::vector<ExternalProperty<HEH, double>>& _gfs,
        const MatXi& _omega,
        const VecXi& _h,
        const double _offset,
        const Color& _c,
        const DrawStyle& _style)
{
    ISM_ASSERT_EQ(_gfs.size(), _h.size());
    ISM_ASSERT_EQ(_omega.rows(), _h.size());
    ISM_ASSERT_EQ(_omega.cols(), _h.size());

    const VecXi h_ortho = _omega * _h;
    auto gf = combine_fields(_mesh, h_ortho, _gfs);
    view_integrable_gradient_field_isolines(_draw, _mesh, gf, _offset, _c, _style);
}

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, double>& _field,
        const Color& _color_from,
        const Color& _color_to)
{
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);

    auto v_colors = m.vertices().make_attribute<tg::color3>();
    const double min = _field.as_eigen().minCoeff();
    const double max = _field.as_eigen().maxCoeff();

    tg::color3 tg_color_from(_color_from[0], _color_from[1], _color_from[2]);
    tg::color3 tg_color_to(_color_to[0], _color_to[1], _color_to[2]);

    for (auto v : m.vertices())
    {
        const double val = _field[VH(v.idx.value)];
        const double lambda = (val - min) / (max - min);
        v_colors[v] = tg::mix(tg_color_from, tg_color_to, lambda);
    }

    gv::view(pos, v_colors);
}

void view_complex_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Complex>& _field)
{
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);

    auto v_colors = m.vertices().make_attribute<tg::color3>();
    for (auto v : m.vertices())
    {
        const VH vh = _mesh.vertex_handle(int(v));
        v_colors[v] = complex_color_tg(_field[vh]);
    }

    gv::view(pos, v_colors);
}

void view_integrable_gradient_field_on_edges(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _scale,
        const double _line_width,
        const Color _color)
{
    for (const auto& eh : _mesh.edges())
    {
        const auto& heh = eh.h(0);
        const auto magnitude = _field[heh];

        const double arrow_length = 2.0 * magnitude / _mesh.calc_edge_length(heh);
        const Vec3d dir = arrow_length * _mesh.calc_edge_vector(heh).normalized();
        const Vec3d& p0 = _mesh.calc_edge_midpoint(eh);
        const Vec3d p1 = p0 + _scale * dir;
        _draw.line(p0, p1, RWTH_BLACK, WidthWorld(_line_width));
        _draw.point(p0, _color, WidthWorld(1.25 * _line_width));
    }
}

void view_integrable_gradient_field_interpolated(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _scale,
        const double _line_width,
        const Color _color)
{
    for (const auto& fh : _mesh.faces())
    {
        Vec3d dir(0, 0, 0);
        const auto f_normal = _mesh.calc_normal(fh);
        const auto f_area = _mesh.calc_face_area(fh);
        for (const auto& heh_in : _mesh.fh_range(fh))
        {
            const auto& heh_out = _mesh.next_halfedge_handle(heh_in);
            const auto& heh_opp = _mesh.next_halfedge_handle(heh_out);
            const auto magnitude = _field[heh_in] - _field[heh_out];

            const auto local_dir_opp = _mesh.calc_edge_vector(heh_opp);
            const auto local_dir = f_normal.cross(local_dir_opp);

            dir += local_dir * magnitude;
        }
        dir /= 2.0 * f_area;

        const auto p0 = _mesh.calc_face_centroid(fh);
        const auto p1 = p0 + _scale * dir;
        _draw.line(p0, p1, RWTH_BLACK, WidthWorld(_line_width));
        _draw.point(p0, _color, WidthWorld(1.25 * _line_width));
    }
}

void view_integrable_gradient_field_interpolated(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, double>& _field,
        const double _scale,
        const double _line_width,
        const Color _color)
{
    GlowDraw draw;
    view_integrable_gradient_field_interpolated(
                draw,
                _mesh,
                _field,
                _scale,
                _line_width,
                _color);
    draw.view();
}

void view_path(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const PrimalPath& _path,
        const Color& _color,
        const DrawStyle& _style)
{
    for (auto heh : _path.hehs)
        _draw.line(heh, _mesh, _color, _style);
}

void view_path(
        const TriMesh& _mesh,
        const PrimalPath& _path,
        const Color& _color,
        const DrawStyle& _style)
{
    GlowDraw draw;
    view_path(draw, _mesh, _path, _color, _style);
    draw.view();
}

void view_path(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const DualPath& _path,
        const Color& _color,
        const DrawStyle& _style)
{
    for (auto heh : _path.hehs)
    {
        const FH& fh0 = _mesh.face_handle(heh);
        const FH& fh1 = _mesh.opposite_face_handle(heh);
        const EH& eh = _mesh.edge_handle(heh);

        const auto& p0 = _mesh.calc_face_centroid(fh0);
        const auto& p1 = _mesh.calc_edge_midpoint(eh);
        const auto& p2 = _mesh.calc_face_centroid(fh1);

        _draw.line(p0, p1, _color, _style);
        _draw.line(p1, p2, _color, _style);
    }
}

void view_path(
        const TriMesh& _mesh,
        const DualPath& _path,
        const Color& _color,
        const DrawStyle& _style)
{
    GlowDraw draw;
    view_path(draw, _mesh, _path, _color, _style);
    draw.view();
}

void view_path_with_arrows(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const PrimalPath& _path,
        const Color& _color,
        const DrawStyle& _style)
{
    for (auto heh : _path.hehs)
    {
        const Vec3d arrow_dir = _mesh.calc_edge_vector(heh).normalized();
        const FH fh_l = _mesh.face_handle(heh);
        const FH fh_r = _mesh.opposite_face_handle(heh);
        const Vec3d n_l = _mesh.calc_face_normal(fh_l);
        const Vec3d n_r = _mesh.calc_face_normal(fh_r);
        const Vec3d t_l = -arrow_dir.cross(n_l).normalized();
        const Vec3d t_r = arrow_dir.cross(n_r).normalized();

        const double head_width = 0.2 * _mesh.calc_edge_vector(heh).norm();
        const double head_height = 0.2 * _mesh.calc_edge_vector(heh).norm();
        Vec3d p_tip = _mesh.calc_edge_midpoint(heh);
        Vec3d p_l = p_tip + head_width * t_l - head_height * arrow_dir;
        Vec3d p_r = p_tip + head_width * t_r - head_height * arrow_dir;

        _draw.line(heh, _mesh, _color, _style);
        _draw.line(p_tip, p_l, _color, _style);
        _draw.line(p_tip, p_r, _color, _style);
    }
}

void view_path_with_arrows(
        GlowDraw& _draw,
        const TriMesh& _mesh,
        const DualPath& _path,
        const Color& _color,
        const DrawStyle& _style)
{
    for (auto heh : _path.hehs)
    {
        const FH& fh0 = _mesh.face_handle(heh);
        const FH& fh1 = _mesh.opposite_face_handle(heh);
        const EH& eh = _mesh.edge_handle(heh);

        const Vec3d p1 = _mesh.calc_edge_midpoint(eh);
        Vec3d p0 = p1;
        Vec3d p2 = p1;
        if (fh0.is_valid())
        {
            p0 = _mesh.calc_face_centroid(fh0);
            _draw.line(p0, p1, _color, _style);
        }
        if (fh1.is_valid())
        {
            p2 = _mesh.calc_face_centroid(fh1);
            _draw.line(p1, p2, _color, _style);
        }

        if (fh0.is_valid() && fh1.is_valid())
        {
            const Vec3d arrow_dir = (p1 - p0).normalized();
            const Vec3d n = _mesh.calc_face_normal(fh0);
            const Vec3d t = arrow_dir.cross(n).normalized();

            const double head_width = 0.2 * (p2 - p1).norm();
            const double head_height = 0.2 * (p2 - p1).norm();
            Vec3d p_tip = p1;
            Vec3d p_l = p_tip - head_width * t - head_height * arrow_dir;
            Vec3d p_r = p_tip + head_width * t - head_height * arrow_dir;

            _draw.line(p_tip, p_l, _color, _style);
            _draw.line(p_tip, p_r, _color, _style);
        }
    }
}

}
