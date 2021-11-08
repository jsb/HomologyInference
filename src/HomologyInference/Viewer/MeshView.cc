/*
 * Author: Patrick Schmidt
 */

#include "MeshView.hh"
#include <HomologyInference/Viewer/ComplexColors.hh>
#include <HomologyInference/Viewer/GlowDraw.hh>
#include <HomologyInference/Utils/Helpers.hh>
#include <HomologyInference/Utils/IO.hh>

namespace HomologyInference
{

gv::SharedGeometricRenderable make_renderable(
        const TriMesh& _mesh)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);

    return gv::make_renderable(pos);
}

gv::detail::raii_config default_style()
{
    return gv::config(gv::no_grid, gv::no_outline, gv::maybe_empty, gv::background_color(tg::color3::white), gv::ssao_power(0.5f));
}

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const glow::viewer::camera_transform& _cam_pos,
        const tg::ivec2& _size,
        const bool _transparent)
{
    make_file_directory(_file_path);

    return gv::config(
                _cam_pos,
                gv::headless_screenshot(
                               _size,
                               64,
                               _file_path.string(),
                               _transparent ? GL_RGBA8 : GL_RGB8));
}

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const tg::ivec2& _size,
        const bool _transparent)
{
    return gv::config(
                gv::headless_screenshot(
                               _size,
                               64,
                               _file_path.string(),
                               _transparent ? GL_RGBA8 : GL_RGB8));
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const glow::SharedTexture2D& _texture)
{
    ISM_ASSERT(_mesh.has_halfedge_texcoords2D());

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set texture coordinates
    auto uvs = m.halfedges().map([&] (auto h)
    {
        const VH vh_from(h.vertex_from().idx.value);
        const VH vh_to(h.vertex_to().idx.value);
        const HEH heh = _mesh.find_halfedge(vh_from, vh_to);
        return tg::pos2(_mesh.texcoord2D(heh)[0], _mesh.texcoord2D(heh)[1]);
    });

    // Flip texture due to different conventions
    configure(*r, gv::textured(uvs, _texture).flip());

    return r;
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const TexCoords& _uvs,
        const glow::SharedTexture2D& _texture)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set texture coordinates
    auto uvs = m.halfedges().map([&] (auto h)
    {
        const VH vh_from(h.vertex_from().idx.value);
        const VH vh_to(h.vertex_to().idx.value);
        const HEH heh = _mesh.find_halfedge(vh_from, vh_to);
        return tg::pos2(_uvs[heh][0], _uvs[heh][1]);
    });

    // Flip texture due to different conventions
    configure(*r, gv::textured(uvs, _texture).flip());

    return r;
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set face colors
    auto tg_color = [] (auto c) { return tg::color(c[0], c[1], c[2], c[3]); };
    auto colors = m.faces().map([&] (auto f) { return tg_color(_colors[FH(f.idx.value)]); });
    configure(*r, colors);

    return r;
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set face colors
    auto tg_color = [] (auto c) { return tg::color(c[0], c[1], c[2], c[3]); };
    auto colors = m.vertices().map([&] (auto v) { return tg_color(_colors[VH(v.idx.value)]); });
    configure(*r, colors);

    return r;
}

TexCoords projected_tex_coords(
        const TriMesh& _mesh,
        const Vec3d& _u_3d,
        const Vec3d& _v_3d)
{
    return TexCoords(_mesh, _mesh.halfedges().to_vector([&] (auto h)
    {
        const auto p = _mesh.point(h.to());
        return Vec2d(_u_3d.dot(p), _v_3d.dot(p));
    }));
}

void projected_tex_coords_to_property(
        TriMesh& _mesh,
        const Vec3d& _u_3d,
        const Vec3d& _v_3d)
{
    auto coords = projected_tex_coords(_mesh, _u_3d, _v_3d);

    _mesh.request_halfedge_texcoords2D();
    for (const auto h : _mesh.halfedges())
        _mesh.set_texcoord2D(h, coords[h]);
}

gv::detail::raii_view_closer view_mesh(
        const gv::SharedRenderable& _r,
        const std::string _caption)
{
    return gv::view(_r, _caption);
}

gv::detail::raii_view_closer view_mesh(
        const TriMesh& _mesh,
        const Color& _color,
        const std::string _caption)
{
    ExternalProperty<FH, Color> colors(_mesh, _color);
    return gv::view(make_renderable(_mesh, colors), _caption);
}

void view_vertex_colors(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set vertex colors
    auto colors = m.vertices().make_attribute<tg::color3>();
    for (auto v : m.vertices())
        colors[v] = tg::color3(_colors[VH(v.idx.value)]);

    gv::view(pos, colors);
}

void view_face_colors(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set face colors
    auto colors = m.faces().make_attribute<tg::color3>();
    for (auto f : m.faces())
    {
        const auto v0 = f.halfedges().first().vertex_from();
        const auto v1 = f.halfedges().first().vertex_to();
        const FH fh = _mesh.find_halfedge(VH(v0.idx.value), VH(v1.idx.value)).face();
        colors[f] = tg::color3(_colors[fh]);
    }

    gv::view(pos, colors);
}

void view_halfedge_colors(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set corner colors
    auto colors = m.halfedges().make_attribute<tg::color3>();
    for (auto he : m.halfedges())
    {
        const auto v0 = he.vertex_from();
        const auto v1 = he.vertex_to();
        const HEH heh = _mesh.find_halfedge(VH(v0.idx.value), VH(v1.idx.value));
        colors[he] = tg::color3(_colors[heh]);
    }

    gv::view(pos, colors);
}

void view_wireframe(
        const TriMesh& _mesh,
        const Color& _color,
        const DrawStyle& _style)
{
    auto v = gv::view();

    GlowDraw draw;
    for (auto e : _mesh.edges())
        draw.line(e, _mesh, _color, _style);

    draw.view();
}

void view_caption(const std::string& _s)
{
    gv::view(gv::make_renderable(std::vector<tg::pos3>()), gv::maybe_empty, _s);
}

}
