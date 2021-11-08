/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <HomologyInference/DualPath.hh>
#include <HomologyInference/PrimalPath.hh>
#include <HomologyInference/Types.hh>
#include <HomologyInference/Utils/Filesystem.hh>
#include <HomologyInference/Viewer/RWTHColors.hh>
#include <HomologyInference/Viewer/IDraw.hh>
#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/view.hh>

namespace HomologyInference
{

/// Make TriMesh renderable. Allows for gv::view(mesh).
gv::SharedGeometricRenderable make_renderable(
        const TriMesh& _mesh);

gv::detail::raii_config default_style();

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const glow::viewer::camera_transform& _cam_pos,
        const tg::ivec2& _size = tg::ivec2(1920, 1080),
        const bool _transparent = true);

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const tg::ivec2& _size = tg::ivec2(1920, 1080),
        const bool _transparent = true);

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const glow::SharedTexture2D& _texture);

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const TexCoords& _uvs,
        const glow::SharedTexture2D& _texture);

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors);

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Color>& _colors);

TexCoords projected_tex_coords(
        const TriMesh& _mesh,
        const Vec3d& _u_3d = Vec3d(1.0, 0.0, 0.0),
        const Vec3d& _v_3d = Vec3d(0.0, 1.0, 0.0));

void projected_tex_coords_to_property(
        TriMesh& _mesh,
        const Vec3d& _u_3d = Vec3d(1.0, 0.0, 0.0),
        const Vec3d& _v_3d = Vec3d(0.0, 1.0, 0.0));

gv::detail::raii_view_closer view_mesh(
        const gv::SharedRenderable& _r,
        const std::string _caption = "");

gv::detail::raii_view_closer view_mesh(
        const TriMesh& _mesh,
        const Color& _color = RWTH_WHITE,
        const std::string _caption = "");

void view_vertex_colors(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Color>& _colors);

void view_face_colors(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors);

void view_halfedge_colors(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, Color>& _colors);

void view_wireframe(
        const TriMesh& _mesh,
        const Color& _color = RWTH_BLUE,
        const DrawStyle& = default_line_style);

void view_caption(
        const std::string& _s);

}
