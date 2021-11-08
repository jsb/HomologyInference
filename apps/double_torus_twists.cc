#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/HomologyInference.hh>
#include <HomologyInference/Utils/IO.hh>
#include <HomologyInference/Utils/MeshNormalization.hh>
#include <HomologyInference/VertexToPointMap.hh>
#include <HomologyInference/Viewer/HomologyInferenceView.hh>
#include <HomologyInference/Viewer/MeshView.hh>
#include <HomologyInference/Viewer/View.hh>

#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <glow-extras/glfw/GlfwContext.hh>

using namespace HomologyInference;

int main()
{
    glow::glfw::GlfwContext ctx;

    TriMesh mesh_A = read_mesh(DATA_PATH / "double_torus_round_iso.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "double_torus_blocky.obj");
    normalize_mesh(mesh_A);
    normalize_mesh(mesh_B);

    const Vec3d cog_A = mesh_A.calc_centroid(OpenMesh::MeshHandle());
    const Vec3d cog_B = mesh_B.calc_centroid(OpenMesh::MeshHandle());
    const VH vh_A_base = closest_vertex(mesh_A, cog_A);
    const VH vh_B_base = closest_vertex(mesh_B, cog_B);
    PrimalLoops hom_basis_A = homology_basis(mesh_A, vh_A_base);
    PrimalLoops hom_basis_B = homology_basis(mesh_B, vh_B_base);

    std::vector<VertexToPointMap> vtpms;
    vtpms.push_back(read_vertex_to_point_map(DATA_PATH / "double_torus_round_iso_on_blocky_id.vtpm",           mesh_A, mesh_B));
    vtpms.push_back(read_vertex_to_point_map(DATA_PATH / "double_torus_round_iso_on_blocky_handle_twist.vtpm", mesh_A, mesh_B));
    vtpms.push_back(read_vertex_to_point_map(DATA_PATH / "double_torus_round_iso_on_blocky_half_twist.vtpm",   mesh_A, mesh_B));
    vtpms.push_back(read_vertex_to_point_map(DATA_PATH / "double_torus_round_iso_on_blocky_full_twist.vtpm",   mesh_A, mesh_B));

    auto cfg_style = default_style();
    for (const auto& vtpm : vtpms)
    {
        auto result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, vtpm);
        HomologyInferenceView view(result);
        view.show_isolines = true;
        view.isoline_width = 0.010f;
        view.show_vtpm_points = true;
        view.vtpm_point_style = WidthWorld(0.003f);
        view.view_interactive();
    }
}
