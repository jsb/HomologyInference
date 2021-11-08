#include <HomologyInference/Genus.hh>
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
    // Load input meshes
    TriMesh mesh_A = read_mesh(DATA_PATH / "shrec07_108.off");
    TriMesh mesh_B = read_mesh(DATA_PATH / "shrec07_105.off");

    // Roughly align the two meshes via the provided landmark annotations
    std::vector<VH> landmarks_A = read_landmarks(DATA_PATH / "shrec07_108.vts");
    std::vector<VH> landmarks_B = read_landmarks(DATA_PATH / "shrec07_105.vts");
    align_rigid(landmarks_A, landmarks_B, mesh_A, mesh_B, true);

    // Generate an input map via nearest-neighbor projections
    VertexToPointMap vtpm = generate_vtpm_by_projection(mesh_A, mesh_B);
    sanitize_vtpm(vtpm);

    // Visualize the nearest-neighbor projections
    {
        auto cfg_style = default_style();
        auto v = gv::view();
        view_mesh(mesh_A, Color(1.0,1.0,1.0,0.5));
        view_mesh(mesh_B, Color(1.0,1.0,1.0,0.5));
        GlowDraw draw;
        for (const auto& vh_A : mesh_A.vertices())
        {
            const Vec3d p_A = mesh_A.point(vh_A);
            const Vec3d p_B = vtpm[vh_A].point(mesh_B);
            draw.line(p_A, p_B, RWTH_BLACK, WidthWorld(0.002f));
        }
        draw.view();
    }

    // Compute homology bases
    const Vec3d centroid_A = mesh_A.calc_centroid(OpenMesh::MeshHandle());
    const Vec3d centroid_B = mesh_B.calc_centroid(OpenMesh::MeshHandle());
    const VH vh_base_A = closest_vertex(mesh_A, centroid_A);
    const VH vh_base_B = closest_vertex(mesh_B, centroid_B);
    PrimalLoops hom_basis_A = homology_basis(mesh_A, vh_base_A);
    PrimalLoops hom_basis_B = homology_basis(mesh_B, vh_base_B);

    // Infer homology map
    HomologyInferenceResult result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, vtpm);

    // Display
    auto cfg_style = default_style();
    HomologyInferenceView view(result);
    view.show_isolines = true;
    view.isoline_width = 0.020f;
    view.show_vtpm_points = true;
    view.vtpm_point_style = WidthWorld(0.003f);
    view.view_interactive();
}
