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

    // Read input meshes
    TriMesh mesh_A = read_mesh(DATA_PATH / "mac_blocky.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "mac_smooth.obj");
    normalize_mesh(mesh_A);
    normalize_mesh(mesh_B);

    // Read input maps
    std::vector<VertexToPointMap> vtpms;
    vtpms.push_back(read_vertex_to_point_map(DATA_PATH / "mac_blocky_on_smooth_noisy.vtpm",   mesh_A, mesh_B));
    vtpms.push_back(read_vertex_to_point_map(DATA_PATH / "mac_blocky_on_smooth_twisted.vtpm", mesh_A, mesh_B));

    // Compute homology bases
    PrimalLoops hom_basis_A = homology_basis(mesh_A, VH(7806));
    PrimalLoops hom_basis_B = homology_basis(mesh_B, VH(35));

    for (const auto& vtpm : vtpms)
    {
        // Infer homology map
        auto result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, vtpm);

        // Display
        HomologyInferenceView view(result);
        view.show_isolines = true;
        view.isoline_width = 0.010f;
        view.show_vtpm_points = true;
        view.vtpm_point_style = WidthWorld(0.002f);
        view.view_interactive();
    }
}
