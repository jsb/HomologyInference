#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/HomologyInference.hh>
#include <HomologyInference/Utils/IO.hh>
#include <HomologyInference/Utils/MeshNormalization.hh>
#include <HomologyInference/VertexToPointMap.hh>
#include <HomologyInference/Viewer/HomologyInferenceView.hh>
#include <HomologyInference/Viewer/MeshView.hh>
#include <HomologyInference/Viewer/View.hh>

#include <glow-extras/glfw/GlfwContext.hh>

using namespace HomologyInference;

int main()
{
    glow::glfw::GlfwContext ctx;

    // Load input meshes
    TriMesh mesh_A = read_mesh(DATA_PATH / "pyramid_frame.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "eiffel_tower.obj");
    normalize_mesh(mesh_A);
    normalize_mesh(mesh_B);

    // Load input map (a sparse set of landmark annotations)
    std::vector<VH> landmarks_A = read_landmarks(DATA_PATH / "pyramid_frame.pinned");
    std::vector<VH> landmarks_B = read_landmarks(DATA_PATH / "eiffel_tower.pinned");
    VertexToPointMap vtpm = generate_vtpm_from_landmarks(mesh_A, mesh_B, landmarks_A, landmarks_B);

    // Compute homology bases
    PrimalLoops hom_basis_A = homology_basis(mesh_A);
    PrimalLoops hom_basis_B = homology_basis(mesh_B);

    // Infer homology map
    HomologyInferenceResult result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, vtpm);

    // Display
    HomologyInferenceView view(result);
    view.show_isolines = true;
    view.isoline_width = 0.010f;
    view.show_vtpm_points = true;
    view.vtpm_point_style = WidthWorld(0.008f);
    view.view_interactive();
}
