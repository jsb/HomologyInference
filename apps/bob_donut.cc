#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/HomologyInference.hh>
#include <HomologyInference/Utils/IO.hh>
#include <HomologyInference/Utils/MeshNormalization.hh>
#include <HomologyInference/VertexToPointMap.hh>
#include <HomologyInference/Viewer/HomologyInferenceView.hh>
#include <HomologyInference/Viewer/MeshView.hh>
#include <HomologyInference/Viewer/View.hh>

#include <glow-extras/glfw/GlfwContext.hh>
#include <GLFW/glfw3.h>

using namespace HomologyInference;

int main()
{
    glow::glfw::GlfwContext ctx;

    // Load input meshes
    TriMesh mesh_A = read_mesh(DATA_PATH / "bob_tri.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "donut.obj");
    normalize_mesh(mesh_A);
    normalize_mesh(mesh_B);

    // Load input map
    VertexToPointMap vtpm = read_vertex_to_point_map(DATA_PATH / "bob_tri_on_donut.vtpm", mesh_A, mesh_B);

    // Compute homology bases
    PrimalLoops loops_A = homology_basis(mesh_A);
    PrimalLoops loops_B = homology_basis(mesh_B);

    // Infer homology map
    HomologyInferenceResult result = infer_homology_map(mesh_A, mesh_B, loops_A, loops_B, vtpm);

    // Display
    auto cfg_style = default_style();
    HomologyInferenceView view(result);
    view.show_isolines = true;
    view.isoline_width = 0.010f;
    view.show_vtpm_points = true;
    view.vtpm_point_style = WidthWorld(0.003f);
    view.view_interactive();
}
