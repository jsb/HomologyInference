#include <HomologyInference/Genus.hh>
#include <HomologyInference/HomologyBasis.hh>
#include <HomologyInference/HomologyInference.hh>
#include <HomologyInference/Types.hh>
#include <HomologyInference/Utils/IO.hh>
#include <HomologyInference/Utils/MeshNormalization.hh>
#include <HomologyInference/VertexToPointMap.hh>
#include <HomologyInference/Viewer/HomologyInferenceView.hh>
#include <HomologyInference/Viewer/MeshView.hh>
#include <HomologyInference/Viewer/View.hh>

#include <OpenMesh/Core/Utils/PropertyManager.hh>

int main()
{
    using namespace HomologyInference;
    glow::glfw::GlfwContext ctx;

    // Load meshes and input map
    TriMesh mesh_A = read_mesh(DATA_PATH / "dog_on_box.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "elephant_on_box.obj");
    normalize_mesh(mesh_A);
    normalize_mesh(mesh_B);
    VertexToPointMap vtpm = read_vertex_to_point_map(DATA_PATH / "dog_on_elephant.vtpm", mesh_A, mesh_B);

    // Compute homology bases
    PrimalLoops hom_basis_A = homology_basis(mesh_A, VH(5086));
    PrimalLoops hom_basis_B = homology_basis(mesh_B, VH(5164));

    // Infer homology map
    HomologyInferenceResult result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, vtpm);

    // Display
    auto cfg_style = default_style();
    HomologyInferenceView view(result);
    view.show_isolines = true;
    view.isoline_offset[0][3] = 0.4f;
    view.isoline_offset[1][3] = 0.4f;
    view.isoline_width = 0.011f;
    view.show_vtpm_points = true;
    view.vtpm_point_style = WidthWorld(0.002f);
    view.view_interactive();
}
