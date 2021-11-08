#include <HomologyInference/Genus.hh>
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

    // Load meshes and input map
    TriMesh mesh_A = read_mesh(DATA_PATH / "pretzel.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "three_hole_smooth.obj");
    normalize_surface_area(mesh_A);
    normalize_surface_area(mesh_B);
    VertexToPointMap vtpm = read_vertex_to_point_map(DATA_PATH / "pretzel_on_three_hole_smooth.vtpm", mesh_A, mesh_B);
    align_rigid(mesh_A, mesh_B, vtpm);
    center_mesh(mesh_A);
    center_mesh(mesh_B);

    // Compute homology bases
    PrimalLoops hom_basis_A = homology_basis(mesh_A);
    PrimalLoops hom_basis_B = homology_basis(mesh_B);

    auto cfg_style = default_style();
    for (const auto& obj : {HomologyInferenceSettings::IntegralMatching, HomologyInferenceSettings::FieldAlignment})
    {
        // Infer homology map
        HomologyInferenceSettings settings;
        settings.objective = obj;
        HomologyInferenceResult result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, vtpm, settings);

        // Visualization
        HomologyInferenceView view(result);
        view.show_isolines = true;
        view.isoline_width = 0.010f;
        view.isoline_offset[0][0] = 0.943f;
        view.isoline_offset[0][1] = 0.257f;
        view.isoline_offset[0][2] = 0.763f;
        view.isoline_offset[0][3] = 0.039f;
        view.isoline_offset[0][4] = 0.261f;
        view.isoline_offset[0][5] = 0.865f;
        if (obj == HomologyInferenceSettings::IntegralMatching)
        {
            result.title = "Inference by integrating mapped cocycles (Eq. (17))";
            view.isoline_offset[1][0] = 0.250f;
            view.isoline_offset[1][1] = 0.575f;
            view.isoline_offset[1][2] = 0.214f;
            view.isoline_offset[1][3] = 0.290f;
            view.isoline_offset[1][4] = 0.509f;
            view.isoline_offset[1][5] = 0.348f;
        }
        if (obj == HomologyInferenceSettings::FieldAlignment)
        {
            result.title = "Inference via least-squares field alignment (Eq. (18))";
            view.isoline_offset[1][0] = 0.500f;
            view.isoline_offset[1][1] = 0.595f;
            view.isoline_offset[1][2] = 0.275f;
            view.isoline_offset[1][3] = 0.225f;
            view.isoline_offset[1][4] = 0.500f;
            view.isoline_offset[1][5] = 0.775f;
        }
        view.show_vtpm_points = true;
        view.vtpm_point_style = WidthWorld(0.003f);
        view.view_interactive();
    }
}
