#include <HomologyInference/FunctionalMaps.hh>
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

    // Read input meshes
    TriMesh mesh_A = read_mesh(DATA_PATH / "komainu_a_25k.obj");
    TriMesh mesh_B = read_mesh(DATA_PATH / "komainu_b_31k.obj");
    normalize_mesh(mesh_A);
    normalize_mesh(mesh_B);

    // Read functional bases of A and B and functional map between them
    MatXd fbasis_A = read_matrix(DATA_PATH / "komainu_a_25k_eigenbasis.mat");
    MatXd fbasis_B = read_matrix(DATA_PATH / "komainu_b_31k_eigenbasis.mat");
    MatXd C = read_matrix(DATA_PATH / "komainu_fmap.mat");

    // Preview functional map
    {
        MatXd xyz_A(mesh_A.n_vertices(), 3);
        for (const auto& vh : mesh_A.vertices())
            xyz_A.row(vh.idx()) = mesh_A.point(vh);
        VecXd min_xyz = xyz_A.colwise().minCoeff();
        VecXd max_xyz = xyz_A.colwise().maxCoeff();

        MatXd color_A(mesh_A.n_vertices(), 3);
        for (int row = 0; row < color_A.rows(); ++row)
            color_A.row(row) = (xyz_A.row(row) - min_xyz).cwiseQuotient(max_xyz - min_xyz);

        MatXd color_B = map_using_fmap(fbasis_A, fbasis_B, C, color_A);

        ExternalProperty<VH, Color> prop_color_A(mesh_A, RWTH_BLACK);
        ExternalProperty<VH, Color> prop_color_B(mesh_B, RWTH_BLACK);
        for (const auto& vh : mesh_A.vertices())
        {
            const VecXd& c = color_A.row(vh.idx());
            prop_color_A[vh] = Color(c[0], c[1], c[2], 1.0);
        }
        for (const auto& vh : mesh_B.vertices())
        {
            const VecXd& c = color_B.row(vh.idx());
            prop_color_B[vh] = Color(c[0], c[1], c[2], 1.0);
        }
        auto cfg_style = default_style();
        auto cfg_grid = gv::grid();
        view_vertex_colors(mesh_A, prop_color_A);
        view_vertex_colors(mesh_B, prop_color_B);
    }

    // Compute homology basis
    PrimalLoops hom_basis_A = homology_basis(mesh_A);
    PrimalLoops hom_basis_B = homology_basis(mesh_B);

    // Infer homology map
    HomologyInferenceResult result = infer_homology_map(mesh_A, mesh_B, hom_basis_A, hom_basis_B, fbasis_A, fbasis_B, C);

    // Display
    auto cfg_style = default_style();
    HomologyInferenceView view(result);
    view.show_isolines = true;
    view.isoline_width = 0.010f;
    view.isoline_offset[0][0] = 0.438f;
    view.isoline_offset[0][1] = 0.611f;
    view.isoline_offset[0][2] = 0.120f;
    view.isoline_offset[0][3] = 0.176f;
    view.isoline_offset[1][0] = 0.784f;
    view.isoline_offset[1][1] = 0.788f;
    view.isoline_offset[1][2] = 0.755f;
    view.isoline_offset[1][3] = 0.634f;
    view.view_interactive();
}
