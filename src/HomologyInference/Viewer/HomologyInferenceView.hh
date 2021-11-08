/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/HomologyInference.hh>
#include <HomologyInference/Viewer/GlowDraw.hh>

namespace HomologyInference
{

struct HomologyInferenceView
{
    explicit HomologyInferenceView(const HomologyInferenceResult& _matching);

    void reset_M_preview();
    void update_M_preview();

    void update_mesh_A();
    void update_mesh_B();

    void update_isolines_A();
    void update_isolines_B();
    void guess_best_isoline_offsets();

    void update_vtpm_lines();
    void update_vtpm_points();

    void render_scene_A();
    void render_scene_B();

    void view_interactive();
    void view_A();
    void view_B();

    // Static data
    const HomologyInferenceResult* matching;

    // State
    MatXi M_preview;
    std::vector<ExternalProperty<HEH, double>> gf_omegaA_A;
    std::vector<ExternalProperty<HEH, double>> gf_omegaB_M_B;

    std::vector<ExternalProperty<HEH, double>> gf_reconstructed;
    std::vector<ExternalProperty<VH, Complex>> igf_reconstructed;

    // Viewer settings
    enum class FieldADisplay : int
    {
        None,
        VF,
        SF,
    };
    constexpr static const char* FieldADisplayLabels[] =
    {
        "None",
        "Basis Cocycle of A",
        "Basis Cocycle of A as Periodic Potential",
    };

    enum class FieldBDisplay : int
    {
        None,
        VF,
        MappedSF,
        MappedVF,
        ReconstructedVF,
    };
    constexpr static const char* FieldBDisplayLabels[] =
    {
        "None",
        "Basis Cocycle of B",
        "Basis Cocycle from A as Periodic Potential on B",
        "Basis Cocycle from A on B",
        "Reconstruction of Cocycle from A using Cocycles of B",
    };

    FieldADisplay field_A_display = FieldADisplay::None;
    FieldBDisplay field_B_display = FieldBDisplay::None;
    int field_A_index = 0;
    int field_B_index = 0;

    enum class VectorFieldDisplay : int
    {
        OnFaces,
        OnEdges,
    };
    constexpr static const char* VectorFieldDisplayLabels[] =
    {
        "OnFaces",
        "OnEdges",
    };
    VectorFieldDisplay vector_field_display = VectorFieldDisplay::OnFaces;
    float vector_field_scale = 0.005;
    float vector_field_line_width = 0.001;

    bool show_isolines = false;
    bool rotate_isolines_by_omega = true;
    float isoline_width = 0.006;
    std::array<std::unique_ptr<bool[]>, 2> isoline_visible; // Nasty: We can't use std::vector<bool> here because ImGui needs bool* for element access.
    std::array<std::vector<float>, 2> isoline_offset;
    std::vector<Color> generator_colors;

    std::vector<Color> color_palette;

    bool show_vtpm_lines = false;
    DrawStyle vtpm_line_style = WidthWorld(0.001);

    bool show_vtpm_points = false;
    DrawStyle vtpm_point_style = WidthWorld(0.004);

    // Cached Renderables
    gv::SharedRenderable r_mesh_A;
    gv::SharedRenderable r_mesh_B;
    GlowDraw d_mesh_A;
    GlowDraw d_mesh_B;

    GlowDraw d_isolines_A;
    GlowDraw d_isolines_B;
    GlowDraw d_vtpm_lines;
    GlowDraw d_vtpm_points_A;
    GlowDraw d_vtpm_points_B;
};

}
