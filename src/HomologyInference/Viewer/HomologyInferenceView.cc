/*
 * Author: Janis Born
 */
#include "HomologyInferenceView.hh"

#include <HomologyInference/CohomologyBasis.hh>
#include <HomologyInference/IntegerMatrix.hh>
#include <HomologyInference/Viewer/ComplexColors.hh>
#include <HomologyInference/Viewer/MeshView.hh>
#include <HomologyInference/Viewer/RWTHColorGenerator.hh>
#include <HomologyInference/Viewer/View.hh>
#include <HomologyInference/Viewer/XYZColors.hh>

#include <imgui/imgui.h>

namespace HomologyInference
{

HomologyInferenceView::HomologyInferenceView(const HomologyInferenceResult& _matching) :
    matching(&_matching)
{
    M_preview = matching->M_homology;
    const int n = M_preview.rows();

    for (int mesh_id = 0; mesh_id < 2; ++mesh_id)
    {
        isoline_visible[mesh_id] = std::make_unique<bool[]>(n);
        for (int i = 0; i < n; ++i)
            isoline_visible[mesh_id][i] = true;
        isoline_offset[mesh_id].assign(n, 0.0f);
    }

    update_M_preview();
    guess_best_isoline_offsets();

    generator_colors = RWTHColorGenerator().generate_next_colors(n);

    color_palette.push_back(RWTH_BLUE);
    color_palette.push_back(RWTH_MAGENTA);
    color_palette.push_back(RWTH_YELLOW);
    color_palette.push_back(RWTH_PETROL);
    color_palette.push_back(RWTH_TEAL);
    color_palette.push_back(RWTH_GREEN);
    color_palette.push_back(RWTH_MAY_GREEN);
    color_palette.push_back(RWTH_ORANGE);
    color_palette.push_back(RWTH_RED);
    color_palette.push_back(RWTH_BORDEAUX);
    color_palette.push_back(RWTH_PURPLE);
    color_palette.push_back(RWTH_LILAC);
}

void HomologyInferenceView::update_M_preview()
{
    const auto& mesh_A = matching->mesh_A;
    const auto& mesh_B = matching->mesh_B;
    const auto& gf_A = matching->gf_A;
    const auto& gf_B = matching->gf_B;
    const auto& omega_A = matching->omega_A;
    const auto& omega_B = matching->omega_B;

    const int n = gf_A.size();
    gf_omegaA_A.clear();
    gf_omegaB_M_B.clear();

    for (int i = 0; i < n; ++i)
    {
        VecXi hom_A = VecXi::Zero(n);
        hom_A[i] = 1;
        VecXi hom_B = M_preview * hom_A;

        VecXi ip_B = hom_B;
        VecXi ip_A = hom_A;
        if (rotate_isolines_by_omega)
            ip_B = omega_B * hom_B;
        if (rotate_isolines_by_omega)
            ip_A = omega_A * hom_A;

        auto field_A = combine_fields(mesh_A, ip_A, gf_A);
        auto field_B = combine_fields(mesh_B, ip_B, gf_B);
        gf_omegaA_A.push_back(field_A);
        gf_omegaB_M_B.push_back(field_B);
    }

    gf_reconstructed.clear();
    igf_reconstructed.clear();
    for (int i = 0; i < n; ++i)
    {
        VecXi cohom_A = VecXi::Zero(n);
        cohom_A[i] = 1;
        VecXi cohom_B = matching->M_cohomology * cohom_A;

        auto rec_gf_B = combine_fields(mesh_B, cohom_B, gf_B);
        auto rec_igf_B = integrate_field_complex(rec_gf_B, mesh_B);
        gf_reconstructed.push_back(rec_gf_B);
        igf_reconstructed.push_back(rec_igf_B);
    }
}

void HomologyInferenceView::update_mesh_A()
{
    const auto& mesh_A = matching->mesh_A;
    const ExternalProperty<HEH, double>* vector_field = nullptr;
    const ExternalProperty<VH, Complex>* scalar_field = nullptr;

    switch (field_A_display)
    {
        case FieldADisplay::None:
        {
            break;
        }
        case FieldADisplay::VF:
        {
            vector_field = &matching->gf_A[field_A_index];
            break;
        }
        case FieldADisplay::SF:
        {
            scalar_field = &matching->igf_A[field_A_index];
            break;
        }
    }

    // Fail-safe if an uninitialized field was selected.
    if (scalar_field && scalar_field->empty())
    {
        scalar_field = nullptr;
        field_A_display = FieldADisplay::None;
    }
    if (vector_field && vector_field->empty())
    {
        vector_field = nullptr;
        field_A_display = FieldADisplay::None;
    }

    // Mesh Renderable
    if (scalar_field)
    {
        ExternalProperty<VH, Color> colors(mesh_A);
        for (const auto& vh : mesh_A.vertices())
            colors[vh] = complex_color((*scalar_field)[vh]);
        r_mesh_A = make_renderable(mesh_A, colors);
    }
    else
    {
        r_mesh_A = make_renderable(mesh_A);
    }

    // Draw
    d_mesh_A.clear();
    if (vector_field)
    {
        if (vector_field_display == VectorFieldDisplay::OnEdges)
        {
            view_integrable_gradient_field_on_edges(
                        d_mesh_A,
                        mesh_A,
                        *vector_field,
                        vector_field_scale,
                        vector_field_line_width);
        }
        else if (vector_field_display == VectorFieldDisplay::OnFaces)
        {
            view_integrable_gradient_field_interpolated(
                        d_mesh_A,
                        mesh_A,
                        *vector_field,
                        vector_field_scale,
                        vector_field_line_width);
        }
    }
}

void HomologyInferenceView::update_mesh_B()
{
    const auto& mesh_B = matching->mesh_B;
    const ExternalProperty<HEH, double>* vector_field = nullptr;
    const ExternalProperty<VH, Complex>* scalar_field = nullptr;

    switch (field_B_display)
    {
        case FieldBDisplay::None:
        {
            break;
        }
        case FieldBDisplay::VF:
        {
            vector_field = &matching->gf_B[field_B_index];
            break;
        }
        case FieldBDisplay::MappedSF:
        {
            scalar_field = &matching->igf_A_on_B[field_B_index];
            break;
        }
        case FieldBDisplay::MappedVF:
        {
            vector_field = &matching->gf_A_on_B[field_B_index];
            break;
        }
        case FieldBDisplay::ReconstructedVF:
        {
            vector_field = &gf_reconstructed[field_B_index];
            break;
        }
    }

    // Fail-safe if an uninitialized field was selected.
    if (scalar_field && scalar_field->empty())
    {
        scalar_field = nullptr;
        field_B_display = FieldBDisplay::None;
    }
    if (vector_field && vector_field->empty())
    {
        vector_field = nullptr;
        field_B_display = FieldBDisplay::None;
    }

    // Mesh Renderable
    if (scalar_field)
    {
        ExternalProperty<VH, Color> colors(mesh_B);
        for (const auto& vh : mesh_B.vertices())
            colors[vh] = complex_color((*scalar_field)[vh]);
        r_mesh_B = make_renderable(mesh_B, colors);
    }
    else
    {
        r_mesh_B = make_renderable(mesh_B);
    }

    // Draw
    d_mesh_B.clear();
    if (vector_field)
    {
        if (vector_field_display == VectorFieldDisplay::OnEdges)
        {
            view_integrable_gradient_field_on_edges(
                        d_mesh_B,
                        mesh_B,
                        *vector_field,
                        vector_field_scale,
                        vector_field_line_width);
        }
        else if (vector_field_display == VectorFieldDisplay::OnFaces)
        {
            view_integrable_gradient_field_interpolated(
                        d_mesh_B,
                        mesh_B,
                        *vector_field,
                        vector_field_scale,
                        vector_field_line_width);
        }
    }
}

void HomologyInferenceView::update_isolines_A()
{
    d_isolines_A.clear();
    const int n = matching->loops_A.size();
    const auto& mesh_A = matching->mesh_A;

    for (int i = 0; i < n; ++i)
        if (isoline_visible[0][i] > 0)
            view_integrable_gradient_field_isolines(
                        d_isolines_A,
                        mesh_A,
                        gf_omegaA_A[i],
                        isoline_offset[0][i],
                        generator_colors[i],
                        WidthWorld(isoline_width));
}

void HomologyInferenceView::update_isolines_B()
{
    d_isolines_B.clear();
    const int n = matching->loops_B.size();
    const auto& mesh_B = matching->mesh_B;

    for (int i = 0; i < n; ++i)
        if (isoline_visible[1][i] > 0)
            view_integrable_gradient_field_isolines(
                        d_isolines_B,
                        mesh_B,
                        gf_omegaB_M_B[i],
                        isoline_offset[1][i],
                        generator_colors[i],
                        WidthWorld(isoline_width));
}

void HomologyInferenceView::guess_best_isoline_offsets()
{
    const int n_generators = matching->loops_A.size();

    for (int mesh_i = 0; mesh_i < 2; ++mesh_i)
    {
        const TriMesh* mesh = nullptr;
        const std::vector<ExternalProperty<HEH, double>>* fields = nullptr;
        if (mesh_i == 0)
        {
            mesh = &matching->mesh_A;
            fields = &gf_omegaA_A;
        }
        if (mesh_i == 1)
        {
            mesh = &matching->mesh_B;
            fields = &gf_omegaB_M_B;
        }

        for (int generator_i = 0; generator_i < n_generators; ++generator_i)
        {
            const ExternalProperty<HEH, double>& field = (*fields)[generator_i];
            isoline_offset[mesh_i][generator_i] = guess_best_isoline_offset(*mesh, field);
        }
    }
}

void HomologyInferenceView::update_vtpm_lines()
{
    const auto& mesh_A = matching->mesh_A;
    const auto& mesh_B = matching->mesh_B;
    const auto& vtpm = matching->vtpm;
    d_vtpm_lines.clear();
    for (const auto& vh_A : mesh_A.vertices())
    {
        const auto& p_A = mesh_A.point(vh_A);
        const auto& bary_B = vtpm[vh_A];
        if (!bary_B.is_valid())
            continue;
        const auto& p_B = bary_B.point(mesh_B);
        d_vtpm_lines.line(p_A, p_B, RWTH_BLACK, vtpm_line_style);
    }
}

void HomologyInferenceView::update_vtpm_points()
{
    const auto& mesh_A = matching->mesh_A;
    const auto& mesh_B = matching->mesh_B;
    const auto& vtpm = matching->vtpm;
    ExternalProperty<VH, Color> vtpm_colors = xyz_colors(mesh_A);

    d_vtpm_points_A.clear();
    d_vtpm_points_B.clear();
    for (const auto& vh_A : mesh_A.vertices())
    {
        Color color = vtpm_colors[vh_A];
        const auto& p_A = mesh_A.point(vh_A);
        const auto& bary_B = vtpm[vh_A];
        if (!bary_B.is_valid())
            continue;
        const auto& p_B = bary_B.point(mesh_B);
        d_vtpm_points_A.point(p_A, color, vtpm_point_style);
        d_vtpm_points_B.point(p_B, color, vtpm_point_style);
    }
}

void HomologyInferenceView::render_scene_A()
{
    auto v = gv::view();

    gv::view(r_mesh_A);
    d_mesh_A.view();

    if (show_isolines)
        d_isolines_A.view();

    if (show_vtpm_lines)
        d_vtpm_lines.view();

    if (show_vtpm_points)
        d_vtpm_points_A.view();
}

void HomologyInferenceView::render_scene_B()
{
    auto v = gv::view();

    gv::view(r_mesh_B);
    d_mesh_B.view();

    if (show_isolines)
        d_isolines_B.view();

    if (show_vtpm_lines)
        d_vtpm_lines.view();

    if (show_vtpm_points)
        d_vtpm_points_B.view();
}

void HomologyInferenceView::view_interactive()
{
    const int n = matching->loops_A.size();

    update_M_preview();
    update_mesh_A();
    update_mesh_B();
    update_isolines_A();
    update_isolines_B();
    update_vtpm_lines();
    update_vtpm_points();

    gv::interactive([&] (double /* dt */)
    {
        bool M_preview_updated = false;

        bool mesh_needs_update[2];
        mesh_needs_update[0] = false;
        mesh_needs_update[1] = false;

        bool field_index_updated[2];
        field_index_updated[0] = false;
        field_index_updated[1] = false;

        bool isolines_need_update[2];
        isolines_need_update[0] = false;
        isolines_need_update[1] = false;

        bool vtpm_lines_need_update = false;
        bool vtpm_points_need_update = false;

        // GUI
        ImGui::Begin("Homology Inference Result");

        if (!matching->title.empty())
            ImGui::Text("%s", matching->title.c_str());

        ImGui::Checkbox("Show Input Map as Points", &show_vtpm_points);
        if (show_vtpm_points)
            vtpm_points_need_update |= ImGui::DragFloat("Input Map Point Size", &vtpm_point_style.width, 0.00005f, 0.0f, 0.2f);

        ImGui::Checkbox("Show Input Map as Lines", &show_vtpm_lines);
        if (show_vtpm_lines)
            vtpm_lines_need_update |= ImGui::DragFloat("Input Map Line Width", &vtpm_line_style.width, 0.00005f, 0.0f, 0.2f);

        if (ImGui::CollapsingHeader("Cocycle Display"))
        {
            mesh_needs_update[0] |= ImGui::Combo("On A", (int*)(&field_A_display), FieldADisplayLabels, std::size(FieldADisplayLabels));
            mesh_needs_update[1] |= ImGui::Combo("On B", (int*)(&field_B_display), FieldBDisplayLabels, std::size(FieldBDisplayLabels));

            if (field_A_display != FieldADisplay::None)
                field_index_updated[0] |= ImGui::InputInt("Cocycle A Index", &field_A_index);

            if (field_B_display != FieldBDisplay::None)
                field_index_updated[1] |= ImGui::InputInt("Cocycle B Index", &field_B_index);

            ImGui::Separator();

            if (ImGui::Combo("Vector Field Display", (int*)(&vector_field_display), VectorFieldDisplayLabels, std::size(VectorFieldDisplayLabels)))
            {
                mesh_needs_update[0] = true;
                mesh_needs_update[1] = true;
            }
        }

        if (ImGui::CollapsingHeader("Homology Map M"))
        {
            ImGui::Columns(M_preview.cols());
            for (int row = 0; row < M_preview.rows(); ++row)
            {
                ImGui::PushID(row);
                for (int col = 0; col < M_preview.cols(); ++col)
                {
                    ImGui::PushID(col);
                    M_preview_updated |= ImGui::InputInt("##M", &M_preview(row, col), 0);
                    ImGui::NextColumn();
                    ImGui::PopID();
                }
                ImGui::PopID();
            }
            ImGui::Columns(1);

            if (is_symplectic_homology_map(M_preview, matching->omega_A, matching->omega_B))
                ImGui::Text("(M is symplectic)");
            else
                ImGui::Text("(M is not symplectic)");

            if (ImGui::Button("Transpose##M"))
            {
                M_preview = M_preview.transpose().eval();
                M_preview_updated = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Invert##M"))
            {
                if (is_unimodular(M_preview))
                {
                    M_preview = integer_unimodular_matrix_inverse(M_preview);
                    M_preview_updated = true;
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Round M_real##M"))
            {
                M_preview = closest_integer_matrix(matching->M_cohomology_unconstrained);
                M_preview_updated = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Clear##M"))
            {
                M_preview.setZero();
                M_preview_updated = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset##M"))
            {
                M_preview = matching->M_homology;
                M_preview_updated = true;
            }
        }

        if (ImGui::CollapsingHeader("Omega_A"))
        {
            const auto& omega_A = matching->omega_A;
            ImGui::Columns(omega_A.cols());
            for (int row = 0; row < omega_A.rows(); ++row)
            {
                for (int col = 0; col < omega_A.cols(); ++col)
                {
                    ImGui::Text("%d", omega_A(row, col));
                    ImGui::NextColumn();
                }
            }
            ImGui::Columns(1);
        }

        if (ImGui::CollapsingHeader("Omega_B"))
        {
            const auto& omega_B = matching->omega_B;
            ImGui::Columns(omega_B.cols());
            for (int row = 0; row < omega_B.rows(); ++row)
            {
                for (int col = 0; col < omega_B.cols(); ++col)
                {
                    ImGui::Text("%d", omega_B(row, col));
                    ImGui::NextColumn();
                }
            }
            ImGui::Columns(1);
        }

        ImGui::Checkbox("Show Representative Cycles", &show_isolines);
        if (show_isolines)
        {
            if (ImGui::DragFloat("Cycle Width", &isoline_width, 0.00005f, 0.0f, 0.2f))
            {
                isolines_need_update[0] = true;
                isolines_need_update[1] = true;
            }

            for (int mesh_id = 0; mesh_id < 2; ++mesh_id)
            {
                const std::string mesh_name = (mesh_id == 0 ? "A" : "B");
                if (ImGui::CollapsingHeader(("Mesh " + mesh_name + " Cycles").c_str()))
                {
                    ImGui::PushID(mesh_id);
                    for (int i = 0; i < n; ++i)
                    {
                        ImGui::PushID(i);
                        const std::string isolines_label = mesh_name + " " + std::to_string(i);

                        isolines_need_update[mesh_id] |= ImGui::Checkbox(isolines_label.c_str(), &isoline_visible[mesh_id][i]);

                        bool color_changed = false;
                        ImGui::SameLine();
                        bool open_popup = ImGui::ColorButton("Color", (ImVec4&)(generator_colors[i][0]));

                        if (ImGui::BeginDragDropTarget())
                        {
                            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(IMGUI_PAYLOAD_TYPE_COLOR_3F))
                            {
                                float* fdata = (float*)payload->Data;
                                generator_colors[i] = Color(fdata[0], fdata[1], fdata[2], 1.0);
                                color_changed = true;
                            }
                            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(IMGUI_PAYLOAD_TYPE_COLOR_4F))
                            {
                                float* fdata = (float*)payload->Data;
                                generator_colors[i] = Color(fdata[0], fdata[1], fdata[2], 1.0);
                                color_changed = true;
                            }
                            ImGui::EndDragDropTarget();
                        }

                        if (open_popup)
                        {
                            ImGui::OpenPopup("CustomPicker");
                        }
                        if (ImGui::BeginPopup("CustomPicker"))
                        {
                            color_changed |= ImGui::ColorPicker4("##Picker", (float*)(&generator_colors[i][0]), ImGuiColorEditFlags_NoSidePreview | ImGuiColorEditFlags_NoSmallPreview);
                            ImGui::Separator();
                            for (int pal_i = 0; pal_i < color_palette.size(); ++pal_i)
                            {
                                ImGui::PushID(pal_i);
                                if (pal_i > 0)
                                    ImGui::SameLine();
                                if (ImGui::ColorButton("##Palette", (ImVec4&)(color_palette[pal_i]), ImGuiColorEditFlags_NoAlpha | ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoTooltip))
                                {
                                    generator_colors[i] = color_palette[pal_i];
                                    color_changed = true;
                                }
                                ImGui::PopID();
                            }
                            ImGui::EndPopup();
                        }

                        isolines_need_update[0] |= color_changed;
                        isolines_need_update[1] |= color_changed;

                        if (isoline_visible[mesh_id][i])
                        {
                            ImGui::SameLine();
                            if (ImGui::DragFloat("", &isoline_offset[mesh_id][i], 0.001f, -FLT_MAX, FLT_MAX))
                            {
                                isoline_offset[mesh_id][i] -= std::floor(isoline_offset[mesh_id][i]);
                                isolines_need_update[mesh_id] = true;
                            }
                        }

                        ImGui::PopID();
                    }
                    ImGui::PopID();
                }
            }
        }

        ImGui::End();

        if (M_preview_updated)
        {
            update_M_preview();
            isolines_need_update[1] = true;
            mesh_needs_update[1] = true;
        }

        if (isolines_need_update[0])
            update_isolines_A();
        if (isolines_need_update[1])
            update_isolines_B();

        if (field_index_updated[0])
        {
            field_A_index = std::clamp(field_A_index, 0, n-1);
            update_vtpm_points();
            mesh_needs_update[0] = true;
        }
        if (field_index_updated[1])
        {
            field_B_index = std::clamp(field_B_index, 0, n-1);
            update_vtpm_points();
            mesh_needs_update[1] = true;
        }

        if (mesh_needs_update[0])
            update_mesh_A();
        if (mesh_needs_update[1])
            update_mesh_B();

        if (vtpm_lines_need_update)
            update_vtpm_lines();

        if (vtpm_points_need_update)
            update_vtpm_points();

        // Rendering
        {
            auto style = default_style();
            auto grid = gv::grid();
            render_scene_A();
            render_scene_B();
        }
    });
}

void HomologyInferenceView::view_A()
{
    update_M_preview();
    update_mesh_A();
    update_isolines_A();
    update_vtpm_lines();
    update_vtpm_points();

    render_scene_A();
}

void HomologyInferenceView::view_B()
{
    update_M_preview();
    update_mesh_B();
    update_isolines_B();
    update_vtpm_lines();
    update_vtpm_points();

    render_scene_B();
}

}
