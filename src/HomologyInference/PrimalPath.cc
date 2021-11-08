/*
 * Author: Janis Born
 */

#include "PrimalPath.hh"

#include <HomologyInference/DualPath.hh>

#include <optional>

namespace HomologyInference
{

template<typename Mesh>
void PrimalPath::shrink(
        const Mesh& _mesh,
        const bool _move_loop_basepoint)
{
    bool good = false;
    do
    {
        // (Optional) cancellations that affect loop start / endpoints.
        if (_move_loop_basepoint && is_closed(_mesh))
        {
            const auto& heh_start = hehs.front();
            const auto& heh_end_opp = _mesh.opposite_halfedge_handle(hehs.back());
            if (heh_start == heh_end_opp)
            {
                hehs.erase(hehs.begin());
                hehs.pop_back();
                good = false;
                ISM_ASSERT(is_closed(_mesh));
                continue;
            }
        }
        // Interior cancellations
        for (int i0 = 0; i0 < hehs.size() - 1; ++i0)
        {
            const int i1 = i0 + 1;
            const auto& heh0 = hehs[i0];
            const auto& heh1 = hehs[i1];
            const auto& heh0_opp = _mesh.opposite_halfedge_handle(heh0);
            if (heh1 == heh0_opp)
            {
                hehs.erase(hehs.begin() + i0, hehs.begin() + i1 + 1);
                good = false;
                break;
            }
        }
        good = true;
    }
    while (!good);
    assert_valid(_mesh);
}
template void PrimalPath::shrink<TriMesh>(const TriMesh&, const bool);
template void PrimalPath::shrink<PolyMesh>(const PolyMesh&, const bool);

template<typename Mesh>
bool PrimalPath::is_closed(const Mesh& _mesh) const
{
    if (hehs.size() == 0)
        return true;
    else
        return vh_start(_mesh) == vh_end(_mesh);
}
template bool PrimalPath::is_closed<TriMesh>(const TriMesh&) const;
template bool PrimalPath::is_closed<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
void PrimalPath::assert_valid(const Mesh& _mesh) const
{
    if (hehs.size() < 2)
        return;

    for (int i0 = 0; i0 < hehs.size() - 1; ++i0)
    {
        const int i1 = i0 + 1;
        const auto& heh0 = hehs[i0];
        const auto& heh1 = hehs[i1];
        ISM_ASSERT(heh0.is_valid());
        ISM_ASSERT(heh1.is_valid());
        const auto& heh0_to   = _mesh.to_vertex_handle(heh0);
        const auto& heh1_from = _mesh.from_vertex_handle(heh1);
        ISM_ASSERT(heh0_to == heh1_from);
    }
}
template void PrimalPath::assert_valid<TriMesh>(const TriMesh&) const;
template void PrimalPath::assert_valid<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
SVH PrimalPath::vh_start(const Mesh& _mesh) const
{
    ISM_ASSERT(!hehs.empty());
    return OpenMesh::make_smart(hehs.front(), _mesh).from();
}
template SVH PrimalPath::vh_start<TriMesh>(const TriMesh&) const;
template SVH PrimalPath::vh_start<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
SVH PrimalPath::vh_end(const Mesh& _mesh) const
{
    ISM_ASSERT(!hehs.empty());
    return OpenMesh::make_smart(hehs.back(), _mesh).to();
}
template SVH PrimalPath::vh_end<TriMesh>(const TriMesh&) const;
template SVH PrimalPath::vh_end<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
double PrimalPath::embedded_length(const Mesh& _mesh) const
{
    double result = 0.0;
    for (const auto& heh : hehs)
        result += _mesh.calc_edge_length(heh);
    return result;
}
template double PrimalPath::embedded_length<TriMesh>(const TriMesh&) const;
template double PrimalPath::embedded_length<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
PrimalPath PrimalPath::reversed(const Mesh& _mesh) const
{
    PrimalPath result = *this;
    std::reverse(begin(result.hehs), end(result.hehs));
    for (auto& heh : result.hehs)
    {
        heh = _mesh.opposite_halfedge_handle(heh);
    }
    return result;
}
template PrimalPath PrimalPath::reversed<TriMesh>(const TriMesh&) const;
template PrimalPath PrimalPath::reversed<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
PrimalPath concatenate(
        const Mesh& _mesh,
        const PrimalPath& _p0,
        const PrimalPath& _p1)
{
    if (_p0.hehs.empty())
        return _p1;
    if (_p1.hehs.empty())
        return _p0;
    ISM_ASSERT(_p0.vh_end(_mesh) == _p1.vh_start(_mesh));
    PrimalPath result;
    std::copy(cbegin(_p0.hehs), cend(_p0.hehs), back_inserter(result.hehs));
    std::copy(cbegin(_p1.hehs), cend(_p1.hehs), back_inserter(result.hehs));
    result.assert_valid(_mesh);
    return result;
}
template PrimalPath concatenate<TriMesh>(const TriMesh&, const PrimalPath&, const PrimalPath&);
template PrimalPath concatenate<PolyMesh>(const PolyMesh&, const PrimalPath&, const PrimalPath&);

template<typename Mesh>
PrimalPath concatenate(
        const Mesh& _mesh,
        const std::vector<PrimalPath>& _segments)
{
    PrimalPath result;
    for (const auto& s : _segments)
    {
        if (s.hehs.empty())
            continue;

        if (!result.hehs.empty())
            ISM_ASSERT(result.vh_end(_mesh) == s.vh_start(_mesh));

        std::copy(cbegin(s.hehs), cend(s.hehs), back_inserter(result.hehs));
    }
    result.assert_valid(_mesh);
    return result;
}
template PrimalPath concatenate<TriMesh>(const TriMesh&, const std::vector<PrimalPath>&);
template PrimalPath concatenate<PolyMesh>(const PolyMesh&, const std::vector<PrimalPath>&);

DualPath
to_dual_path(
        const TriMesh& _mesh,
        const PrimalPath& _path)
{
    _path.assert_valid(_mesh);
    const bool primal_closed = _path.is_closed(_mesh);
    const bool primal_starts_on_boundary = _mesh.is_boundary(_path.vh_start(_mesh));
    const bool primal_ends_on_boundary = _mesh.is_boundary(_path.vh_end(_mesh));

    if (!primal_closed)
    {
        // Non-closed paths must start and end at the boundary (for now).
        ISM_ASSERT(primal_starts_on_boundary);
        ISM_ASSERT(primal_ends_on_boundary);
    }

    std::vector<HEH> corridor;
    if (!primal_closed)
        corridor.push_back(HEH(-1));
    for (const auto& heh : _path.hehs)
        corridor.push_back(heh);
    if (!primal_closed)
        corridor.push_back(HEH(-1));

    const auto& left_fan = [&](SHEH heh0, SHEH heh1) -> std::optional<DualPath>
    {
        bool starts_on_boundary = false;
        bool ends_on_boundary = false;
        if (!heh0.is_valid())
        {
            starts_on_boundary = true;
            ISM_ASSERT(heh1.from().is_boundary());
            heh0 = heh1.from().halfedge().opp();
        }
        if (!heh1.is_valid())
        {
            ends_on_boundary = true;
            ISM_ASSERT(heh0.to().is_boundary());
            heh1 = heh0.to().halfedge().prev().opp();
        }

        ISM_ASSERT_EQ(heh0.to(), heh1.from());
        DualPath result;
        if (starts_on_boundary)
            result.hehs.push_back(heh0.opp());

        SHEH heh_current = heh0.opp();
        while (true)
        {
            if (heh_current == heh1)
                break;
            if (heh_current.opp().is_boundary())
                return std::nullopt;

            heh_current = heh_current.opp().next();

            if (heh_current == heh1)
                break;

            result.hehs.push_back(heh_current);
        }
        if (ends_on_boundary)
            result.hehs.push_back(heh1);

        if (!starts_on_boundary)
            result.fh_start = heh0.face();
        if (!ends_on_boundary)
            result.fh_end = heh1.face();

        return result;
    };

    const auto& right_fan = [&](SHEH heh0, SHEH heh1) -> std::optional<DualPath>
    {
        bool starts_on_boundary = false;
        bool ends_on_boundary = false;
        if (!heh0.is_valid())
        {
            starts_on_boundary = true;
            ISM_ASSERT(heh1.from().is_boundary());
            heh0 = heh1.from().halfedge().prev();
        }
        if (!heh1.is_valid())
        {
            ends_on_boundary = true;
            ISM_ASSERT(heh0.to().is_boundary());
            heh1 = heh0.to().halfedge();
        }

        ISM_ASSERT_EQ(heh0.to(), heh1.from());
        DualPath result;
        if (starts_on_boundary)
            result.hehs.push_back(heh0);

        SHEH heh_current = heh0.opp();
        while (true)
        {
            if (heh_current == heh1)
                break;
            if (heh_current.is_boundary())
                return std::nullopt;

            heh_current = heh_current.prev().opp();

            if (heh_current == heh1)
                break;

            result.hehs.push_back(heh_current.opp());
        }
        if (ends_on_boundary)
            result.hehs.push_back(heh1.opp());

        if (!starts_on_boundary)
            result.fh_start = heh0.opp().face();
        if (!ends_on_boundary)
            result.fh_end = heh1.opp().face();

        return result;
    };

    DualPath result;
    {
        FH prev_fh;
        const int last_i = primal_closed ? corridor.size() : corridor.size() - 1;
        for (int i = 0; i < last_i; ++i)
        {
            const SHEH& heh0 = OpenMesh::make_smart(corridor[i], _mesh);
            const SHEH& heh1 = OpenMesh::make_smart(corridor[(i + 1) % corridor.size()], _mesh);

            if (heh0.is_valid() && heh1.is_valid())
                ISM_ASSERT_EQ(heh0.to(), heh1.from());

            std::optional<DualPath> lf = left_fan(heh0, heh1);
            std::optional<DualPath> rf = right_fan(heh0, heh1);
            ISM_ASSERT(lf.has_value() || rf.has_value());

            std::optional<DualPath> best_fan;
            double best_fan_length = INF_DOUBLE;
            for (const auto& fan : {lf, rf})
            {
                if (fan)
                {
                    const double fan_length = fan->embedded_length(_mesh);
                    if (fan_length < best_fan_length)
                    {
                        best_fan_length = fan_length;
                        best_fan = fan;
                    }
                }
            }
            ISM_ASSERT(best_fan.has_value());

            if (prev_fh.is_valid())
            {
                if (best_fan->fh_start != prev_fh)
                {
                    if (heh0.face() == prev_fh)
                        result.hehs.push_back(heh0);
                    else if (heh0.opp().face() == prev_fh)
                        result.hehs.push_back(heh0.opp());
                    else
                        ISM_ERROR_throw("");
                }
            }
            std::copy(best_fan->hehs.begin(), best_fan->hehs.end(), std::back_inserter(result.hehs));

            prev_fh = best_fan->fh_end;
        }
    }

    result.fh_start = _mesh.face_handle(result.hehs.front());
    result.fh_end = _mesh.opposite_face_handle(result.hehs.back());

    const bool dual_starts_on_boundary = _mesh.is_boundary(result.hehs.front());
    const bool dual_ends_on_boundary = _mesh.is_boundary(_mesh.opposite_halfedge_handle(result.hehs.back()));
    if (primal_closed)
    {
        ISM_ASSERT(result.is_closed(_mesh));
    }
    else
    {
        ISM_ASSERT_EQ(primal_starts_on_boundary, dual_starts_on_boundary);
        ISM_ASSERT_EQ(primal_ends_on_boundary, dual_ends_on_boundary);
    }
    result.assert_valid(_mesh);
    return result;
}

}
