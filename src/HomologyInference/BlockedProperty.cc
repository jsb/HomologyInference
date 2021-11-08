/*
 * Author: Janis Born
 */
#include "BlockedProperty.hh"

namespace HomologyInference
{

template<typename Mesh>
bool is_blocked(const Mesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const EH _eh)
{
    return _blocked_ehs[_eh];
}
template bool is_blocked<TriMesh>(const TriMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const EH _eh);
template bool is_blocked<PolyMesh>(const PolyMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const EH _eh);

template<typename Mesh>
bool is_blocked(const Mesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const HEH _heh)
{
    const auto& eh = _mesh.edge_handle(_heh);
    return is_blocked(_mesh, _blocked_ehs, eh);
}
template bool is_blocked<TriMesh>(const TriMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const HEH _heh);
template bool is_blocked<PolyMesh>(const PolyMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const HEH _heh);

template<typename Mesh>
bool is_blocked(const Mesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const VH _vh)
{
    for (const auto& eh : _mesh.ve_range(_vh))
        if (is_blocked(_mesh, _blocked_ehs, eh))
            return true;
    return false;
}
template bool is_blocked<TriMesh>(const TriMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const VH _vh);
template bool is_blocked<PolyMesh>(const PolyMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const VH _vh);

template<typename Mesh>
bool is_blocked(const Mesh& _mesh, const BlockedVertexProperty& _blocked_vhs, const VH _vh)
{
    return _blocked_vhs[_vh];
}

template<typename Mesh>
bool is_blocked(const Mesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const BlockedVertexProperty& _blocked_vhs, const VH _vh)
{
    return is_blocked(_mesh, _blocked_vhs, _vh) || is_blocked(_mesh, _blocked_ehs, _vh);
}
template bool is_blocked<TriMesh>(const TriMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const BlockedVertexProperty& _blocked_vhs, const VH _vh);
template bool is_blocked<PolyMesh>(const PolyMesh& _mesh, const BlockedEdgeProperty& _blocked_ehs, const BlockedVertexProperty& _blocked_vhs, const VH _vh);

template<typename Mesh>
HEH sector_start(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _heh_in_sector)
{
    ISM_ASSERT(_mesh.is_valid_handle(_heh_in_sector));
    const VH& vh = _mesh.from_vertex_handle(_heh_in_sector);

    // No sectors at the current vertex? Just return the halfedge itself.
    if (!is_blocked(_mesh, _blocked_ehs, vh))
        return _heh_in_sector;

    // Rotate clockwise until a sector boundary is encountered.
    HEH current_heh = _heh_in_sector;
    while (!is_blocked(_mesh, _blocked_ehs, current_heh))
        current_heh = _mesh.cw_rotated_halfedge_handle(current_heh);
    return current_heh;
}
template HEH sector_start<TriMesh>(const TriMesh&, const BlockedEdgeProperty&, const HEH);
template HEH sector_start<PolyMesh>(const PolyMesh&, const BlockedEdgeProperty&, const HEH);

template<typename Mesh>
HEH sector_start(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const VH _vh,
        const FH _sector_fh)
{
    // Find outgoing halfedge indident to fh
    SHEH h_out;
    for (auto h : _mesh.fh_range(_sector_fh))
    {
        if (h.from() == _vh)
        {
            h_out = h;
            break;
        }
    }

    // Rotate h_out cw until a sector boundary is encountered
    return sector_start(_mesh, _blocked_ehs, h_out);
}
template HEH sector_start<TriMesh>(const TriMesh&, const BlockedEdgeProperty&, const VH, const FH);
template HEH sector_start<PolyMesh>(const PolyMesh&, const BlockedEdgeProperty&, const VH, const FH);

template<typename Mesh>
std::set<FH> faces_in_sector(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _sector_heh)
{
    ISM_ASSERT(_mesh.is_valid_handle(_sector_heh));
    const VH& vh = _mesh.from_vertex_handle(_sector_heh);
    std::set<FH> result;

    // If there are no sectors at the current vertex ...
    if (!is_blocked(_mesh, _blocked_ehs, vh))
    {
        // ... return all surrounding faces.
        for (const FH& fh : _mesh.vf_range(vh))
            result.insert(fh);
        return result;
    }
    // Otherwise, we expect _sector_heh to be at a sector boundary.
    ISM_ASSERT(is_blocked(_mesh, _blocked_ehs, _sector_heh));
    HEH current_heh = _sector_heh;
    do
    {
        result.insert(_mesh.face_handle(current_heh));
        current_heh = _mesh.ccw_rotated_halfedge_handle(current_heh);
    }
    while (!is_blocked(_mesh, _blocked_ehs, current_heh));
    return result;
}
template std::set<FH> faces_in_sector<TriMesh>(const TriMesh&, const BlockedEdgeProperty&, const HEH);
template std::set<FH> faces_in_sector<PolyMesh>(const PolyMesh&, const BlockedEdgeProperty&, const HEH);

template<typename Mesh>
std::vector<SHEH> sectors(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const VH _vh)
{
    return sectors_starting_at(_mesh, _blocked_ehs, *_mesh.cvoh_begin(_vh));
}
template std::vector<SHEH> sectors<TriMesh>(const TriMesh&, const BlockedEdgeProperty&, const VH);
template std::vector<SHEH> sectors<PolyMesh>(const PolyMesh&, const BlockedEdgeProperty&, const VH);

template<typename Mesh>
std::vector<SHEH> sectors_starting_at(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _heh_start)
{
    ISM_ASSERT(_mesh.is_valid_handle(_heh_start));
    auto h = OpenMesh::make_smart(_heh_start, _mesh);

    // If there are no sectors return provided halfedge
    if (!is_blocked(_mesh, _blocked_ehs, h.from()))
        return { h };

    std::vector<SHEH> result;
    do
    {
        if (is_blocked(_mesh, _blocked_ehs, h))
            result.push_back(h);

        h = h.opp().next();
    }
    while (h != _heh_start);

    return result;
}
template std::vector<SHEH> sectors_starting_at<TriMesh>(const TriMesh&, const BlockedEdgeProperty&, const HEH);
template std::vector<SHEH> sectors_starting_at<PolyMesh>(const PolyMesh&, const BlockedEdgeProperty&, const HEH);

}
