/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/Types.hh>
#include <HomologyInference/Utils/ExternalProperty.hh>

#include <set>

namespace HomologyInference
{

using BlockedEdgeProperty   = ExternalProperty<EH, bool>;
using BlockedVertexProperty = ExternalProperty<VH, bool>;

template<typename Mesh>
bool is_blocked(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const EH _eh);

template<typename Mesh>
bool is_blocked(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _heh);

template<typename Mesh>
bool is_blocked(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const VH _vh);

template<typename Mesh>
bool is_blocked(
        const Mesh& _mesh,
        const BlockedVertexProperty& _blocked_vhs,
        const VH _vh);

template<typename Mesh>
bool is_blocked(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const BlockedVertexProperty& _blocked_vhs,
        const VH _vh);

template<typename Mesh>
HEH sector_start(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _sector_heh);

template<typename Mesh>
HEH sector_start(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const VH _vh,
        const FH _sector_fh);

template<typename Mesh>
std::set<FH> faces_in_sector(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _sector_heh);

/// Returns cw-most outgoing halfedge in each sector.
/// Returns any halfedge if there are no sectors.
template<typename Mesh>
std::vector<SHEH> sectors(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const VH _vh);

/// Returns cw-most outgoing halfedge in each sector of from-vertex.
/// Returns any halfedge if there are no sectors.
/// Starts enumeration at the provided halfedge.
template<typename Mesh>
std::vector<SHEH> sectors_starting_at(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked_ehs,
        const HEH _heh_start);

}
