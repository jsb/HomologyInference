/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

struct DualPath
{
    DualPath() = default;
    DualPath(
            const FH _fh_start,
            const FH _fh_end);

    int length() const;

    bool is_closed(
            const TriMesh& _mesh) const;

    SFH face(
            const int _idx,
            const TriMesh& _mesh) const;

    double embedded_length(
            const TriMesh& _mesh) const;

    DualPath reversed(
            const TriMesh& _mesh) const;

    void assert_start_end_valid(
            const TriMesh& _mesh) const;

    void assert_valid(
            const TriMesh& _mesh) const;

    FH fh_start = FH(-1);
    FH fh_end = FH(-1);
    std::vector<HEH> hehs; // Halfedges crossed, incident to previous face.
};

DualPath concatenate(
        const std::vector<DualPath>& _paths);

bool self_intersecting(
        const DualPath& _path,
        const TriMesh& _mesh);

/// Return path prefix before first self-intersection (exclusive)
DualPath until_self_intresection(
        const DualPath& _path,
        const TriMesh& _mesh);

}
