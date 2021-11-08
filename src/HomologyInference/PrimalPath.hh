/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

struct DualPath;

struct PrimalPath
{
    std::vector<HEH> hehs;

    /// Removes successive edges that go back and forth in opposite directions.
    /// If the path is closed, _move_loop_basepoint enables removing redundant start / end elements,
    /// effectively moving the basepoint.
    template<typename Mesh>
    void shrink(const Mesh& _mesh, const bool _move_loop_basepoint = false);

    template<typename Mesh>
    bool is_closed(const Mesh& _mesh) const;

    template<typename Mesh>
    void assert_valid(const Mesh& _mesh) const;

    template<typename Mesh>
    SVH vh_start(const Mesh& _mesh) const;

    template<typename Mesh>
    SVH vh_end(const Mesh& _mesh) const;

    template<typename Mesh>
    double embedded_length(const Mesh& _mesh) const;

    template<typename Mesh>
    PrimalPath reversed(const Mesh& _mesh) const;
};

template<typename Mesh>
PrimalPath concatenate(
        const Mesh& _mesh,
        const PrimalPath& _p0,
        const PrimalPath& _p1);

template<typename Mesh>
PrimalPath concatenate(
        const Mesh& _mesh,
        const std::vector<PrimalPath>& _segments);

DualPath
to_dual_path(
        const TriMesh& _mesh,
        const PrimalPath& _path);

}
