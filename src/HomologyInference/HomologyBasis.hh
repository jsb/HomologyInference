/*
 * Author: Janis Born
 */
#pragma once

#include <HomologyInference/BlockedProperty.hh>
#include <HomologyInference/DualPath.hh>
#include <HomologyInference/PrimalPath.hh>
#include <HomologyInference/Types.hh>

namespace HomologyInference
{

// Primal: Edges going from vertex to vertex.
// Dual:   Edges going from face to face.
// Chain:  Unordered collection of oriented edges.
// Path:   Ordered, continuous chain of edges.
// Loop:   Closed path.

// PrimalLoop: A continuous chain of primal edges, going from vertex to vertex.
using PrimalLoop = PrimalPath;
using PrimalLoops = std::vector<PrimalLoop>;
using PrimalPaths = std::vector<PrimalPath>;

// DualLoop: A continuous chain of dual edges, going from face to face.
using DualLoop = DualPath;
using DualLoops = std::vector<DualLoop>;
using DualPaths = std::vector<DualPath>;

using PrimalCycle = ExternalProperty<HEH, int>;
using PrimalCycles = std::vector<PrimalCycle>;

using DualCycle = ExternalProperty<HEH, int>;
using DualCycles = std::vector<DualCycle>;

struct PrimalSpanningTree
{
    std::vector<VH> roots;
    ExternalProperty<VH, HEH> parent;
    ExternalProperty<VH, double> distance_to_root;
    ExternalProperty<EH, bool> included;
};

struct DualSpanningTree
{
    ExternalProperty<EH, bool> crossed;
};

template<typename Mesh>
PrimalSpanningTree primal_spanning_tree(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked = BlockedEdgeProperty());

template<typename Mesh>
PrimalSpanningTree primal_spanning_tree(
        const Mesh& _mesh,
        const VH _root,
        const BlockedEdgeProperty& _blocked = BlockedEdgeProperty());

template<typename Mesh>
PrimalSpanningTree primal_spanning_tree(
        const Mesh& _mesh,
        const std::vector<VH>& _roots,
        BlockedEdgeProperty _blocked = BlockedEdgeProperty());

/// Forms a cut graph that includes a given set of vertices.
ExternalProperty<EH, bool> seeded_cut_graph(
        const TriMesh& _mesh,
        const std::vector<VH>& _seeds);

/// Forms a cut graph from all primal edges not crossed by a given DST.
ExternalProperty<EH, bool> primal_cut_graph(
        const TriMesh& _mesh,
        const DualSpanningTree& _dst);

int cut_graph_valence(
        const TriMesh& _mesh,
        ExternalProperty<EH, bool>& _cut_graph,
        const VH _vh);

/// Removes all pendant edges and vertices from _cut_graph.
void simplify_cut_graph(
        const TriMesh& _mesh,
        ExternalProperty<EH, bool>& _cut_graph);

/// Returns a modified version of _cut_graph without pendant edges and vertices.
ExternalProperty<EH, bool> simplified_cut_graph(
        const TriMesh& _mesh,
        const ExternalProperty<EH, bool>& _cut_graph);

template<typename Mesh>
double bridge_cost(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const HEH _bridge);

template<typename Mesh>
double bridge_cost(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const EH _bridge);

/// Returns a loop in the spanning tree obtained by inserting _bridge.
/// _bridge will be the first element of this loop.
template<typename Mesh>
PrimalLoop bridge_loop(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const HEH _bridge);

/// Returns a loop in the spanning tree obtained by inserting _bridge.
template<typename Mesh>
PrimalLoop bridge_loop(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const EH _bridge);

/// Returns a path between two roots of the spanning tree obtained by inserting _bridge.
PrimalPath bridge_path(
        const TriMesh& _mesh,
        const PrimalSpanningTree& _pst,
        const HEH _bridge);

/// Returns a path between two roots of the spanning tree obtained by inserting _bridge.
PrimalPath bridge_path(
        const TriMesh& _mesh,
        const PrimalSpanningTree& _pst,
        const EH _bridge);

/// Returns a max-cost spanning tree on the faces of _mesh.
/// _cost can be Inf for certain edges, implying those edges are blocked and will not be included in the DST.
template<typename Mesh>
DualSpanningTree dual_spanning_tree(
        const Mesh& _mesh,
        const ExternalProperty<EH, double>& _cost);

/// Returns a DST that does not intersect any edges of the given PST.
template<typename Mesh>
DualSpanningTree dual_spanning_tree(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst);

/// Returns a DST that does not intersect any edges of the given PST or of the given additional list of blocked edges.
template<typename Mesh>
DualSpanningTree dual_spanning_tree(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const ExternalProperty<EH, bool>& _blocked_eh);

template<typename Mesh>
std::vector<EH> find_bridges(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const DualSpanningTree& _dst);

PrimalLoops homology_basis(
        const TriMesh& _mesh,
        const VH _root = VH(0));

PrimalCycle to_cycle(
        const TriMesh& _mesh,
        const PrimalLoop& _l);

PrimalCycles to_cycles(
        const TriMesh& _mesh,
        const PrimalLoops& _ls);

DualCycle to_cycle(
        const TriMesh& _mesh,
        const DualPath& _l);

DualCycles to_cycles(
        const TriMesh& _mesh,
        const DualPaths& _ls);

bool is_on_cycle(
        const TriMesh& _mesh,
        const PrimalCycle& _cycle,
        const VH _vh);

bool is_on_any_cycle(
        const TriMesh& _mesh,
        const PrimalCycles& _cycles,
        const VH _vh);

/// Number of incident edges that are part of the cycle
int cycle_valence(
        const TriMesh& _mesh,
        const PrimalCycle& _cycle,
        const VH _vh);

double intersection_number_double(
        const TriMesh& _mesh,
        const PrimalLoop& _a,
        const PrimalCycle& _b);

int intersection_number(
        const TriMesh& _mesh,
        const PrimalLoop& _a,
        const PrimalCycle& _b);

int intersection_number(
        const TriMesh& _mesh,
        const PrimalLoop& _a,
        const PrimalLoop& _b);

MatXi intersection_form(
        const TriMesh& _mesh,
        const PrimalLoops& _loops);

}
