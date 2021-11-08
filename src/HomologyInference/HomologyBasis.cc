/*
 * Author: Janis Born
 */

#include "HomologyBasis.hh"

#include <HomologyInference/Genus.hh>
#include <HomologyInference/Utils/Out.hh>
#include <HomologyInference/Utils/Round.hh>
#include <HomologyInference/Utils/Helpers.hh>
#include <HomologyInference/Utils/UnionFind.hh>

#include <queue>

namespace HomologyInference
{

template<typename Mesh>
PrimalSpanningTree
primal_spanning_tree(
        const Mesh& _mesh,
        const BlockedEdgeProperty& _blocked)
{
    return primal_spanning_tree(_mesh, VH(0), _blocked);
}
template PrimalSpanningTree primal_spanning_tree<TriMesh>(const TriMesh&, const BlockedEdgeProperty&);
template PrimalSpanningTree primal_spanning_tree<PolyMesh>(const PolyMesh&, const BlockedEdgeProperty&);

template<typename Mesh>
PrimalSpanningTree
primal_spanning_tree(
        const Mesh& _mesh,
        const VH _root,
        const BlockedEdgeProperty& _blocked)
{
    return primal_spanning_tree(_mesh, std::vector<VH>{_root}, _blocked);
}
template PrimalSpanningTree primal_spanning_tree<TriMesh>(const TriMesh&, const VH, const BlockedEdgeProperty&);
template PrimalSpanningTree primal_spanning_tree<PolyMesh>(const PolyMesh&, const VH, const BlockedEdgeProperty&);

template<typename Mesh>
PrimalSpanningTree
primal_spanning_tree(
        const Mesh& _mesh,
        const std::vector<VH>& _roots,
        BlockedEdgeProperty _blocked)
{
    ISM_ASSERT(_blocked.empty() || _blocked.size_okay(_mesh));
    if (_blocked.empty())
        _blocked.init(_mesh, false);

    PrimalSpanningTree pst;
    pst.roots = _roots;
    pst.parent.init(_mesh, Mesh::InvalidHalfedgeHandle);
    pst.distance_to_root.init(_mesh, INF_DOUBLE);
    pst.included.init(_mesh, false);

    struct Candidate
    {
        VH vh;
        double cost;
        bool operator>(const Candidate _rhs) const { return cost > _rhs.cost; }
    };

    std::priority_queue<Candidate, std::vector<Candidate>, std::greater<>> q;
    for (const auto& root : pst.roots)
    {
        q.push({root, 0.0});
        pst.distance_to_root[root] = 0.0;
        pst.parent[root] = Mesh::InvalidHalfedgeHandle;
    }

    // Build the tree.
    while (!q.empty())
    {
        const Candidate c = q.top();
        q.pop();
        for (const HEH& heh : _mesh.voh_range(c.vh))
        {
            const EH& eh = _mesh.edge_handle(heh);
            if (_blocked[eh])
                continue;

            const HEH& heh_opp = _mesh.opposite_halfedge_handle(heh);
            const VH& to_vh = _mesh.to_vertex_handle(heh);
            const double edge_length = _mesh.calc_edge_length(eh);
            const double new_distance_to_root = pst.distance_to_root[c.vh] + edge_length;
            if (new_distance_to_root < pst.distance_to_root[to_vh])
            {
                pst.distance_to_root[to_vh] = new_distance_to_root;
                pst.parent[to_vh] = heh_opp;
                q.push({to_vh, new_distance_to_root});
            }
        }
    }

    // Mark all edges that are part of the tree.
    for (const VH& vh : _mesh.vertices())
    {
        if (pst.parent[vh].is_valid())
        {
            const HEH heh = pst.parent[vh];
            const EH eh = _mesh.edge_handle(heh);
            pst.included[eh] = true;
        }
    }

    return pst;
}
template PrimalSpanningTree primal_spanning_tree<TriMesh>(const TriMesh&, const std::vector<VH>&, BlockedEdgeProperty);
template PrimalSpanningTree primal_spanning_tree<PolyMesh>(const PolyMesh&, const std::vector<VH>&, BlockedEdgeProperty);

ExternalProperty<EH, bool>
seeded_cut_graph(
        const TriMesh& _mesh,
        const std::vector<VH>& _seeds)
{
    PrimalSpanningTree pst = primal_spanning_tree(_mesh, _seeds);
    DualSpanningTree dst = dual_spanning_tree(_mesh, pst);
    ExternalProperty<EH, bool> result(_mesh, false);

    // Find bridges: Edges that are neither in the PST nor crossed by the DST.
    std::vector<EH> bridges = find_bridges(_mesh, pst, dst);
    for (const auto& bridge_eh : bridges)
    {
        // For each bridge, reconstruct a path connecting two seedpoints.
        PrimalPath path = bridge_path(_mesh, pst, bridge_eh);
        for (const auto& heh : path.hehs)
        {
            const auto& eh = _mesh.edge_handle(heh);
            result[eh] = true;
        }
    }
    return result;
}

ExternalProperty<EH, bool>
primal_cut_graph(
        const TriMesh& _mesh,
        const DualSpanningTree& _dst)
{
    ISM_ASSERT(_dst.crossed.size_okay(_mesh));
    ExternalProperty<EH, bool> result(_mesh);
    for (const auto& eh : _mesh.edges())
    {
        result[eh] = !_dst.crossed[eh];
    }
    return result;
}

int
cut_graph_valence(
        const TriMesh& _mesh,
        ExternalProperty<EH, bool>& _cut_graph,
        const VH _vh)
{
    ISM_ASSERT(_cut_graph.size_okay(_mesh));
    int valence = 0;
    for (const auto& eh : _mesh.ve_range(_vh))
    {
        if (_cut_graph[eh])
        {
            ++valence;
        }
    }
    return valence;
}

void
simplify_cut_graph(
        const TriMesh& _mesh,
        ExternalProperty<EH, bool>& _cut_graph)
{
    ISM_ASSERT(_cut_graph.size_okay(_mesh));
    for (const auto& vh : _mesh.vertices())
    {
        VH current_vh = vh;
        while (cut_graph_valence(_mesh, _cut_graph, current_vh) == 1)
        {
            // Find the incident pendant halfedge, pointing from current_vh to the only adjacent cut graph vertex.
            HEH pendant_halfedge = TriMesh::InvalidHalfedgeHandle;
            for (const auto& heh : _mesh.voh_range(current_vh))
            {
                if (_cut_graph[heh.edge()])
                {
                    pendant_halfedge = heh;
                    break;
                }
            }
            ISM_ASSERT(pendant_halfedge.is_valid());
            const EH pendant_edge = _mesh.edge_handle(pendant_halfedge);

            // Remove the pendant edge from the cut graph and proceed to the adjacent vertex.
            _cut_graph[pendant_edge] = false;
            current_vh = _mesh.to_vertex_handle(pendant_halfedge);
        }
    }
}

ExternalProperty<EH, bool>
simplified_cut_graph(
        const TriMesh& _mesh,
        const ExternalProperty<EH, bool>& _cut_graph)
{
    ExternalProperty<EH, bool> cut_graph_copy = _cut_graph;
    simplify_cut_graph(_mesh, cut_graph_copy);
    return cut_graph_copy;
}

template<typename Mesh>
double
bridge_cost(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const HEH _bridge)
{
    const VH vh0 = _mesh.from_vertex_handle(_bridge);
    const VH vh1 = _mesh.to_vertex_handle(_bridge);
    const double edge_length = _mesh.calc_edge_length(_bridge);
    return edge_length + _pst.distance_to_root[vh0] + _pst.distance_to_root[vh1];
}
template double bridge_cost<TriMesh>(const TriMesh&, const PrimalSpanningTree&, const HEH);
template double bridge_cost<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&, const HEH);

template<typename Mesh>
double
bridge_cost(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const EH _bridge)
{
    const HEH heh = _mesh.halfedge_handle(_bridge, 0);
    return bridge_cost(_mesh, _pst, heh);
}
template double bridge_cost<TriMesh>(const TriMesh&, const PrimalSpanningTree&, const EH);
template double bridge_cost<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&, const EH);

template<typename Mesh>
PrimalLoop
bridge_loop(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const HEH _bridge)
{
    ISM_ASSERT(_pst.roots.size() == 1);
    PrimalLoop loop;

    // The bridge itself.
    loop.hehs.push_back(_bridge);

    // From the "to"-vertex of the bridge to the root.
    {
        auto current_vh = _mesh.to_vertex_handle(_bridge);
        while (_pst.parent[current_vh].is_valid())
        {
            const auto heh = _pst.parent[current_vh];
            const auto next_vh = _mesh.to_vertex_handle(heh);
            loop.hehs.push_back(heh);
            current_vh = next_vh;
        }
    }

    // From the "from"-vertex of the bridge to the root.
    // (We follow this path in reverse order and then reverse the result).
    {
        const size_t segment_start_pos = loop.hehs.size();
        auto current_vh = _mesh.from_vertex_handle(_bridge);
        while (_pst.parent[current_vh].is_valid())
        {
            const auto heh = _pst.parent[current_vh];
            const auto heh_opp = _mesh.opposite_halfedge_handle(heh);
            const auto prev_vh = _mesh.to_vertex_handle(heh);
            loop.hehs.push_back(heh_opp);
            current_vh = prev_vh;
        }
        const size_t segment_end_pos = loop.hehs.size();
        std::reverse(loop.hehs.begin() + segment_start_pos, loop.hehs.begin() + segment_end_pos);
    }

    ISM_ASSERT(!loop.hehs.empty());
    ISM_ASSERT(loop.is_closed(_mesh));
    loop.assert_valid(_mesh);
    return loop;
}
template PrimalLoop bridge_loop<TriMesh>(const TriMesh&, const PrimalSpanningTree&, const HEH);
template PrimalLoop bridge_loop<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&, const HEH);

template<typename Mesh>
PrimalLoop
bridge_loop(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const EH _bridge)
{
    const HEH heh = _mesh.halfedge_handle(_bridge, 0);
    return bridge_loop(_mesh, _pst, heh);
}
template PrimalLoop bridge_loop<TriMesh>(const TriMesh&, const PrimalSpanningTree&, const EH);
template PrimalLoop bridge_loop<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&, const EH);

PrimalPath
bridge_path(
        const TriMesh& _mesh,
        const PrimalSpanningTree& _pst,
        const HEH _bridge)
{
    PrimalLoop loop;

    // From the "from"-vertex of the bridge to the root.
    // (We follow this path in reverse order and then reverse the result).
    {
        auto current_vh = _mesh.from_vertex_handle(_bridge);
        while (_pst.parent[current_vh].is_valid())
        {
            const auto heh = _pst.parent[current_vh];
            const auto heh_opp = _mesh.opposite_halfedge_handle(heh);
            const auto prev_vh = _mesh.to_vertex_handle(heh);
            loop.hehs.push_back(heh_opp);
            current_vh = prev_vh;
        }
        std::reverse(loop.hehs.begin(), loop.hehs.end());
    }

    // The bridge itself.
    loop.hehs.push_back(_bridge);

    // From the "to"-vertex of the bridge to the root.
    {
        auto current_vh = _mesh.to_vertex_handle(_bridge);
        while (_pst.parent[current_vh].is_valid())
        {
            const auto heh = _pst.parent[current_vh];
            const auto next_vh = _mesh.to_vertex_handle(heh);
            loop.hehs.push_back(heh);
            current_vh = next_vh;
        }
    }

    ISM_ASSERT(!loop.hehs.empty());
    loop.assert_valid(_mesh);
    return loop;
}

PrimalPath
bridge_path(
        const TriMesh& _mesh,
        const PrimalSpanningTree& _pst,
        const EH _bridge)
{
    const HEH heh = _mesh.halfedge_handle(_bridge, 0);
    return bridge_path(_mesh, _pst, heh);
}

template<typename Mesh>
DualSpanningTree
dual_spanning_tree(
        const Mesh& _mesh,
        const ExternalProperty<EH, double>& _cost)
{
    ISM_ASSERT(_cost.size_okay(_mesh));

    DualSpanningTree dst;
    dst.crossed.init(_mesh, false);

    struct Candidate
    {
        EH eh;
        double cost;
        bool operator<(const Candidate _rhs) const { return cost < _rhs.cost; }
    };
    std::priority_queue<Candidate> candidates;
    for (const EH& eh : _mesh.edges())
    {
        if (_mesh.is_boundary(eh))
            continue;
        if (!std::isfinite(_cost[eh]))
            continue;
        Candidate c;
        c.eh = eh;
        c.cost = _cost[eh];
        candidates.push(c);
    }

    UnionFind uf(_mesh.n_faces());
    while (!candidates.empty())
    {
        const auto& c = candidates.top();
        const auto& fh_a = _mesh.face_handle(_mesh.halfedge_handle(c.eh, 0));
        const auto& fh_b = _mesh.face_handle(_mesh.halfedge_handle(c.eh, 1));
        if (!uf.equivalent(fh_a.idx(), fh_b.idx())) {
            dst.crossed[c.eh] = true;
            uf.merge(fh_a.idx(), fh_b.idx());
        }
        candidates.pop();
    }

    return dst;
}
template DualSpanningTree dual_spanning_tree<TriMesh>(const TriMesh&, const ExternalProperty<EH, double>&);
template DualSpanningTree dual_spanning_tree<PolyMesh>(const PolyMesh&, const ExternalProperty<EH, double>&);

template<typename Mesh>
DualSpanningTree
dual_spanning_tree(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst)
{
    ExternalProperty<EH, bool> nothing_blocked(_mesh, false);
    return dual_spanning_tree(_mesh, _pst, nothing_blocked);
}
template DualSpanningTree dual_spanning_tree<TriMesh>(const TriMesh&, const PrimalSpanningTree&);
template DualSpanningTree dual_spanning_tree<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&);

template<typename Mesh>
DualSpanningTree
dual_spanning_tree(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const ExternalProperty<EH, bool>& _blocked_eh)
{
    ISM_ASSERT(_blocked_eh.size_okay(_mesh));
    ExternalProperty<EH, double> cost(_mesh);
    for (const auto& eh : _mesh.edges())
    {
        if (_pst.included[eh] || _blocked_eh[eh])
            cost[eh] = INF_DOUBLE;
        else
            cost[eh] = bridge_cost(_mesh, _pst, eh);
    }
    return dual_spanning_tree(_mesh, cost);
}
template DualSpanningTree dual_spanning_tree<TriMesh>(const TriMesh&, const PrimalSpanningTree&, const ExternalProperty<EH, bool>&);
template DualSpanningTree dual_spanning_tree<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&, const ExternalProperty<EH, bool>&);

template<typename Mesh>
std::vector<EH>
find_bridges(
        const Mesh& _mesh,
        const PrimalSpanningTree& _pst,
        const DualSpanningTree& _dst)
{
    std::vector<EH> bridges;
    for (const auto& eh : _mesh.edges())
    {
        if (!_pst.included[eh] && !_dst.crossed[eh])
        {
            bridges.push_back(eh);
        }
    }
    return bridges;
}
template std::vector<EH> find_bridges<TriMesh>(const TriMesh&, const PrimalSpanningTree&, const DualSpanningTree&);
template std::vector<EH> find_bridges<PolyMesh>(const PolyMesh&, const PrimalSpanningTree&, const DualSpanningTree&);

PrimalLoops
homology_basis(
        const TriMesh& _mesh,
        const VH _root)
{
    PrimalLoops result;

    // Interdigitating spanning trees.
    PrimalSpanningTree pst = primal_spanning_tree(_mesh, _root);
    DualSpanningTree dst = dual_spanning_tree(_mesh, pst);

    // Find bridges: Edges that are neither in the PST nor crossed by the DST.
    std::vector<EH> bridges = find_bridges(_mesh, pst, dst);
    ISM_ASSERT_EQ((int)bridges.size(), 2 * genus(_mesh));

    // For each bridge, reconstruct a loop by following paths from the endpoints back to the root.
    for (const auto& bridge_eh : bridges)
    {
        PrimalLoop loop = bridge_loop(_mesh, pst, bridge_eh);
        ISM_ASSERT(!loop.hehs.empty());
        result.push_back(loop);
    }

    return result;
}

PrimalCycle
to_cycle(
        const TriMesh& _mesh,
        const PrimalLoop& _l)
{
    ExternalProperty<HEH, int> result(_mesh);
    for (const auto& heh : _l.hehs)
    {
        result[heh] += 1;
        result[_mesh.opposite_halfedge_handle(heh)] -= 1;
    }
    return result;
}

PrimalCycles
to_cycles(
        const TriMesh& _mesh,
        const PrimalLoops& _ls)
{
    PrimalCycles result;
    for (const auto& l : _ls)
    {
        result.push_back(to_cycle(_mesh, l));
    }
    return result;
}

DualCycle
to_cycle(
        const TriMesh& _mesh,
        const DualPath& _p)
{
    ExternalProperty<HEH, int> result(_mesh);
    for (const auto& heh : _p.hehs)
    {
        result[heh] += 1;
        result[_mesh.opposite_halfedge_handle(heh)] -= 1;
    }
    return result;
}

DualCycles
to_cycles(
        const TriMesh& _mesh,
        const DualPaths& _ps)
{
    DualCycles result;
    for (const auto& l : _ps)
    {
        result.push_back(to_cycle(_mesh, l));
    }
    return result;
}

bool
is_on_cycle(
        const TriMesh& _mesh,
        const PrimalCycle& _cycle,
        const VH _vh)
{
    for (const auto& heh : _mesh.voh_range(_vh))
        if (_cycle[heh] != 0)
            return true;
    return false;
}

bool
is_on_any_cycle(
        const TriMesh& _mesh,
        const PrimalCycles& _cycles,
        const VH _vh)
{
    for (const auto& cycle : _cycles)
        if (is_on_cycle(_mesh, cycle, _vh))
            return true;
    return false;
}

/// Number of incident edges that are part of the cycle
int
cycle_valence(
        const TriMesh& _mesh,
        const PrimalCycle& _cycle,
        const VH _vh)
{
    int valence = 0;

    for (const auto& heh : _mesh.voh_range(_vh))
    {
        if (_cycle[heh] != 0)
            ++valence;
    }

    return valence;
}

double
intersection_number_double(
        const TriMesh& _mesh,
        const PrimalLoop& _a,
        const PrimalCycle& _b)
{
    if (_a.hehs.empty())
        return 0;

    int result2 = 0;

    // Measure intersections at all inner vertices of the path
    int end_index = _a.hehs.size() - 1;
    if (_a.is_closed(_mesh))
        end_index = _a.hehs.size();

    for (int i0 = 0; i0 < end_index; ++i0)
    {
        const int i1 = (i0 + 1) % _a.hehs.size();
        const auto& heh0 = _a.hehs[i0];
        const auto& heh1 = _a.hehs[i1];
        ISM_ASSERT_EQ(_mesh.to_vertex_handle(heh0), _mesh.from_vertex_handle(heh1));

        if (heh1 == _mesh.opposite_halfedge_handle(heh0))
            continue;

        const auto heh_start = heh1;
        auto heh_current = heh_start;
        int sign = 1;
        do
        {
            if (heh_current == heh1)
                sign = 1;
            else if (heh_current == _mesh.opposite_halfedge_handle(heh0))
                sign = -1;
            else
                result2 += sign * _b[heh_current];

            // Rotate CCW around common vertex
            heh_current = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(heh_current));
        }
        while (heh_current != heh_start);
    }

    // If the path is not closed: Measure intersections at start and end vertex
    if (!_a.is_closed(_mesh))
    {
        const auto v_path_start = _a.vh_start(_mesh);
        const auto v_path_end = _a.vh_end(_mesh);
        ISM_ASSERT_NEQ(v_path_start, v_path_end);

        for (auto v : { v_path_start, v_path_end })
        {
            // Endpoint not on a cycle. Nothing to do.
            if (!is_on_cycle(_mesh, _b, v))
                continue;

            // Path endpoints with other than two incident cycle edges are not supported.
            // TODO: How to implement valence > 2 correctly?
            ISM_ASSERT_EQ(cycle_valence(_mesh, _b, v), 2);

            const auto heh_start = (v == v_path_start) ?
                        _a.hehs.front()
                      : _mesh.opposite_halfedge_handle(_a.hehs.back());

            // Last edge runs along the cycle. Nothing to do.
            if (_b[heh_start] != 0)
                continue;

            // Find first cycle halfedge edge in ccw order
            const int sign = (v == v_path_start) ? -1 : 1;
            auto heh_current = heh_start;
            do
            {
                if (_b[heh_current] != 0)
                {
                    result2 += sign * _b[heh_current];
                    break;
                }

                // Rotate CW around vertex
                heh_current = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(heh_current));
            }
            while (heh_current != heh_start);
        }
    }

    if (is_on_cycle(_mesh, _b, _a.vh_start(_mesh)) && is_on_cycle(_mesh, _b, _a.vh_end(_mesh)))
        ISM_ASSERT_EQ(result2 % 2, 0)
    else if (!is_on_cycle(_mesh, _b, _a.vh_start(_mesh)) && !is_on_cycle(_mesh, _b, _a.vh_end(_mesh)))
        ISM_ASSERT_EQ(result2 % 2, 0)
    else
        ISM_ASSERT_NEQ(result2 % 2, 0)

    return (double)result2 / 2.0;
}

int
intersection_number(
        const TriMesh& _mesh,
        const PrimalLoop& _a,
        const PrimalCycle& _b)
{
    if (_a.hehs.empty())
        return 0;

    // Loop _a may only start / end on a cycle if it is closed.
    if (!_a.is_closed(_mesh))
    {
        ISM_ASSERT(!is_on_cycle(_mesh, _b, _a.vh_start(_mesh)));
        ISM_ASSERT(!is_on_cycle(_mesh, _b, _a.vh_end(_mesh)));
    }

    const double ip_d = intersection_number_double(_mesh, _a, _b);
    ISM_ASSERT(is_integer(ip_d));

    return to_integer(ip_d);
}

int
intersection_number(
        const TriMesh& _mesh,
        const PrimalLoop& _a,
        const PrimalLoop& _b)
{
    const PrimalCycle b_cycle = to_cycle(_mesh, _b);
    return intersection_number(_mesh, _a, b_cycle);
}

MatXi
intersection_form(
        const TriMesh& _mesh,
        const PrimalLoops& _loops)
{
    const int n = _loops.size();
    MatXi result(n, n);
    for (int row = 0; row < n; ++row)
    {
        for (int col = 0; col < n; ++col)
        {
            result(row, col) = intersection_number(_mesh, _loops[row], _loops[col]);
        }
    }
    return result;
}

}
