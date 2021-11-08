/*
 * Author: Patrick Schmidt
 */

#include "DualPath.hh"

namespace HomologyInference
{

DualPath::DualPath(const FH _fh_start, const FH _fh_end) :
    fh_start(_fh_start),
    fh_end(_fh_end)
{
}

int DualPath::length() const
{
    ISM_ASSERT(fh_start.is_valid());
    ISM_ASSERT(fh_end.is_valid());
    return (int)hehs.size() + 1;
}

bool DualPath::is_closed(const TriMesh& _mesh) const
{
    return (fh_start == fh_end);
}

SFH DualPath::face(
        const int _idx,
        const TriMesh& _mesh) const
{
    ISM_ASSERT_GEQ(_idx, 0);
    if (_idx < (int)hehs.size())
        return OpenMesh::make_smart(_mesh.face_handle(hehs[_idx]), _mesh);
    else if (_idx == (int)hehs.size())
        return OpenMesh::make_smart(fh_end, _mesh);
    else
        ISM_ERROR_throw("Index out of range.");
}

double DualPath::embedded_length(const TriMesh& _mesh) const
{
    double result = 0.0;
    for (const auto& heh : hehs)
    {
        // Approximate the intrinsic length of the dual edge by two segements:
        // Face centroid (p0) <-> edge midpoint (p1) <-> face centroid (p2).
        const FH fh_from = _mesh.face_handle(heh);
        const EH eh = _mesh.edge_handle(heh);
        const FH fh_to = _mesh.opposite_face_handle(heh);

        const Vec3d p1 = _mesh.calc_edge_midpoint(eh);

        Vec3d p0 = p1;
        if (fh_from.is_valid())
            p0 = _mesh.calc_face_centroid(fh_from);

        Vec3d p2 = p1;
        if (fh_to.is_valid())
            p2 = _mesh.calc_face_centroid(fh_to);

        const double edge_length = (p1 - p0).norm() + (p2 - p1).norm();
        result += edge_length;
    }
    return result;
}

DualPath DualPath::reversed(const TriMesh& _mesh) const
{
    DualPath path_rev(fh_end, fh_start);

    const int n = (int)hehs.size();
    path_rev.hehs.resize(n);
    for (int i = 0; i < n; ++i)
        path_rev.hehs[i] = _mesh.opposite_halfedge_handle(hehs[n-i-1]);

    assert_valid(_mesh);
    return path_rev;
}

void DualPath::assert_start_end_valid(const TriMesh& _mesh) const
{
    if (hehs.empty())
    {
        ISM_ASSERT_EQ(fh_start, fh_end);
    }
    else
    {
        if (fh_start.is_valid())
            ISM_ASSERT_EQ(_mesh.face_handle(hehs.front()), fh_start);
        if (fh_end.is_valid())
            ISM_ASSERT_EQ(_mesh.opposite_face_handle(hehs.back()), fh_end);
    }
}

void DualPath::assert_valid(const TriMesh& _mesh) const
{
    assert_start_end_valid(_mesh);
    FH current_fh = fh_start;
    for (const HEH heh : hehs)
    {
        if (_mesh.is_boundary(heh))
        {
            ISM_ASSERT(!current_fh.is_valid());
        }
        else
        {
            ISM_ASSERT(_mesh.face_handle(heh) == current_fh);
        }
        current_fh = _mesh.opposite_face_handle(heh);
    }
    ISM_ASSERT(current_fh == fh_end);
}

DualPath concatenate(const std::vector<DualPath>& _paths)
{
    ISM_ASSERT_GEQ((int)_paths.size(), 2);

    int n_hehs = 0;
    for (const auto& path : _paths)
        n_hehs += path.hehs.size();

    DualPath result(_paths.front().fh_start, _paths.back().fh_end);
    result.hehs.reserve(n_hehs);

    FH prev_end(-1);
    for (const auto& path : _paths)
    {
        ISM_ASSERT(path.fh_start.is_valid());
        ISM_ASSERT(path.fh_end.is_valid());
        ISM_ASSERT(path.fh_start == prev_end || !prev_end.is_valid());

        result.hehs.insert(result.hehs.end(), path.hehs.begin(), path.hehs.end());

        prev_end = path.fh_end;
    }

    return result;
}

bool self_intersecting(
        const DualPath& _path,
        const TriMesh& _mesh)
{
    ExternalProperty<FH, bool> visited(_mesh, false);

    for (int i = 0; i < _path.length(); ++i)
    {
        const auto f = _path.face(i, _mesh);
        ISM_ASSERT(_mesh.is_valid_handle(f));

        if (visited[f])
            return true;

        visited[f] = true;
    }

    return false;
}

DualPath until_self_intresection(
        const DualPath& _path,
        const TriMesh& _mesh)
{
    ExternalProperty<FH, bool> visited(_mesh, false);

    for (int i = 0; i < _path.length(); ++i)
    {
        const auto f = _path.face(i, _mesh);
        ISM_ASSERT(_mesh.is_valid_handle(f));

        if (visited[f])
        {
            ISM_ASSERT_GEQ(i, 1);

            // Face i is the first self-intersecting face.
            // Return prefix with last face being i-1.
            DualPath prefix(_path.fh_start, _path.face(i-1, _mesh));
            prefix.hehs.reserve(i-1);
            for (int j = 0; j < i-1; ++j)
                prefix.hehs.push_back(_path.hehs[j]);

            prefix.assert_valid(_mesh);
            return prefix;
        }

        visited[f] = true;
    }

    // No self-intersection, return input path.
    return _path;
}

}
