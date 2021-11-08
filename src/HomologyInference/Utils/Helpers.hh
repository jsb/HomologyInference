/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

inline OpenMesh::Vec2d to_acg2(const Vec2d& _v)
{
    ISM_ASSERT_EQ(_v.rows(), 2);
    return OpenMesh::Vec2d(_v[0], _v[1]);
}

template <typename VecT>
inline OpenMesh::Vec3d to_acg3(const VecT& _v)
{
    return OpenMesh::Vec3d(_v[0], _v[1], _v[2]);
}

inline Vec2d to_eigen(const OpenMesh::Vec2d& _v)
{
    return Vec2d(_v[0], _v[1]);
}

inline Vec3d to_eigen(const OpenMesh::Vec3d& _v)
{
    return Vec3d(_v[0], _v[1], _v[2]);
}

template <typename T, typename Vec2Type>
inline Vec2<T> to_vec2(const Vec2Type& _v)
{
    return Vec2<T>(_v[0], _v[1]);
}

template <typename T, typename Vec3Type>
inline Vec3<T> to_vec3(const Vec3Type& _v)
{
    return Vec3<T>(_v[0], _v[1], _v[2]);
}

template <typename T, typename Vec4Type>
inline Vec4<T> to_vec4(const Vec4Type& _v)
{
    return Vec4<T>(_v[0], _v[1], _v[2], _v[3]);
}

template <typename T, typename Mat2Type>
inline Mat2<T> to_mat2(const Mat2Type& _M)
{
    Mat2<T> M;
    M << _M(0, 0), _M(0, 1),
         _M(1, 0), _M(1, 1);

    return M;
}

/// Construct matrix from column vectors
template <typename T>
Mat2<T> mat2(
        const Vec2<T>& _a,
        const Vec2<T>& _b)
{
    Mat2<T> M;
    M << _a, _b;

    return M;
}

template <typename T>
inline Vec2<T> drop_z(const Vec3<T>& _v)
{
    return Vec2<T>(_v[0], _v[1]);
}

template <typename T>
inline T sqr(const T& _x)
{
    return _x * _x;
}

inline bool equals(const VecXd& _v1, const VecXd& _v2)
{
    if (_v1.size() != _v2.size())
        return false;

    for (int i = 0; i < _v1.size(); ++i)
    {
        if (_v1[i] != _v2[i])
            return false;
    }

    return true;
}

inline SparseMatrix Id(const size_t _n)
{
    SparseMatrix Id(_n, _n);
    Id.setIdentity();

    return Id;
}

int rank(const SparseMatrix& _M);

inline void handles(
        const TriMesh& _mesh, const HEH _heh,
        VH& _out_vh_a, VH& _out_vh_b, VH& _out_vh_c)
{
    auto heh = _heh;
    _out_vh_a = _mesh.to_vertex_handle(heh);
    heh = _mesh.next_halfedge_handle(heh);
    _out_vh_b = _mesh.to_vertex_handle(heh);
    heh = _mesh.next_halfedge_handle(heh);
    _out_vh_c = _mesh.to_vertex_handle(heh);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh,
        SVH& _out_v_a, SVH& _out_v_b, SVH& _out_v_c)
{
    auto h = OpenMesh::make_smart(_heh, _mesh);
    _out_v_a = h.to();
    h = h.next();
    _out_v_b = h.to();
    h = h.next();
    _out_v_c = h.to();
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        VH& _out_vh_a, VH& _out_vh_b, VH& _out_vh_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_vh_a, _out_vh_b, _out_vh_c);
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        SVH& _out_v_a, SVH& _out_v_b, SVH& _out_v_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_v_a, _out_v_b, _out_v_c);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh,
        VH& _out_vh_a, VH& _out_vh_b, VH& _out_vh_c,
        HEH& _out_heh_a, HEH& _out_heh_b, HEH& _out_heh_c)
{
    _out_heh_a = _heh;
    _out_vh_a = _mesh.to_vertex_handle(_out_heh_a);
    _out_heh_b = _mesh.next_halfedge_handle(_out_heh_a);
    _out_vh_b = _mesh.to_vertex_handle(_out_heh_b);
    _out_heh_c = _mesh.next_halfedge_handle(_out_heh_b);
    _out_vh_c = _mesh.to_vertex_handle(_out_heh_c);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh,
        SVH& _out_v_a, SVH& _out_v_b, SVH& _out_v_c,
        SHEH& _out_h_a, SHEH& _out_h_b, SHEH& _out_h_c)
{
    _out_h_a = OpenMesh::make_smart(_heh, _mesh);
    _out_v_a = _out_h_a.to();
    _out_h_b = _out_h_a.next();
    _out_v_b = _out_h_b.to();
    _out_h_c = _out_h_b.next();
    _out_v_c = _out_h_c.to();
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        VH& _out_vh_a, VH& _out_vh_b, VH& _out_vh_c,
        HEH& _out_heh_a, HEH& _out_heh_b, HEH& _out_heh_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_vh_a, _out_vh_b, _out_vh_c, _out_heh_a, _out_heh_b, _out_heh_c);
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        SVH& _out_v_a, SVH& _out_v_b, SVH& _out_v_c,
        SHEH& _out_h_a, SHEH& _out_h_b, SHEH& _out_h_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_v_a, _out_v_b, _out_v_c, _out_h_a, _out_h_b, _out_h_c);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh_to_a,
        HEH& _out_heh_a, HEH& _out_heh_b, HEH& _out_heh_c)
{
    _out_heh_a = _heh_to_a;
    _out_heh_b = _mesh.next_halfedge_handle(_out_heh_a);
    _out_heh_c = _mesh.next_halfedge_handle(_out_heh_b);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh_to_a,
        SHEH& _out_h_a, SHEH& _out_h_b, SHEH& _out_h_c)
{
    _out_h_a = OpenMesh::make_smart(_heh_to_a, _mesh);
    _out_h_b = _out_h_a.next();
    _out_h_c = _out_h_b.next();
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        HEH& _out_heh_a, HEH& _out_heh_b, HEH& _out_heh_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_heh_a, _out_heh_b, _out_heh_c);
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        SHEH& _out_h_a, SHEH& _out_h_b, SHEH& _out_h_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_h_a, _out_h_b, _out_h_c);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh_to_a,
        EH& _out_eh_a, EH& _out_eh_b, EH& _out_eh_c)
{
    HEH heh_a, heh_b, heh_c;
    handles(_mesh, _heh_to_a, heh_a, heh_b, heh_c);
    _out_eh_a = _mesh.edge_handle(heh_a);
    _out_eh_b = _mesh.edge_handle(heh_b);
    _out_eh_c = _mesh.edge_handle(heh_c);
}

inline void handles(
        const TriMesh& _mesh, const HEH _heh_to_a,
        SEH& _out_e_a, SEH& _out_e_b, SEH& _out_e_c)
{
    auto h = OpenMesh::make_smart(_heh_to_a, _mesh);
    _out_e_a = h.edge();
    h = h.next();
    _out_e_b = h.edge();
    h = h.next();
    _out_e_c = h.edge();
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        EH& _out_eh_a, EH& _out_eh_b, EH& _out_eh_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_eh_a, _out_eh_b, _out_eh_c);
}

inline void handles(
        const TriMesh& _mesh, const FH _fh,
        SEH& _out_e_a, SEH& _out_e_b, SEH& _out_e_c)
{
    const auto heh = _mesh.halfedge_handle(_fh);
    handles(_mesh, heh, _out_e_a, _out_e_b, _out_e_c);
}

inline void points(
        const TriMesh& _mesh, const HEH _heh,
        Vec3d& _out_p_a, Vec3d& _out_p_b, Vec3d& _out_p_c)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _heh, vh_a, vh_b, vh_c);

    _out_p_a = _mesh.point(vh_a);
    _out_p_b = _mesh.point(vh_b);
    _out_p_c = _mesh.point(vh_c);
}

inline void points(
        const TriMesh& _mesh, const FH _fh,
        Vec3d& _out_p_a, Vec3d& _out_p_b, Vec3d& _out_p_c)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _fh, vh_a, vh_b, vh_c);

    _out_p_a = _mesh.point(vh_a);
    _out_p_b = _mesh.point(vh_b);
    _out_p_c = _mesh.point(vh_c);
}

template <int d, typename T>
void points(
        const TriMesh& _mesh, const ExternalProperty<VH, Vec<d, T>>& _embedding, const HEH _heh,
        Vec<d, T>& _out_p_a, Vec<d, T>& _out_p_b, Vec<d, T>& _out_p_c)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _heh, vh_a, vh_b, vh_c);

    _out_p_a = _embedding[vh_a];
    _out_p_b = _embedding[vh_b];
    _out_p_c = _embedding[vh_c];
}

template <int d, typename T>
void points(
        const TriMesh& _mesh, const ExternalProperty<VH, Vec<d, T>>& _embedding, const FH _fh,
        Vec<d, T>& _out_p_a, Vec<d, T>& _out_p_b, Vec<d, T>& _out_p_c)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _fh, vh_a, vh_b, vh_c);

    _out_p_a = _embedding[vh_a];
    _out_p_b = _embedding[vh_b];
    _out_p_c = _embedding[vh_c];
}

double total_area(
        const TriMesh& _mesh);

std::pair<Vec3d, Vec3d> bounding_box(
        const TriMesh& _mesh);

double bounding_box_diagonal(
        const TriMesh& _mesh);

Color log_color(
        const double _val,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color);

bool incident(
        const TriMesh& _mesh,
        const VH _vh,
        const FH _fh);

template <typename T>
bool pairwise_distinct(
        const std::vector<T>& _v)
{
    auto v = _v; // copy
    std::sort(v.begin(), v.end());

    return std::unique(v.begin(), v.end()) == v.end();
}

VH closest_vertex(
        const TriMesh& _mesh,
        const Vec3d& _p);

}
