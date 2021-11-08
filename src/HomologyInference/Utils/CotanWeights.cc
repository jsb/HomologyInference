/*
 * Author: Patrick Schmidt
 */

#include <HomologyInference/Utils/CotanWeights.hh>

namespace HomologyInference
{

double cotan_weight(
        const TriMesh& _mesh,
        const HEH _heh)
{
    HEH h0 = _heh;
    HEH h1 = _mesh.opposite_halfedge_handle(_heh);

    if (_mesh.is_boundary(h0) || _mesh.is_boundary(h1))
        ISM_ERROR_throw("Cotan weight not defined for boundary edge.")

    double weight = 0.0;

    {
        Vec3d e0, e1;
        h0 = _mesh.next_halfedge_handle(h0);
        _mesh.calc_edge_vector(h0, e0);
        e0 *= -1.0;
        h0 = _mesh.next_halfedge_handle(h0);
        _mesh.calc_edge_vector(h0, e1);

        weight += (e0.dot(e1)) / e0.cross(e1).norm();
    }

    {
        Vec3d e0, e1;
        h1 = _mesh.next_halfedge_handle(h1);
        _mesh.calc_edge_vector(h1, e0);
        e0 *= -1.0;
        h1 = _mesh.next_halfedge_handle(h1);
        _mesh.calc_edge_vector(h1, e1);

        weight += (e0.dot(e1)) / e0.cross(e1).norm();
    }

    if (!std::isfinite(weight))
        ISM_ERROR_throw("Cotan weight infinite.")

    return weight;
}

}
