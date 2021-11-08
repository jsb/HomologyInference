#include "XYZColors.hh"

#include <HomologyInference/Viewer/RWTHColors.hh>

namespace HomologyInference
{

ExternalProperty<VH, Color>
xyz_colors(
        const TriMesh& _mesh,
        const Vec3d& _x,
        const Vec3d& _y,
        const Vec3d& _z,
        const ExternalProperty<VH, BarycentricPoint>& _vtpm)
{
    ExternalProperty<VH, Color> result(_mesh, RWTH_BLACK);

    // Only color vertices at which vtpm is defined
    const int n = _vtpm.empty() ?
                _mesh.n_vertices()
              : _mesh.vertices().count_if([&] (auto v) { return _vtpm[v].is_valid(); });

    MatXd xyz(n, 3);
    {
        int row = 0;
        for (const auto& vh : _mesh.vertices())
        {
            if (_vtpm.empty() || _vtpm[vh].is_valid())
            {
                const VecXd& p = _mesh.point(vh);
                xyz(row, 0) = _x.dot(p);
                xyz(row, 1) = _y.dot(p);
                xyz(row, 2) = _z.dot(p);
                ++row;
            }
        }
    }

    const VecXd min_xyz = xyz.colwise().minCoeff();
    const VecXd max_xyz = xyz.colwise().maxCoeff();
    const VecXd xyz_range = max_xyz - min_xyz;

    {
        int row = 0;
        for (const auto& vh : _mesh.vertices())
        {
            if (_vtpm.empty() || _vtpm[vh].is_valid())
            {
                Color c = Color(0,0,0,1);
                for (int col = 0; col < 3; ++col)
                    c[col] = (xyz(row, col) - min_xyz[col]) / xyz_range[col];
                result[vh] = c;
                ++row;
            }
        }
    }

    return result;
}

}
