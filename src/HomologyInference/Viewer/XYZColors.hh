#pragma once

#include <HomologyInference/Types.hh>
#include <HomologyInference/Utils/BarycentricPoint.hh>

namespace HomologyInference
{

ExternalProperty<VH, Color>
xyz_colors(
        const TriMesh& _mesh,
        const Vec3d& _x = Vec3d(1,0,0),
        const Vec3d& _y = Vec3d(0,1,0),
        const Vec3d& _z = Vec3d(0,0,1),
        const ExternalProperty<VH, BarycentricPoint>& _vtpm = ExternalProperty<VH, BarycentricPoint>());

}
