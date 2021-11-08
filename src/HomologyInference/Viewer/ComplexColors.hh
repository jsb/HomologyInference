#pragma once

#include <HomologyInference/Types.hh>
#include <typed-geometry/types/color.hh>
#include <complex>

namespace HomologyInference
{

/// Maps angle to hue, magnitude to lightness (black towards 0, white towards Infty).
Color complex_color(const Complex& _x);

/// Maps angle to hue, magnitude to lightness (black towards 0, white towards Infty).
/// Result as a tg::color3.
tg::color3 complex_color_tg(const Complex& _x);

}
