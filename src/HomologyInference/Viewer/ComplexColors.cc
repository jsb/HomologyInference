#include "ComplexColors.hh"

#include <glow-extras/colors/color.hh>

namespace HomologyInference
{

Color complex_color(const Complex& _x)
{
    auto tg_color = complex_color_tg(_x);
    return Color(tg_color.r, tg_color.g, tg_color.b, 1.0);
}

tg::color3 complex_color_tg(const Complex& _x)
{
    const double angle = std::arg(_x);
    const double magnitude = std::abs(_x);

    const double angle_deg = angle * 360 / (2 * M_PI);
    const double lightness = std::atan(magnitude) * 2 / M_PI;

    auto glow_color = glow::colors::color::from_hsl(angle_deg, 1.0, lightness);
    return tg::color3(glow_color.r, glow_color.g, glow_color.b);
}

}
