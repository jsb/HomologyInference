#pragma once

#include <HomologyInference/Types.hh>

namespace HomologyInference
{

class RWTHColorGenerator
{
public:
    RWTHColorGenerator(int _n_skip = 0, int _cycle_size = -1)
        : current_index(_n_skip),
          cycle_size(_cycle_size) { }

    Color generate_next_color();
    std::vector<Color> generate_next_colors(int n);

private:
    int current_index = 0;
    int cycle_size = -1;
};

}
