/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <string>
#include <HomologyInference/Extern/StopWatch.hh>
#include "Out.hh"

namespace HomologyInference
{

class Timer
{
public:
    Timer(
            const std::string& _name,
            const bool _silent = false)
        : name(_name),
          silent(_silent)
    {
        sw.start();
    }

    ~Timer()
    {
        if (!silent)
            ISM_INFO(ANSI_FG_WHITE << "[TIMER] " << name << " took " << sw.stop() / 1000.0 << "s.")
    }

private:
    const std::string name;
    StopWatch sw;
    bool silent;
};

}
