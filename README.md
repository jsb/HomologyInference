# HomologyInference

This is a prototype implementation of the paper [Surface Map Homology Inference](https://www.graphics.rwth-aachen.de/publication/03335/) by Janis Born, Patrick Schmidt, Marcel Campen, and Leif Kobbelt.

![Teaser](teaser.png)

This repository contains:
* The core `HomologyInference` library.
* Example applications that replicate figures from the SGP 2021 paper.

## Prerequisites

### Dependencies

There are several third-party dependencies.
On a Debian-based Linux system, you can install all of them using the following command:

    sudo apt install cmake libeigen3-dev libgl1-mesa-dev mesa-utils libglfw3 libglfw3-dev libxinerama-dev libxcursor-dev libxi-dev git g++

### Gurobi

This project depends on the **Gurobi Optimizer** library.
Please first obtain a valid Gurobi license.
Gurobi offers [free academic licenses](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).

The installation and setup process for Linux is described in the [Software Installation Guide](https://www.gurobi.com/documentation/9.1/quickstart_linux/index.html).

Depending on your compiler, it may be necessary to recompile the Gurobi C++ interface to avoid ABI incompatibilities.
Follow [these instructions](https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-).

## Build Instructions

Make sure to check out this repository with all submodule dependencies.
Either initially clone with

    git clone --recursive https://github.com/jsb/HomologyInference.git

or check out the submodules afterwards using

    git submodule update --init --recursive

Before building, please set the `GUROBI_HOME` environment variable of your shell to the path containing the `bin` directory from your Gurobi installation, e. g.

    export GUROBI_HOME=/home/username/opt/gurobi911/linux64

Navigate to the `HomologyInference` directory and use the following commands to build the project:

    mkdir build
    cd build
    cmake ..
    make -j4

A script that automates the above steps (except the installation of Gurobi) has been provided for **Ubuntu Linux** (tested on Ubuntu 20.04 LTS).
It will download the required dependencies and automatically build the project.
Simply run

    ./install.sh
    
If any problems occur during the build process, please consult the Troubleshooting section below.

## Run Instructions

The `build` directory should now contain several binary files:

* `bob_donut` (Figure 5)
* `dog_elephant` (Figure 17)
* `double_torus_twists` (Figure 2)
* `functional_maps` (Figure 8)
* `landmarks` (Figure 9)
* `mother_and_child` (Figure 1)
* `nn_projection` (Figure 7)
* `pretzel_three_hole` (Figure 11)

Each application demonstrates a map inference example and will open one or several viewer widgets.

Homology classes are indicated by representative cycles (colored loops).
Their shape can be modified by clicking and dragging the offset values next to the color swatches in the *Mesh A Cycles* and *Mesh B Cycles* toolboxes.

Viewer navigation:

* Drag left mouse button to rotate
* Double-click to set pivot point
* Scroll mouse wheel to zoom

Press Escape to close the widget. Some demos will open several widgets in succession.

## Troubleshooting

### Build Errors

The CMake configuration will fail if the `GUROBI_HOME` environment variable is not set (correctly).
On Linux, please follow the instructions in the [Quick Start Guide](https://www.gurobi.com/documentation/9.1/quickstart_linux/software_installation_guid.html) to set the correct path to your Gurobi home directory.
Alternatively, you can pass this path to CMake on the command line via an additional `-DGUROBI_HOME=...` argument.

Build errors of the form `undefined reference to GRBModel::addVar`... may indicate ABI incompatibilities of the supplied Gurobi C++ library with your compiler.
In this case, simply recompile the library according to [these instructions](https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-).

### Run-Time Errors

If you get the error message `terminate called after throwing an instance of 'GRBException'`, your Gurobi license might be expired or not properly installed.

## Authors and Contributors

* [Janis Born](https://www.graphics.rwth-aachen.de/person/97/)
* [Patrick Schmidt](https://www.graphics.rwth-aachen.de/person/232/)
* [Marcel Campen](http://graphics.cs.uos.de/)
* [Sasa Lukic](https://www.graphics.rwth-aachen.de/person/300/)
* [Leif Kobbelt](https://www.graphics.rwth-aachen.de/person/3/)

## License

This project is released under the **MIT License**.
