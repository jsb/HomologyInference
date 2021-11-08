#!/bin/bash

# Install dependencies in a first step
sudo apt install cmake libeigen3-dev libgl1-mesa-dev mesa-utils libglfw3 libglfw3-dev libxinerama-dev libxcursor-dev libxi-dev

# Perform a check to see whether GUROBI_HOME was set
if [ -z ${GUROBI_HOME+x} ]
then
    echo "GUROBI_HOME is unset. Please set this environment variable accordingly."
    exit 1
else
    echo "GUROBI_HOME is set to '$GUROBI_HOME'. Proceed with the installation."
fi

# install instructions
mkdir build
cd build
cmake ..
make -j6
MAKE_EXIT_CODE=$?
if [ $MAKE_EXIT_CODE -eq 0 ]
then
    echo "Build process successful."
else
    echo "Build process failed. Please consult the Troubleshooting section in the README.md file."
fi

