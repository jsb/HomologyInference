# Set the GUROBI_HOME variable to the directory containing Gurobi's lib and include directories.
set(GUROBI_HOME "$ENV{GUROBI_HOME}" CACHE PATH "Root directory of your Gurobi distribution")

if(NOT EXISTS ${GUROBI_HOME})
    message("GUROBI_HOME not set.")
else()
    string(REGEX MATCH "gurobi[0-9][0-9][0-9]" GUROBI_RELEASE_NAME ${GUROBI_HOME})
    string(REGEX MATCH "[0-9][0-9][0-9]" GUROBI_PATCH_VERSION ${GUROBI_RELEASE_NAME})
    string(SUBSTRING ${GUROBI_PATCH_VERSION} 0 2 GUROBI_MINOR_VERSION)
    set(GUROBI_LIBRARY_NAME "gurobi${GUROBI_MINOR_VERSION}")

    message("GUROBI_HOME: ${GUROBI_HOME}")
    message("GUROBI_RELEASE_NAME: ${GUROBI_RELEASE_NAME}")
    message("GUROBI_PATCH_VERSION: ${GUROBI_PATCH_VERSION}")
    message("GUROBI_MINOR_VERSION: ${GUROBI_MINOR_VERSION}")
    message("GUROBI_LIBRARY_NAME: ${GUROBI_LIBRARY_NAME}")

    find_library(
        GUROBI_LIBRARY
        NAMES gurobi ${GUROBI_LIBRARY_NAME}
        PATHS "${GUROBI_HOME}"
        PATH_SUFFIXES lib
    )

    find_library(
        GUROBI_CPP_LIBRARY
        NAMES gurobi_c++
        PATHS "${GUROBI_HOME}"
        PATH_SUFFIXES lib
    )

    find_path(
        GUROBI_INCLUDE_DIR
        NAMES gurobi_c.h
        PATHS "${GUROBI_HOME}"
        PATH_SUFFIXES include
    )

    find_path(
        GUROBI_CPP_INCLUDE_DIR
        NAMES gurobi_c++.h
        PATHS "${GUROBI_HOME}"
        PATH_SUFFIXES include
    )

    message("GUROBI_LIBRARY: ${GUROBI_LIBRARY}")
    message("GUROBI_CPP_LIBRARY: ${GUROBI_CPP_LIBRARY}")
    message("GUROBI_INCLUDE_DIR: ${GUROBI_INCLUDE_DIR}")
    message("GUROBI_CPP_INCLUDE_DIR: ${GUROBI_CPP_INCLUDE_DIR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    Gurobi
    DEFAULT_MSG
    GUROBI_LIBRARY
    GUROBI_INCLUDE_DIR
    GUROBI_CPP_LIBRARY
    GUROBI_CPP_INCLUDE_DIR
)

if(GUROBI_FOUND)
    # Gurobi library
    add_library(gurobi SHARED IMPORTED)
    set_target_properties(gurobi PROPERTIES IMPORTED_LOCATION ${GUROBI_LIBRARY})

    #target_include_directories(gurobi INTERFACE ${GUROBI_INCLUDE_DIR})
    set_property(TARGET gurobi APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${GUROBI_INCLUDE_DIR})

    # Gurobi C++ library
    add_library(gurobi_c++ STATIC IMPORTED)
    set_target_properties(gurobi_c++ PROPERTIES IMPORTED_LOCATION ${GUROBI_CPP_LIBRARY})

    #target_include_directories(gurobi INTERFACE ${GUROBI_CPP_INCLUDE_DIR})
    #target_link_libraries(gurobi_c++ INTERFACE gurobi)
    set_property(TARGET gurobi_c++ APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${GUROBI_CPP_INCLUDE_DIR})
    set_property(TARGET gurobi_c++ APPEND PROPERTY INTERFACE_LINK_LIBRARIES gurobi)
endif()

mark_as_advanced(
    GUROBI_LIBRARY
    GUROBI_INCLUDE_DIR
    GUROBI_CPP_LIBRARY
    GUROBI_CPP_INCLUDE_DIR
)
