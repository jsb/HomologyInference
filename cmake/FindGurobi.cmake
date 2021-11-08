# Set the GUROBI_HOME variable to the directory containing Gurobi's lib and include directories.
set(GUROBI_HOME "$ENV{GUROBI_HOME}" CACHE PATH "Root directory of your Gurobi distribution")

find_library(
  GUROBI_LIBRARY
  NAMES gurobi gurobi81 gurobi90 gurobi91
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

MESSAGE("GUROBI_HOME: ${GUROBI_HOME}")
MESSAGE("GUROBI_LIBRARY: ${GUROBI_LIBRARY}")
MESSAGE("GUROBI_CPP_LIBRARY: ${GUROBI_CPP_LIBRARY}")
MESSAGE("GUROBI_INCLUDE_DIR: ${GUROBI_INCLUDE_DIR}")
MESSAGE("GUROBI_CPP_INCLUDE_DIR: ${GUROBI_CPP_INCLUDE_DIR}")


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
