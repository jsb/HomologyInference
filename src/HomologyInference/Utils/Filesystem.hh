/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <experimental/filesystem>

namespace HomologyInference
{

namespace fs = std::experimental::filesystem;
inline fs::path SOURCE_PATH = fs::path(SOURCE_PATH_STR);
inline fs::path DATA_PATH = fs::path(DATA_PATH_STR);
inline fs::path OUTPUT_PATH = fs::path(OUTPUT_PATH_STR);

}
