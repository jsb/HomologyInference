/*
 * Author: Janis Born
 */
#pragma once

#include <vector>

namespace HomologyInference {

/// Represents a partition of a set of items (each identified by an integer index) into equivalence classes.
/// The equivalence relation between items can be extended using merge(..., ...).
/// This equivalence will be propagated transitively.
struct UnionFind
{
    explicit UnionFind(int num_items = 0);

    /// Initializes a set of num_items items, each in its own equivalence class.
    void reset(int num_items);
    void clear();

    /// Adds a new item and returns its representative.
    /// Initially, this item is equivalent to no other item.
    int add();

    /// Merges the equivalence classes of items x and y.
    void merge(int x, int y);

    bool equivalent(int x, int y);
    bool equivalent(int x, int y) const;

    int representative(int x);
    int representative(int x) const;

    std::vector<int> parents;
};

}
