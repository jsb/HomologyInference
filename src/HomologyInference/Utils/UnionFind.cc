/*
 * Author: Janis Born
 */
#include "UnionFind.hh"

#include <HomologyInference/Utils/Out.hh>
#include <numeric>

namespace HomologyInference {

UnionFind::UnionFind(int num_items)
{
    reset(num_items);
}

void UnionFind::reset(int num_items)
{
    parents.resize(num_items);
    std::iota(parents.begin(), parents.end(), 0);
}

void UnionFind::clear()
{
    reset(0);
}

int UnionFind::add()
{
    parents.push_back(parents.size());
    return parents.back();
}

int UnionFind::representative(int x)
{
    ISM_ASSERT_GEQ(x, 0);
    ISM_ASSERT_L(x, parents.size());
    if (x != parents[x]) {
        parents[x] = representative(parents[x]);
    }
    return parents[x];
}

int UnionFind::representative(int x) const
{
    ISM_ASSERT_GEQ(x, 0);
    ISM_ASSERT_L(x, parents.size());
    if (x == parents[x]) {
        return x;
    }
    else {
        return representative(parents[x]);
    }
}

void UnionFind::merge(int x, int y)
{
    int rootx = representative(x);
    int rooty = representative(y);
    parents[rootx] = rooty;
}

bool UnionFind::equivalent(int x, int y)
{
    return representative(x) == representative(y);
}

bool UnionFind::equivalent(int x, int y) const
{
    return representative(x) == representative(y);
}

}
