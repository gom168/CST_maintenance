#include "gadget.h"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <set>
#include <time.h>

#include "../defs.h"

namespace gadget {
void RepeatWith(const char symbol, const int repeat) {
  ASSERT(repeat > 0);
  for (int i = 0; i < repeat; ++i) {
    printf("%c", symbol);
  }
  putchar('\n');
}
std::vector<std::vector<int>> ReadGraph(const char* const path,
                                        int* const n, int* const m) {
  auto file = fopen(path, "r");
  fscanf(file, "%d %d", n, m);
  ASSERT(*n > 0 && *m > 0);
  std::vector<std::vector<int>> graph(*n);
  // for debug purpose
  std::set<std::pair<int, int>> debug;
  // read the edges
  while (true) {
    int v1, v2;
    fscanf(file, "%d %d", &v1, &v2);
    if (feof(file)) break;
    // for debug purpose
    ASSERT(0 <= v1 && v1 < *n && 0 <= v2 && v2 < *n && v1 != v2);
    if (v1 > v2) std::swap(v1, v2);
    ASSERT(debug.count(std::make_pair(v1, v2)) == 0);
    debug.insert(std::make_pair(v1, v2));

    graph[v1].push_back(v2);
    graph[v2].push_back(v1);
  }
  ASSERT(static_cast<decltype(debug.size())>(*m) == debug.size());
  fclose(file);
  return graph;
}
std::vector<std::pair<int, int>> ReadTempEdgesS(const char* const path,
                                                int* const n, int* const m) {
  auto file = fopen(path, "r");
  fscanf(file, "%d %d", n, m);
  printf("%d %d\n", *n, *m);
  ASSERT(*n > 0 && *m > 0);
  std::vector<std::pair<int, int>> edges;
  // for debug purpose
  std::set<std::pair<int, int>> debug;
  int prev_time = std::numeric_limits<int>::min();
  // read the edges
  while (true) {
    int v1, v2, tm;
    fscanf(file, "%d %d %d", &v1, &v2, &tm);
    if (feof(file)) break;
    // for debug purpose
    ASSERT(0 <= v1 && v1 < *n && 0 <= v2 && v2 < *n && v1 != v2);
    ASSERT(prev_time <= tm); // time is non-decreasing
    if (v1 > v2) std::swap(v1, v2);
    const auto e = std::make_pair(v1, v2);
    ASSERT(debug.count(e) == 0);
    debug.insert(e);
    // insert the edge
    edges.push_back(e);
    prev_time = tm;
    printf("%d %d %d\n", v1, v2, tm);
  }
  printf("%d\n", debug.size());
  ASSERT(static_cast<int>(debug.size()) == *m);
  ASSERT(static_cast<int>(edges.size()) == *m);
  fclose(file);
  return edges;
}
std::vector<std::pair<int, int>> ReadEdgesS(const char* const path,
                                            int* const n, int* const m) {
  srand(time(NULL));

  auto file = fopen(path, "r");
  fscanf(file, "%d %d", n, m);
  ASSERT(*n > 0 && *m > 0);
  std::vector<std::pair<int, int>> edges;
  // for debug purpose
  std::set<std::pair<int, int>> debug;
  // read the edges
  while (true) {
    int v1, v2;
    fscanf(file, "%d %d", &v1, &v2);
    if (feof(file)) break;
    // for debug purpose
    ASSERT(0 <= v1 && v1 < *n && 0 <= v2 && v2 < *n && v1 != v2);
    if (v1 > v2) std::swap(v1, v2);
    const auto e = std::make_pair(v1, v2);
    ASSERT(debug.count(e) == 0);
    debug.insert(e);
    // insert the edge
    edges.push_back(e);
  }
  ASSERT(static_cast<int>(debug.size()) == *m);
  ASSERT(static_cast<int>(edges.size()) == *m);
  // randomly shuffle the edges
  for (size_t e = 0; e < edges.size(); ++e) {
    const size_t e2 = e + rand() % (edges.size() - e);
    std::swap(edges[e], edges[e2]);
  }
  fclose(file);
  return edges;
}
}  // namespace gadget
