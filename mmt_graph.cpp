#include "mmt_graph.h"

MMTGraph::MMTGraph(const int pncount, const int pecount, int **pelist) : n(pncount), m(pecount), adjList(pncount,std::unordered_set<nodeid>()) {
  for (size_t i = 0; i < pecount; i++) {
    adjList[(*pelist)[2*i]].insert((*pelist)[2*i + 1]);
    adjList[(*pelist)[2*i + 1]].insert((*pelist)[2*i]);
  }
}

bool MMTGraph::isAdj(const nodeid u, const nodeid v) const {
  assert(u != v && isValid(u) && isValid(v));
  return adjList[u].find(v) != adjList[u].end();
}

const std::unordered_set<nodeid>* MMTGraph::getNeighbors(const nodeid u) const {
  assert(isValid(u));
  return &adjList[u];
}

int MMTGraph::getDegree(const nodeid u) const {
  assert(isValid(u));
  return adjList[u].size();
}

void MMTGraph::toString(int maxLines, bool real) const
{
  std::cout << "MMTGraph.toString() of " << this << '\n';
  assert(n != 0 && m != 0);
  std::cout << "n = " << n << " : m = " << m << '\n';
  for (auto u = 0; u < n; u++) {
    if(maxLines == 0) {
      std::cout << "\t..." << '\n';
      return;
    }
    for (auto &v : adjList[u]) {
      if (real || u < v) {
         std::cout << u << " " << v << '\t';
      }
    }
    std::cout << "("<< getDegree(u) <<")" << '\n';
    maxLines--;
  }
}

bool MMTGraph::isValid(const nodeid u) const {
    return u >= 0 && u < n;
}
