#ifndef MMT_GRAPH_H_
#define MMT_GRAPH_H_

#include <iostream>
#include <vector>
#include <unordered_set>

using nodeid = uint32_t;

class MMTGraph {
public:
  const int n, m;

  MMTGraph(const int pncount, const int pecount, int **pelist) ;

  bool isAdj(const nodeid u, const nodeid v) const ;

  const std::unordered_set<nodeid>* getNeighbors(const nodeid u) const ;

  int getDegree(const nodeid u) const ;

  void toString(int maxLines = 5, bool real = true) const ;

private:
  std::vector<std::unordered_set<nodeid> > adjList;

  bool isValid(const nodeid u) const;
};

#endif
