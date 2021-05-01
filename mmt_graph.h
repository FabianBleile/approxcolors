#ifndef MMT_GRAPH_H_
#define MMT_GRAPH_H_

extern "C" {
  #include "mmt_read.h"
}

#include <iostream>
#include <vector>
#include <unordered_set>

using nodeid = uint32_t;

class MMTGraph {
public:
  int n, m;

  MMTGraph(int argc, char **av) ;

  bool isAdj(const nodeid u, const nodeid v) const ;

  const std::vector<nodeid>* getNeighbors(const nodeid u) const ;

  int getDegree(const nodeid u) const ;

  void toString(int maxLines = 5, bool real = true) const ;

private:
  std::vector<std::vector<nodeid> > adjList;

  bool isValid(const nodeid u) const;
};

#endif
