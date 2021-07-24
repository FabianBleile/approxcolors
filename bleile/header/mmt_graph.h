#ifndef MMT_GRAPH_H_
#define MMT_GRAPH_H_

extern "C" {
  #include "bleile/utils/mmt_read.h"
}

#include <iostream>
#include <vector>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>

using nodeid = int;

class MMTGraph {
public:
  int n, m;
  float dens;

  std::string instance;

  MMTGraph(int argc, char **av) ;

  MMTGraph(int ncount, int ecount, int *elist) ;

  MMTGraph(MMTGraph * input) ;

  void initFromElist(int ncount, int ecount, int *elist);

  bool isAdj(const nodeid u, const nodeid v) const ;

  const std::vector<nodeid>* getNeighbors(const nodeid u) const ;

  int getDegree(const nodeid u) const ;

  void toString(int maxLines = 5, bool real = true) const ;

private:
  std::vector<std::vector<nodeid> > adjList;

  bool isValid(const nodeid u) const;
};

#endif
