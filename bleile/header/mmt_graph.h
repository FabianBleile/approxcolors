#ifndef MMT_GRAPH_H_
#define MMT_GRAPH_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>

using nodeid = int;

class Graph {
public:
  int n = 0, m = 0;
  float dens;

  std::string instance;

  Graph(int argc, char **av) ;

  Graph(int ncount, int ecount, int *elist) ;

  Graph(Graph * input) ;

  bool writeToElist(char *f, int *pncount, int *pecount, int **pelist);

  void initFromElist(int *elist);

  bool isAdj(const nodeid u, const nodeid v) const ;

  const std::vector<nodeid>* getNeighbors(const nodeid u) const ;

  int getDegree(const nodeid u) const ;

  void toString(int maxLines = 5, bool real = true) const ;

private:
  std::vector<std::vector<nodeid> > adjList;

  bool isValidNode(const nodeid u) const;
};

#endif
