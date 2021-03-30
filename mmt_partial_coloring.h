#ifndef MMT_PARTIAL_COLORING_H_
#define MMT_PARTIAL_COLORING_H_

#include "mmt_graph.h"

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <map>
#include <queue>
#include <tuple>
#include <stdexcept>

using color = uint32_t;
using measure = int;


class MMTPartialColoring {
public:
  // empty constructor
  MMTPartialColoring(const int k, const MMTGraph * igraph, int L, int T);

  void crossover(const MMTPartialColoring* S1, const MMTPartialColoring* S2);

  void tabuSearch();

  void greedy();

  void dsatur();

  void toString(int maxLines = 10) const ;

private:
  const int k;

  const MMTGraph * graph;

  const int L, T;

  std::unordered_map<nodeid,color> colors;

  std::unordered_set<nodeid> uncolored;

  std::vector<std::unordered_set<nodeid> > color_classes;

  measure evaluate() const ;

  bool isValidColor(color value) const ;

  void setColor(nodeid u, color c) ;

  void moveToColor(nodeid u, color c) ;

  void lockColoring() ;

  std::tuple<const MMTPartialColoring*, std::vector<int>*, int*, const MMTPartialColoring*, std::vector<int>*, int* > selectParent(const MMTPartialColoring* s1, const MMTPartialColoring* s2, std::vector<int>* s1_c, std::vector<int>* s2_c, int* s1_n, int* s2_n, int cur_color);

  int dsatur_selectMaxNode(std::vector<std::pair<int, int> >& degs) const ;

  void dsatur_updateSatDeg(nodeid u, std::vector<std::pair<int, int> >& degs);

  int findMinAvailableColor(nodeid u);

  struct UInt32PairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const ;
  };
};

#endif
