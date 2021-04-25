#ifndef MMT_PARTIAL_COLORING_H_
#define MMT_PARTIAL_COLORING_H_

#include "mmt_graph.h"

extern "C" {
  #include "hungarian.h"
}

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
  MMTPartialColoring(const int k, MMTGraph * igraph, int L, int T);

  bool crossover(MMTPartialColoring* S1, MMTPartialColoring* S2);

  bool tabuSearch();

  bool greedy();

  bool priorityGreedy(const std::vector<int>& v);

  bool dsatur();

  void lockColoring() ;

  measure evaluate() const ;

  int distanceTo(MMTPartialColoring* S, bool exact = false) ;

  int getNumColors() const ;

  void toString(int maxLines = 7) const ;

  std::unordered_set<nodeid> uncolored;

  std::vector<std::unordered_set<nodeid> > color_classes;

  int k;

private:

  MMTGraph * graph;

  int L, T;

  std::unordered_map<nodeid,color> colors;

  bool isValidColor(color value) const ;

  void setColor(nodeid u, color c) ;

  void moveToColor(nodeid u, color c) ;

  std::tuple<const MMTPartialColoring*, std::vector<int>*, int*, const MMTPartialColoring*, std::vector<int>*, int* > selectParent(const MMTPartialColoring* s1, const MMTPartialColoring* s2, std::vector<int>* s1_c, std::vector<int>* s2_c, int* s1_n, int* s2_n, int cur_color);

  int dsatur_selectMaxNode(std::vector<std::pair<int, int> >& degs) const ;

  void dsatur_updateSatDeg(nodeid u, std::vector<std::pair<int, int> >& degs);

  int findMinAvailableColor(nodeid u);

  int approxDistance(std::vector<std::vector<int> >& matIntersec);

  int exactDistance(std::vector<std::vector<int> >& matIntersec);

  struct UInt32PairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const ;
  };
};

#endif
