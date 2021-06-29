#ifndef MMT_PARTIAL_COLORING_H_
#define MMT_PARTIAL_COLORING_H_

#include "bleile/header/mmt_graph.h"
#include "bleile/utils/pair_hasher.hpp"

extern "C" {
  #include "bleile/utils/hungarian.h"
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

class PartialColoring {
public:
  PartialColoring(const int k, MMTGraph * graph);

  MMTGraph * graph;
  int k;
  std::vector<color> colors;
  std::unordered_set<nodeid> uncolored;

  bool greedy();

  bool dsatur();

  int distanceTo(PartialColoring* S, bool exact = false) ;

  measure evaluate() const ;

  void setK(int k);

  void setColor(nodeid u, color c) ;

  void moveToColor(nodeid u, color c) ;

  int findMinAvailableColor(nodeid u);

  int getNumColors() const ;

  void toString(int maxLines = 14) const ;

private:
  int dsatur_selectMaxNode(const std::vector<nodeid>& shuffled_nodes, std::vector<int>& satdegree, std::vector<int>& freedegree) const ;

  void dsatur_updateDeg(nodeid u, std::vector<int>& satdegree, std::vector<int>& freedegree);

  int approxDistance(std::vector<std::vector<double> >& matIntersec, int num_uncolored);

  int exactDistance(std::vector<std::vector<double> >& matIntersec, int num_uncolored);

};


class MMTPartialColoring : public PartialColoring {
public:
  // empty constructor
  MMTPartialColoring(const int k, MMTGraph * graph, int L, int T);

  // TODO : constructor from superclass

  bool crossover(MMTPartialColoring* S1, MMTPartialColoring* S2);

  bool tabuSearch();

  bool tabuSearchSimplified();

  bool priorityGreedy(const std::vector<int>& v);

  void lockColoring() ;

  void checkColoring() ;

  std::vector<std::unordered_set<nodeid> > color_classes;

private:

  int L, T;

  std::tuple<const MMTPartialColoring*,std::vector<int>*, int*,
        const MMTPartialColoring*,std::vector<int>*, int* >
                                  selectParent(const MMTPartialColoring* s1,
                                                const MMTPartialColoring* s2,
                                                std::vector<int>* s1_c,
                                                std::vector<int>* s2_c,
                                                int* s1_n, int* s2_n, int cur_color
                                              );
};

#endif
