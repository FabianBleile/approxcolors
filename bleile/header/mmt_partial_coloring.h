#ifndef MMT_PARTIAL_COLORING_H_
#define MMT_PARTIAL_COLORING_H_

#include "../header/mmt_graph.h"
#include "../utils/pair_hasher.hpp"

extern "C" {
  #include "../utils/hungarian.h"
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

class PartialCol {
public:
  PartialCol(const int k, Graph * graph);

  Graph * graph;
  int k;
  std::vector<color> colors;
  std::unordered_set<nodeid> uncolored;

  int fitness;

  bool greedy(const vector<nodeid>& v = vector<nodeid>());

  bool dsatur();

  int distanceTo(PartialCol* S, bool exact = false) ;

  int evaluate() ;

  void setColor(nodeid u, color c) ;

  void moveToColor(nodeid u, color c) ;

  int getMinAvailableColor(nodeid u);

  int getNumColors() const ;

  void toString(int maxLines = 3) const ;

  bool operator < (const PartialCol& S) const ;

private:
  int dsatur_selectMaxNode(const std::vector<nodeid>& shuffled_nodes, std::vector<int>& satdegree, std::vector<int>& freedegree) const ;

  void dsatur_updateDeg(nodeid u, std::vector<int>& satdegree, std::vector<int>& freedegree);

  int approxDistance(std::vector<std::vector<double> >& matIntersec, int num_uncolored);

  int exactDistance(std::vector<std::vector<double> >& matIntersec, int num_uncolored);

};


class EvolPartialCol : public PartialCol {
public:
  // empty constructor
  EvolPartialCol(const int k, Graph * graph);

  // TODO : constructor from superclass

  bool crossover(EvolPartialCol& S1, EvolPartialCol& S2);

  bool tabuSearch(int L, int T);

  void setK(int k);

  int getOptimalColor(nodeid u, std::vector<std::vector<int>>& tabuList, int it);

  bool priorityGreedy(const std::vector<int>& v);

  void buildColorClasses() ;

  bool checkColoring() ;

  int getFitness();

  std::vector<std::unordered_set<nodeid> > color_classes;

private:

  std::tuple<const EvolPartialCol*,std::vector<int>*, int*,
        const EvolPartialCol*,std::vector<int>*, int* >
                                  selectParent(const EvolPartialCol* s1,
                                                const EvolPartialCol* s2,
                                                std::vector<int>* s1_c,
                                                std::vector<int>* s2_c,
                                                int* s1_n, int* s2_n, int cur_color
                                              );
};

#endif
