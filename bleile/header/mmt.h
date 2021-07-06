#ifndef MMT_H_
#define MMT_H_

extern "C" {
  #include "lp.h"
  #include "color_defs.h"
}

#include "mmt_graph.h"
#include "mmt_partial_coloring.h"
#include "bleile/utils/pair_hasher.hpp"

#include "bleile/utils/color_cluster.cpp"

#include <time.h>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>

class MMT {
public:

  MMTGraph graph;

  MMT(MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, bool setBounds = false);

  enum status {
    UNSOLVED = 1,
    INIT_GREEDY = 2,
    INIT_DSATUR = 4,
    INIT_TABU = 8,
    EA = 16,
    COl_OPT = 32,
    EA_TIME_OUT = 64
  };

  void start();

  void PHASE0_EAInit();

  void PHASE1_EAOptimizer();

  status EADecision(int k);

  void PHASE2_ColumnOptimization();

  MMTPartialColoring* getColoring();

  std::stringstream streamLogs();

private:

  struct LogData {
    status status = UNSOLVED;
    int totTimeInSec = 0;
    int lastItTimeInSec = 0;
    int totNumOffsprings = 0;
    int lastItNumOffsprings = 0;
    int UB;
    int LB = 2;
    int colOpt = -1;
  };

  int L, T, time_limit_sec, pool_size, measure_best_solution;
  const int N;
  MMTPartialColoring cur_best_coloring;
  const double priority_noise = 0.5;
  int R; // pool spacing

  LogData logger;

  const int numColOpt = 1000;
  std::queue< std::unordered_set<nodeid> > columns;

  bool initPool(int k, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, int pool_size);

  void insertPool(MMTPartialColoring& new_individual, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority);

  void updatePool(MMTPartialColoring& new_individual, int old_individual,
    std::vector<MMTPartialColoring>& pool, std::vector<int>& priority);

  void insertDistance(std::vector<std::vector<int> >& dist, std::vector<int>& distOffspringToPool);

  void updateDistance(std::vector<std::vector<int> >& dist, std::vector<int>& distOffspringToPool, int elimIndv);

  void printPoolDistance(std::vector<MMTPartialColoring>& pool, bool expanded = false);

  void printPoolFitness(std::vector<MMTPartialColoring>& pool);

  // stable sets of each newly best partial coloring is added to columns
  void addStableSets(MMTPartialColoring* new_best);
};

#endif
