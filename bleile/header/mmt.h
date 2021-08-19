#ifndef MMT_H_
#define MMT_H_

#include "mmt_graph.h"
#include "mmt_partial_coloring.h"
#include "../utils/pair_hasher.hpp"

#include <time.h>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>

typedef struct COLORset {
    int count;
    int age;
    int *members;
    struct COLORset *next;
} COLORset;

int COLORbleile(int ncount, int ecount, int *elist, int *ncolors,
                 COLORset **colorclasses, int L, int T, int time_limit, int lb);

class MMT {
public:

  Graph graph;

  enum status {
    UNSOLVED = 1,
    INIT_GREEDY = 2,
    INIT_DSATUR = 4,
    INIT_TABU = 8,
    EA = 16,
    COl_OPT = 32,
    EA_TIME_OUT = 64
  };

  struct kLogData {
    int k;
    float timeStampEnd;
    status status;
  };

  struct LogData {
    status status = UNSOLVED;
    int totTimeInSec = 0;
    int lastItTimeInSec = 0;
    int totNumOffsprings = 0;
    int lastItNumOffsprings = 0;
    int UB;
    int LB = 2;
    std::vector<kLogData> kLogData;
  };

  MMT(Graph * graph, int L, int T, int time_limit_sec, int pool_size = 10, bool set_bounds = false, int lb = 2);

  void start();

  void evolInit();

  void evolOptimize();

  status evolDecision(int k, std::vector<EvolPartialCol>& pool, clock_t T);

  EvolPartialCol* getColoring();

  std::stringstream streamLogs();

  LogData logger;

private:

  int L, T, R, timeLimit, PS;
  int update_limit, delta_L, delta_PS, delta_R;
  const int N;
  EvolPartialCol best_col;

  std::queue< std::unordered_set<nodeid> > columns;

  status initPool(int k, std::vector<EvolPartialCol>& pool, std::vector<int>& priority, bool diversify);

  void insertPool(EvolPartialCol& new_individual, std::vector<EvolPartialCol>& pool, std::vector<int>& priority);

  void updatePool(EvolPartialCol& new_individual, EvolPartialCol* old_individual,std::vector<EvolPartialCol>& pool, std::vector<int>& priority);

  void removePool(int indexRemoveIndv, std::vector<EvolPartialCol>& pool, std::vector<int>& priority);

  std::vector<int> getWorstIndvs(std::vector<EvolPartialCol>& pool, int returnSize);

  void adoptBounds();

  void printPoolDistance(std::vector<EvolPartialCol>& pool, bool expanded = false);

  void printPoolFitness(std::vector<EvolPartialCol>& pool);
};

#endif
