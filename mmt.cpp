extern "C" {
  #include "lp.h"
  #include "color_defs.h"
}

#include "mmt_graph.h"
#include "mmt_partial_coloring.h"

#include "color_cluster.cpp"

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

  MMT(int argc, char **av, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5)
   : graph(argc,av), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3),
   cur_best_coloring(MMTPartialColoring(graph.n, &graph, L, T)), pGreedy(pGreedy), N(graph.n)

   {
    // compute lower bound
    // compute upper bound
    logger.UB = N+1;
    measure_best_solution = N*N;

    std::cout << "N = " << N << '\n';
  }

  void start(){
    clock_t t = clock();
    PHASE0_EAInit();
    PHASE1_EAOptimizer();
    if(logger.UB != logger.LB) {
      // PHASE2_ColumnOptimization();
    } else {
      cur_best_coloring.greedy();
    }
    cur_best_coloring.toString();
    logger.totTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
  }

  void PHASE0_EAInit(){
    MMTPartialColoring init = MMTPartialColoring(logger.UB, &graph, L, T);
    init.dsatur();
    init.tabuSearch();
    cur_best_coloring = init;
    logger.UB = init.getNumColors();
  }

  void PHASE1_EAOptimizer() {
    while (logger.UB > logger.LB) {
      EADecision(logger.UB-1);
      if (logger.status == EA_TIME_OUT) break;
      logger.UB = cur_best_coloring.getNumColors();
      std::cout << "UB updated to " << logger.UB << " with status " << logger.status << "\n";
    }
  }

  void EADecision(int k) {
    // reset lastItNumOffsprings
    logger.lastItNumOffsprings = 0;
    // poolSimilarity : collects properties for every individual in the pool
    //                  on which basis two individuals are considered similar
    //                  ( #unclored vertices , fitness )
    std::unordered_set<std::pair<int, measure>, UInt32PairHash> poolSimilarity;

    // priority vector : for every vertex keeps track of the total number this
    //                   vertex is left uncolored in the current pool
    std::vector<int> priority(graph.n, -pool_size);

    // pool : stores the current partial colorings
    //        init default pool with empty partial solutions
    std::vector<MMTPartialColoring> pool;

    // apply different initialization algorithms on the pool
    // 1/3 SEQ , 1/3 DSATUR , 1/3 TABU SEARCH

    // DSATUR Block
    int dsatur_block_size = pool_size/3;
    for (size_t i = 0; i < dsatur_block_size; i++) {
      MMTPartialColoring dsatur = MMTPartialColoring(k, &graph, L, T);
      if(dsatur.dsatur() || dsatur.tabuSearch()) {
        logger.status = INIT_DSATUR;
        cur_best_coloring = dsatur;
        return;
      }
      insertPool(dsatur, pool, poolSimilarity, priority);
    }

    // SEQ Block
    int seq_block_size = pool_size/3;
    for (size_t i = 0; i < seq_block_size; i++) {
      MMTPartialColoring dsatur = MMTPartialColoring(k, &graph, L, T);
      if(dsatur.dsatur() || dsatur.tabuSearch()) {
        logger.status = INIT_GREEDY;
        cur_best_coloring = dsatur;
        return;
      }
      insertPool(dsatur, pool, poolSimilarity, priority);
    }

    // TABU SEARCH Block
    int tabusearch_block_size = pool_size - seq_block_size - dsatur_block_size;
    for (size_t i = 0; i < tabusearch_block_size; i++) {
      MMTPartialColoring tabusearch = MMTPartialColoring(k, &graph, L, T);
      if(tabusearch.tabuSearch()) {
        logger.status = INIT_TABU;
        cur_best_coloring = tabusearch;
        return;
      }
      insertPool(tabusearch, pool, poolSimilarity, priority);
    }

    // printPoolDistance(pool);

    clock_t t = clock();
    size_t iter = 0;
    while (((float) clock() - t)/CLOCKS_PER_SEC < time_limit_sec) {
      auto parent_1 = std::next(std::begin(pool), (int) rand() % pool.size());
      auto parent_2 = std::next(std::begin(pool), (int) rand() % pool.size());

      MMTPartialColoring offspring(k, &graph, L, T);
      // generate offspring and if it is not already a solution improve by calling tabuSearch on it
      if(offspring.crossover(&(*parent_1), &(*parent_2)) || offspring.tabuSearch()) {
        logger.status = EA;
        cur_best_coloring = offspring;
        // printPoolDistance(pool);
        return;
      }

      // std::cout << (offspring->distanceTo(&(*parent_1)) + offspring->distanceTo(&(*parent_2)))/2 << ',';
      /*

        offspring similar to a solution in the pool?
        if yes then trigger priorityGreedy with probability pGreedy
        (offspring is considered similar to another partial coloring if each
        fitness and #uncolered vertices are equal)

      */

      if (poolSimilarity.find(std::make_pair(offspring.uncolored.size(), offspring.evaluate())) != poolSimilarity.end()) {
        // there is a similar individual in the pool
        // srand( (unsigned)time( NULL ) );

        if ((float) rand()/RAND_MAX < pGreedy) {
          // drop offspring and generate new partial coloring with priorityGreedy()
          offspring = MMTPartialColoring(k, &graph, L, T);

          if(offspring.priorityGreedy(priority) || offspring.tabuSearch()) {
            cur_best_coloring = offspring;
            return;
          }
        }
      }

      // delete worst parent and insert child to pool
      if (parent_1->evaluate() <= parent_2->evaluate()) {
        updatePool(offspring, &(*parent_2), pool, poolSimilarity, priority);
      } else {
        updatePool(offspring, &(*parent_1), pool, poolSimilarity, priority);
      }
      iter++;
      if (!(iter % 500)) {
        // std::cout << "Anzahl Iterationen " << iter << '\t'; printPoolFitness(pool);
        // printPoolDistance(pool);
      }

      logger.totNumOffsprings++;
      logger.lastItNumOffsprings++;
    }
    logger.status = EA_TIME_OUT;
    return;
  }

  void PHASE2_ColumnOptimization(){
    // std::cout << "Number of stable sets : " << columns.size() << '\n';
    COLORlp * lp = (COLORlp *) NULL;

    // number of stable sets - max : numColOpt
    int varCount = columns.size();
    // number of rows = number of nodes of G
    const int nodeCount = N;
    // init
    int rval = COLORlp_init (&lp, "colorme");

    // add empty row for every node with row >= 1
    for (int i = 0; i < nodeCount; i++) {
        rval = COLORlp_addrow (lp, 0, (int *) NULL, (double *) NULL, COLORlp_GREATER_EQUAL,
                               1.0, (char*) NULL);
    }

    // add variables/ stable sets
    double coeff[nodeCount];
    std::fill_n(coeff, nodeCount, 1.0);

    while (!columns.empty()) {
      std::unordered_set<nodeid> stable_set = columns.front();
      columns.pop();
      int cind[nodeCount];
      int i = 0;
      for (auto& node : stable_set) {
        cind[i] = node;
        i++;
      }
      // add variable between with objective coeff 1.0 and between 0.0 and 1.0
      rval = COLORlp_addcol (lp, stable_set.size(), cind, coeff,
                    /*obj*/ 1.0,/*lb*/ 0.0,/*ub*/ 1.0,
                    /*sense*/ COLORlp_CONTINUOUS, (char*) NULL);
    }


    // rval = COLORlp_write (lp, "mmt_output.txt");

    int cnt = 0;
    while (COLORlp_get_rowcount (lp) > 0) {
      int rval = COLORlp_optimize(lp);
      // get solution
      double x[varCount];
      // solution available?
      int status = -1;
      COLORlp_get_status(lp, &status);
      if(status != 1) return;
      COLORlp_x (lp, x);
      // retrieve max value
      int argmax = 0;
      for (size_t i = 1; i < varCount; i++)
        if (x[i] > x[argmax])
          argmax = i;

      // argmax-th stable set from lp
      int * colind;
      int colcnt;
      // std::cout << "argmax = " << argmax << " and varCount = " << varCount << '\n';
      rval = COLORlp_get_column(lp, argmax, &colcnt, &colind);

      int temp_colind[colcnt];
      for (size_t i = 0; i < colcnt; i++) {
        temp_colind[i] = colind[i];
      }

      // remove all rows which occur with non zero coeffs in column
      rval = COLORlp_deleterows(lp, colcnt, temp_colind);

      // remove column
      rval = COLORlp_deletecol (/*COLORlp*/ lp, /*first_cind*/ argmax);
      varCount--;

      cnt++;
    }

    std::cout << "NUMBER OF COLORS USED : " << cnt << '\n';
    logger.colOpt = cnt;
    if (cnt < logger.UB) {
      logger.status = COl_OPT;
      // return solution
    }
  }

  MMTPartialColoring* getColoring(){
    return &cur_best_coloring;
  }

  std::stringstream streamLogs(){
    std::stringstream logs;
    logs << logger.status << ',';
    logs << logger.totTimeInSec << ',' << logger.lastItTimeInSec << ',';
    logs << logger.totNumOffsprings << ',' << logger.lastItNumOffsprings << ',';
    logs << logger.UB << ',' << logger.LB << ',';
    logs << logger.colOpt;
    return logs;
  }
  MMTGraph graph;

private:

  enum status {
    UNSOLVED = 1,
    INIT_GREEDY = 2,
    INIT_DSATUR = 4,
    INIT_TABU = 8,
    EA = 16,
    COl_OPT = 32,
    EA_TIME_OUT = 64
  };

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

  struct UInt32PairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const {
        assert(sizeof(std::size_t)>=8);  //Ensure that std::size_t, the type of the hash, is large enough
        //Shift first integer over to make room for the second integer. The two are
        //then packed side by side.
        return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
    }
  };

  int L, T, time_limit_sec, pool_size, measure_best_solution;
  const int N;
  MMTPartialColoring cur_best_coloring;
  double pGreedy;
  const double priority_noise = 0.5;

  LogData logger;

  const int numColOpt = 1000;
  std::queue< std::unordered_set<nodeid> > columns;

  void insertPool(MMTPartialColoring& new_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update priority
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    pool.push_back(new_individual);

    // update columns
    if(new_individual.evaluate() < measure_best_solution){
      measure_best_solution = new_individual.evaluate();
      addStableSets(&new_individual);
    }
  }

  void updatePool(MMTPartialColoring& new_individual, MMTPartialColoring* old_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.erase(std::make_pair(old_individual->uncolored.size(), old_individual->evaluate()));
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update
    for (const auto & uncol_v : old_individual->uncolored) priority[uncol_v]--;
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    *old_individual = new_individual;

    // update columns
    if(new_individual.evaluate() < measure_best_solution){
      measure_best_solution = new_individual.evaluate();
      addStableSets(&new_individual);
    }
  }

  void printPoolDistance(std::vector<MMTPartialColoring>& pool, bool expanded = false){
    assert(pool.size() != 0);
    std::cout << "pool distances : " << '\t';
    int sum = 0, size = pool.size();
    for (int i = 0; i < size; i++) {
      for (int j = i+1; j < size; j++) {
        int approx = pool[i].distanceTo(&pool[j], false);
        int exact = pool[i].distanceTo(&pool[j],true);
        if (expanded) {
          std::cout << approx << " | " << exact << '\t';
        } else {
          sum += exact;
        }
      }
      if (expanded) std::cout << '\n';
    }
    std::cout << "avg = " << sum / ((size*(size+1))/2) << '\n';
  }

  void printPoolFitness(std::vector<MMTPartialColoring>& pool){
    assert(pool.size() != 0);
    measure sum = 0;
    measure best = std::numeric_limits<measure>::max();
    for (auto& individual : pool) {
      measure temp = individual.evaluate();
      sum += temp;
      best = std::min(best, temp);
    }
    std::cout << "best = " << best << "; average = " << sum / pool.size() << '\n';
  }

  // stable sets of each newly best partial coloring is added to columns
  void addStableSets(MMTPartialColoring* new_best){
    new_best->lockColoring();
    for (auto& stable_set : new_best->color_classes) {
      columns.push(stable_set);
      if (columns.size() > numColOpt) columns.pop();
    }
  }
};

void documentation(char *instance, MMT* mmt){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ',';
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

int main(int argc, char **av) {

  MMT mmt(argc, av, /*L*/ 1000,/*T*/ 25, /*time limit*/ 100, /*pool size*/ 20, /*pgreedy*/0.2);

  mmt.start();

  documentation(av[1], &mmt);


  // MMTGraph graph(argc,av);
  //
  // int k = 240;
  //
  // MMTPartialColoring test(k, &graph, 750, 45);
  // test.greedy();
  // test.toString();
  //
  // std::cout << "TEST" << '\n';
  // // std::vector<int> hist0(graph.n, 0);
  // std::vector<int> hist1(graph.n, 0);
  // for (size_t i = 0; i < 100000; i++) {
  //   PartialColoring temp0(k, &graph);
  //   temp0.greedy();
  //   PartialColoring temp1(k, &graph);
  //   temp1.greedy();
  //   // hist0[temp0.distanceTo(&temp1, false)]++;
  //   hist1[temp0.distanceTo(&temp1, true)]++;
  // }
  //
  // for (size_t i = 0; i < graph.n; i++) {
  //   if (hist1[i] > 0) {
  //     std::cout << i << ',' << hist1[i] << ' ';
  //   }
  // }
  // std::cout << '\n';

  // PartialColoringCluster A(10*graph.n, graph.n, k, &graph);
  //
  // for (size_t i = 0; i < 10000; i++) {
  //   PartialColoring temp0(k, &graph);
  //   temp0.greedy();
  //   A.feed(temp0);
  // }
  //
  // std::cout << "FEEDED " << '\n';
  //
  // for (size_t i = 0; i < 10000; i++) {
  //   PartialColoring temp0(k, &graph);
  //   temp0.greedy();
  //   A.test(temp0);
  // }
  //
  // std::cout << "TESTED" << '\n';
  //
  // A.writeToFile();
  //
  // A.toString();

  return 0;
};
