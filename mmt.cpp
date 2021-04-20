extern "C" {
  #include "lp.h"
  #include "color_defs.h"
}

#include "mmt_graph.h"
#include "mmt_partial_coloring.h"

#include <time.h>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <array>

class MMT {
public:
  MMT(int argc, char **av, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5)
   : graph(argc,av), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3),
   cur_best_coloring(MMTPartialColoring(graph.n, &graph, L, T)), pGreedy(pGreedy), N(graph.n)

   {
    // compute lower bound
    // compute upper bound
    UB = 15;
    LB = 13;
    measure_best_solution = N*N;

    std::cout << "N = " << N << '\n';
  }

  void start(){
    PHASE1_EAOptimizer();
    if(UB != LB) {
      PHASE2_ColumnOptimization();
    } else {
      cur_best_coloring.greedy();
    }
    cur_best_coloring.toString();
  }

  void PHASE1_EAOptimizer() {
    clock_t t = clock();
    while (((float) clock() - t)/CLOCKS_PER_SEC < time_limit_sec && UB > LB) {
      cur_best_coloring = EADecision(UB);
      cur_best_coloring.toString();
      if (UB == LB) {
        std::cout << "YAAY in " << ((float) clock() - t)/CLOCKS_PER_SEC << " secs" << '\n';
        return;
      }
      UB--;
    }
  }

  MMTPartialColoring EADecision(int UB) {
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

    // SEQ Block
    int seq_block_size = pool_size/3;
    for (size_t i = 0; i < seq_block_size; i++) {
      MMTPartialColoring seq = MMTPartialColoring(UB-1, &graph, L, T);
      if(seq.greedy() || seq.tabuSearch()) return seq;
      insertPool(seq, pool, poolSimilarity, priority);
    }

    // DSATUR Block
    int dsatur_block_size = pool_size/3;
    for (size_t i = 0; i < dsatur_block_size; i++) {
      MMTPartialColoring dsatur = MMTPartialColoring(UB-1, &graph, L, T);
      if(dsatur.dsatur() || dsatur.tabuSearch()) return dsatur;
      insertPool(dsatur, pool, poolSimilarity, priority);
    }

    // TABU SEARCH Block
    int tabusearch_block_size = pool_size - seq_block_size - dsatur_block_size;
    for (size_t i = 0; i < tabusearch_block_size; i++) {
      MMTPartialColoring tabusearch = MMTPartialColoring(UB-1, &graph, L, T);
      if(tabusearch.tabuSearch()) return tabusearch;
      insertPool(tabusearch, pool, poolSimilarity, priority);
    }

    clock_t t = clock();
    size_t iter = 0;
    while (((float) clock() - t)/CLOCKS_PER_SEC < time_limit_sec) {
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::shuffle(pool.begin(), pool.end(), std::default_random_engine(seed));

      MMTPartialColoring offspring(UB-1, &graph, L, T);
      // generate offspring and if it is not already a solution improve by calling tabuSearch on it
      if(offspring.crossover(&pool[0], &pool[1]) || offspring.tabuSearch()) return offspring;
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
          offspring = MMTPartialColoring(UB-1, &graph, L, T);

          if(offspring.priorityGreedy(priority) || offspring.tabuSearch()) return offspring;
        }
      }

      // delete worst parent and insert child to pool
      if (pool[0].evaluate() <= pool[1].evaluate()) {
        updatePool(offspring, pool[1], pool, poolSimilarity, priority);
      } else {
        updatePool(offspring, pool[0], pool, poolSimilarity, priority);
      }
      iter++;
      if (!(iter % 2000)) {
        std::cout << "Anzahl Iterationen " << iter << " mit durchschnittsfitness " << averageFitness(pool) << '\n';
      }
    }
    return cur_best_coloring;
  }

  void PHASE2_ColumnOptimization(){
    std::cout << "Number of stable sets : " << columns.size() << '\n';
    COLORlp * lp = (COLORlp *) NULL;

    // number of stable sets - max : numColOpt
    const int varCount = columns.size();
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


    rval = COLORlp_write (lp, "mmt_output.txt");

    int cnt = 0;
    while (COLORlp_get_rowcount (lp) > 0) {
      int rval = COLORlp_optimize(lp);
      // get solution
      double x[varCount];
      rval = COLORlp_x (lp, x);
      // retrieve max value
      int argmax = 0;
      for (size_t i = 1; i < varCount; i++)
        if (x[i] > x[argmax])
          argmax = i;

      // argmax-th stable set from lp
      int * colind;
      int colcnt;
      rval = COLORlp_get_column(lp, argmax, &colcnt, &colind);

      std::cout << "argmax = " << argmax << " and colcnt = " << colcnt << '\n';
      int temp_colind[colcnt];
      for (size_t i = 0; i < colcnt; i++) {
        std::cout << colind[i] << '\n';
        temp_colind[i] = colind[i];
      }

      // remove all rows which occur with non zero coeffs in column
      rval = COLORlp_deleterows(lp, colcnt, temp_colind);

      // remove column
      rval = COLORlp_deletecols (/*COLORlp*/ lp, /*first_cind*/ argmax, /*last_cind*/ argmax);

      cnt++;
    }

    std::cout << "NUMBER OF COLORS USED : " << cnt << '\n';
  }

  MMTPartialColoring* getColoring(){
    return &cur_best_coloring;
  }

private:
  int L, T, time_limit_sec, pool_size, measure_best_solution;
  MMTGraph graph;
  const int N;
  int UB, LB;
  MMTPartialColoring cur_best_coloring;
  double pGreedy;
  const double priority_noise = 0.5;

  const int numColOpt = 1000;
  std::queue< std::unordered_set<nodeid> > columns;

  struct UInt32PairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const {
        assert(sizeof(std::size_t)>=8);  //Ensure that std::size_t, the type of the hash, is large enough
        //Shift first integer over to make room for the second integer. The two are
        //then packed side by side.
        return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
    }
  };

  // Hashing tecnique for stable set hashing
  // struct NBitStringHash {
  //   std::size_t operator()(bool[N] nstring) const {
  //       return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
  //   }
  // };

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

  void updatePool(MMTPartialColoring& new_individual, MMTPartialColoring& old_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.erase(std::make_pair(old_individual.uncolored.size(), old_individual.evaluate()));
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update
    for (const auto & uncol_v : old_individual.uncolored) priority[uncol_v]--;
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    old_individual = new_individual;

    // update columns
    if(new_individual.evaluate() < measure_best_solution){
      measure_best_solution = new_individual.evaluate();
      addStableSets(&new_individual);
    }
  }

  measure averageFitness(std::vector<MMTPartialColoring>& pool){
    assert(pool.size() != 0);
    measure sum = 0;
    for (auto& individual : pool) {
      sum += individual.evaluate();
    }
    return sum / pool.size();
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

void testLP(){
  COLORlp * lp = (COLORlp *) NULL;

  const int varCount = 5;
  const int nodeCount = 3;

  int rval = COLORlp_init (&lp, "colorme");

  // add cons
  for (int i = 0; i < nodeCount; i++) {
      rval = COLORlp_addrow (lp, 0, (int *) NULL, (double *) NULL, COLORlp_GREATER_EQUAL,
                             1.0, (char*) NULL);
  }

  double coeff[nodeCount];
  std::fill_n(coeff, nodeCount, 1.0);
  int cind[nodeCount];
  std::iota(std::begin(cind), std::end(cind), 0);

  for (int i = 0; i < varCount; i++) {
    // add variable between with objective coeff 1.0 and between 0.0 and 1.0
    rval = COLORlp_addcol (lp, nodeCount,
                 cind, coeff,/*obj*/ 1.0,/*lb*/ 0.0,/*ub*/ 1.0,
                 /*sense*/ COLORlp_CONTINUOUS, (char*) NULL);
  }

  rval = COLORlp_deletecols (/*COLORlp*/ lp, /*first_cind*/ 2, /*last_cind*/ 2);

  int dellist[2] = {0,2};
  //rval = COLORlp_deleterows(/*COLORlp*/ lp, /*numdel*/ (int) 2, /*dellist*/ (int **) &dellist);

  rval = COLORlp_write (lp, "mmt_output.txt");
}

int main(int argc, char **av) {

  MMT mmt(argc, av, /*L*/ 500,/*T*/ 100, /*time limit*/ 5, /*pool size*/ 100, 0.1);

  mmt.start();
  // mmt.testing();

  // testLP();

  return 0;
}
