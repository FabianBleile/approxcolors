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

  MMT(Graph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5)
   : graph(graph), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3),
   best_col(EvolPartialCol(graph->n, graph, L, T)), pGreedy(pGreedy), N(graph->n),
   clustering()

   {
     this->graph = *graph;
    // compute lower bound
    // compute upper bound
    logger.UB = graph->n;
    int_best_solution = N*N;

    std::cout << "N = " << N << '\n';
  }

  void start(){
    clock_t t = clock();
    evolInit();
    evolOptimize();
    if(logger.UB != logger.LB) {
      // PHASE2_ColumnOptimization();
    } else {
      best_col.greedy();
    }
    best_col.toString();
    logger.totTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
  }

  void evolInit(){
    EvolPartialCol init = EvolPartialCol(logger.UB, &graph, L, T);
    init.dsatur();
    init.tabuSearch();
    best_col = init;
    if (init.uncolored.empty()) {
      logger.UB = init.getNumColors();
    }
  }

  void evolOptimize() {
    while (logger.UB > logger.LB) {
      evolDecision(logger.UB-1);
      if (logger.status == EA_TIME_OUT) break;
      logger.UB = best_col.getNumColors();
      std::cout << "UB updated to " << logger.UB << " with status " << logger.status << "\n";
    }
  }

  void evolDecision(int k) {
    // bool analyse_pool = clustering.k == k ? true : false;
    // std::cout << "analyse_pool ? " << analyse_pool << '\n';
    // reset lastItNumOffsprings
    logger.lastItNumOffsprings = 0;
    // poolSimilarity : collects properties for every individual in the pool
    //                  on which basis two individuals are considered similar
    //                  ( #unclored vertices , fitness )
    std::unordered_set<std::pair<int, int>, UInt32PairHash> poolSimilarity;

    // priority vector : for every vertex keeps track of the total number this
    //                   vertex is left uncolored in the current pool
    std::vector<int> priority(graph.n, -pool_size);

    // pool : stores the current partial colorings
    //        init default pool with empty partial solutions
    std::vector<EvolPartialCol> pool;

    // apply different initialization algorithms on the pool
    // 1/3 SEQ , 1/3 DSATUR , 1/3 TABU SEARCH

    // DSATUR Block
    int dsatur_block_size = pool_size/3;
    for (size_t i = 0; i < dsatur_block_size; i++) {
      EvolPartialCol dsatur = EvolPartialCol(k, &graph, L, T);
      if(dsatur.dsatur() || dsatur.tabuSearch()) {
        logger.status = INIT_DSATUR;
        best_col = dsatur;
        return;
      }
      insertPool(dsatur, pool, poolSimilarity, priority);
    }

    // SEQ Block
    int seq_block_size = pool_size/3;
    for (size_t i = 0; i < seq_block_size; i++) {
      EvolPartialCol dsatur = EvolPartialCol(k, &graph, L, T);
      if(dsatur.dsatur() || dsatur.tabuSearch()) {
        logger.status = INIT_GREEDY;
        best_col = dsatur;
        return;
      }
      insertPool(dsatur, pool, poolSimilarity, priority);
    }

    // TABU SEARCH Block
    int tabusearch_block_size = pool_size - seq_block_size - dsatur_block_size;
    for (size_t i = 0; i < tabusearch_block_size; i++) {
      EvolPartialCol tabusearch = EvolPartialCol(k, &graph, L, T);
      if(tabusearch.tabuSearch()) {
        logger.status = INIT_TABU;
        best_col = tabusearch;
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

      EvolPartialCol offspring(k, &graph, L, T);
      // generate offspring and if it is not already a solution improve by calling tabuSearch on it
      if(offspring.crossover(&(*parent_1), &(*parent_2)) || offspring.tabuSearch()) {
        logger.status = EA;
        best_col = offspring;
        // printPoolDistance(pool);
        return;
      }
      // if (analyse_pool) {
      //   clustering.test(offspring, true);
      // }

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
          offspring = EvolPartialCol(k, &graph, L, T);

          if(offspring.priorityGreedy(priority) || offspring.tabuSearch()) {
            best_col = offspring;
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
      // if (analyse_pool) {
      //   clustering.test(offspring, false);
      // }
      iter++;
      logger.totNumOffsprings++;
      logger.lastItNumOffsprings++;
    }
    // printPoolDistance(pool);
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

  EvolPartialCol* getColoring(){
    return &best_col;
  }

  void setClustering(PartialColCluster& clustering){
    this->clustering = clustering;
  }

  void streamClustering(const char * filename){
    assert(clustering.k > 0);
    clustering.writeAnalysisToFile(filename);
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
  Graph graph;

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

  int L, T, time_limit_sec, pool_size, int_best_solution;
  const int N;
  EvolPartialCol best_col;
  double pGreedy;
  const double priority_noise = 0.5;
  PartialColCluster clustering;

  LogData logger;

  const int numColOpt = 1000;
  std::queue< std::unordered_set<nodeid> > columns;

  void insertPool(EvolPartialCol& new_individual, std::vector<EvolPartialCol>& pool, std::unordered_set<std::pair<int, int>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update priority
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    pool.push_back(new_individual);

    // update columns
    if(new_individual.evaluate() < int_best_solution){
      int_best_solution = new_individual.evaluate();
      addStableSets(&new_individual);
    }
  }

  void updatePool(EvolPartialCol& new_individual, EvolPartialCol* old_individual, std::vector<EvolPartialCol>& pool, std::unordered_set<std::pair<int, int>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.erase(std::make_pair(old_individual->uncolored.size(), old_individual->evaluate()));
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update
    for (const auto & uncol_v : old_individual->uncolored) priority[uncol_v]--;
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    *old_individual = new_individual;

    // update columns
    if(new_individual.evaluate() < int_best_solution){
      int_best_solution = new_individual.evaluate();
      addStableSets(&new_individual);
    }
  }

  void printPoolDistance(std::vector<EvolPartialCol>& pool, bool expanded = false){
    assert(pool.size() != 0);
    std::cout << "pool distances : " << '\n';
    int sum = 0, size = pool.size();
    for (int i = 0; i < size; i++) {
      for (int j = i+1; j < size; j++) {
        // int approx = pool[i].distanceTo(&pool[j], false);
        int exact = pool[i].distanceTo(&pool[j],true);
        if (expanded) {
          std::cout << exact << '\t'; // << pool[i].uncolored.size() << '\t';
        }
        sum += exact;
      }
      if (expanded) std::cout << '\n';
    }
    std::cout << "avg = " << sum / (size*(size-1)/2) << '\n';
  }

  void printPoolFitness(std::vector<EvolPartialCol>& pool){
    assert(pool.size() != 0);
    int sum = 0;
    int best = std::numeric_limits<int>::max();
    for (auto& individual : pool) {
      int temp = individual.evaluate();
      sum += temp;
      best = std::min(best, temp);
    }
    std::cout << "best = " << best << "; average = " << sum / pool.size() << '\n';
  }

  // stable sets of each newly best partial coloring is added to columns
  void addStableSets(EvolPartialCol* new_best){
    new_best->buildColorClasses();
    for (auto& stable_set : new_best->color_classes) {
      columns.push(stable_set);
      if (columns.size() > numColOpt) columns.pop();
    }
  }
};

void documentation(char *instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ',' << i << '/' << imax << ',';
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

void createCluster(const char * filename, Graph * graph, int k, int num_center, int L = 100, int T = 10){

  PartialColCluster A(num_center, graph->n, k, graph);

  int it  = num_center*10000;

  std::cout << std::endl;

  for (size_t i = 0; i < it; i++) {
    PartialCol basetemp0(k, graph);
    basetemp0.greedy();
    // temp0.tabuSearch();
    // PartialCol basetemp0 = temp0;
    A.feed(basetemp0);
    // EvolPartialCol temp1(k, &graph, L, T);
    // temp1.dsatur();
    // temp1.tabuSearch();
    // PartialCol basetemp1 = temp1;
    // A.feed(basetemp1);
    if((i % (it/500)) == 0){
      cout << "\e[A\r\e[0K"<< ((float) i/it)*100 << '%'<<endl;;
    }
  }

  A.writeToFile(filename);
}

void testCluster(const char * filename, Graph * graph, int k, int num_center, int L = 500, int T = 10){
  PartialColCluster B(filename, graph);

  int it  = num_center*10000;

  std::cout << std::endl;

  for (size_t i = 0; i < it; i++) {
    PartialCol basetemp0(k, graph);
    // basetemp0.greedy();
    // temp0.tabuSearch();
    // PartialCol basetemp0 = temp0;
    B.test(basetemp0, true);
    // EvolPartialCol temp1(k, graph, L, T);
    // temp1.dsatur();
    // temp1.tabuSearch();
    // PartialCol basetemp1 = temp1;
    // B.test(basetemp1, true);
    if((i % (it/500)) == 0){
      cout << "\e[A\r\e[0K"<< ((float) i/it)*100 << '%'<<endl;;
    }
  }

  B.writeAnalysisToFile(filename);
}

int main(int argc, char **av) {

  int PS = 20, L = 5000, T = 45, time_limit = 30;

  // string filename = "DSC250.5.clu";
  // const char * arrfilename = filename.c_str();
  // int k = 28;
  //
  Graph graph(argc,av);
  //
  // EvolPartialCol col(k, &graph, L, T);
  // col.dsatur();
  // col.tabuSearch();

  // createCluster(arrfilename, &graph, 28, 50);
  //
  // testCluster(arrfilename, &graph, k, 50);

  // PartialColCluster B(arrfilename, &graph);

  // B.toString(true);

  for (size_t i = 1; i <= 1; i++) {
    MMT mmt(&graph, L, T, time_limit, PS, 0.1);
  //
  //   mmt.setClustering(B);
  //
    mmt.start();
  //
  //   mmt.streamClustering(arrfilename);
    documentation(av[1], &mmt, i, 4);
  }

  //
  //
  // PartialColCluster A(PS, graph.n, k, &graph);
  //
  // for (size_t i = 0; i < PS*50; i++) {
  //   PartialCol temp0(k, &graph);
  //   temp0.greedy();
  //   A.feed(temp0);
  //   PartialCol temp1(k, &graph);
  //   temp1.dsatur();
  //   A.feed(temp1);
  // }
  //
  // A.writeToFile(arrfilename);
  //
  // A.toString();

  // for (size_t i = 0; i < 10000; i++) {
  //   PartialCol temp0(k, &graph);
  //   temp0.greedy();
  //   A.test(temp0);
  // }
  //
  // std::cout << "TESTED" << '\n';

  return 0;
};
