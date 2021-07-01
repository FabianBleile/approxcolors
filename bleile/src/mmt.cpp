#include "bleile/header/mmt.h"

MMT::MMT(MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size, double pGreedy)
 : graph(graph), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3),
 cur_best_coloring(MMTPartialColoring(graph->n, graph, L, T)), pGreedy(pGreedy), N(graph->n)

 {
  // compute lower bound
  // compute upper bound
  logger.UB = N+1;
  measure_best_solution = N*N;

  std::cout << "N = " << N << '\n';
}

void MMT::start(){
  clock_t t = clock();
  PHASE0_EAInit();
  PHASE1_EAOptimizer();
  if(logger.UB != logger.LB) {
    // PHASE2_ColumnOptimization();
  } else {
    cur_best_coloring.greedy();
  }
  cur_best_coloring.checkColoring();
  cur_best_coloring.toString();
  logger.totTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
}

void MMT::PHASE0_EAInit(){
  MMTPartialColoring init = MMTPartialColoring(logger.UB, &graph, L, T);
  init.dsatur();
  init.tabuSearch();
  cur_best_coloring = init;
  logger.UB = init.getNumColors();
}

void MMT::PHASE1_EAOptimizer() {
  while (logger.UB > logger.LB) {
    EADecision(logger.UB-1);
    if (logger.status == EA_TIME_OUT) break;
    logger.UB = cur_best_coloring.getNumColors();
    std::cout << "UB updated to " << logger.UB << " with status " << logger.status << "\n";
  }
}

void MMT::EADecision(int k) {
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
    if(dsatur.greedy() || dsatur.tabuSearch()) {
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
    if (!(iter % 2000)) {
      std::cout << "Anzahl Iterationen " << iter << '\t';
      //printPoolFitness(pool);
      printPoolDistance(pool, false);
    }

    auto parent_1 = std::next(std::begin(pool), (int) rand() % pool.size());
    auto parent_2 = std::next(std::begin(pool), (int) rand() % pool.size());
    while (parent_2 == parent_1) {
      parent_2 = std::next(std::begin(pool), (int) rand() % pool.size());
    }

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
    //std::cout << offspring.evaluate() << " : " << parent_1->evaluate() << " : " << parent_2->evaluate() << " ";
    if (parent_1->evaluate() <= parent_2->evaluate()) {
      //std::cout << "(2) | \t |";
      updatePool(offspring, &(*parent_2), pool, poolSimilarity, priority);
    } else {
      //std::cout << "(1) | \t |";
      updatePool(offspring, &(*parent_1), pool, poolSimilarity, priority);
    }
    iter++;

    logger.totNumOffsprings++;
    logger.lastItNumOffsprings++;
  }
  logger.status = EA_TIME_OUT;
  return;
}

void MMT::PHASE2_ColumnOptimization(){
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

MMTPartialColoring* MMT::getColoring(){
  return &cur_best_coloring;
}

std::stringstream MMT::streamLogs(){
  std::stringstream logs;
  logs << logger.status << ',';
  logs << logger.totTimeInSec << ',' << logger.lastItTimeInSec << ',';
  logs << logger.totNumOffsprings << ',' << logger.lastItNumOffsprings << ',';
  logs << logger.UB << ',' << logger.LB << ',';
  logs << logger.colOpt;
  return logs;
}

void MMT::insertPool(MMTPartialColoring& new_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
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

void MMT::updatePool(MMTPartialColoring& new_individual, MMTPartialColoring* old_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
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

void MMT::printPoolDistance(std::vector<MMTPartialColoring>& pool, bool expanded){
  assert(pool.size() != 0);
  std::cout << "pool distances : " << '\n';
  int sum = 0, size = pool.size();
  for (int i = 0; i < size; i++) {
    for (int j = i; j < size; j++) {
      int approx = pool[i].distanceTo(&pool[j], false);
      int exact = pool[i].distanceTo(&pool[j],true);
      if (expanded) {
        //std::cout << approx << " | " << exact << '\t';
        std::cout << exact << "|\t";
      }
      sum += exact;
    }
    if (expanded) std::cout << '\n';
  }
  std::cout << "avg = " << sum / ((size*(size+1))/2) << '\n';
}

void MMT::printPoolFitness(std::vector<MMTPartialColoring>& pool){
  measure sum = 0;
  measure best = std::numeric_limits<measure>::max();
  for (auto& individual : pool) {
    measure temp = individual.evaluate();
    std::cout << temp << " ; " << individual.uncolored.size() << "\t";
    sum += temp;
    best = temp < best ? temp : best;
  }
  std::cout << "|\tbest = " << best << "; average = " << sum / pool.size() << '\n';

  for (auto& individual : pool) {
    auto indv_colors = individual.colors;
    for (size_t i = 0; i < graph.n; i++) {
      if (indv_colors[i] == 0) {
        std::cout << i << ' ';
      }
    }
    std::cout << '\n';
  }
}

// stable sets of each newly best partial coloring is added to columns
void MMT::addStableSets(MMTPartialColoring* new_best){
  new_best->lockColoring();
  for (auto& stable_set : new_best->color_classes) {
    columns.push(stable_set);
    if (columns.size() > numColOpt) columns.pop();
  }
}
