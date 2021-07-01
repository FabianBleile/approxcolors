#include "bleile/header/mmt.h"

MMT::MMT(MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size, double pGreedy)
 : graph(graph), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3),
 cur_best_coloring(MMTPartialColoring(graph->n, graph, L, T)), pGreedy(pGreedy), N(graph->n)

 {
   assert(pool_size >= 3);
   // compute lower bound
   // compute upper bound
   logger.UB = N+1;
   measure_best_solution = N*N;
   R = N / 10 + 1; // pool spacing

   std::cout << "N = " << N << '\n';
}

void MMT::start(){
  clock_t t = clock();
  PHASE0_EAInit();
  PHASE1_EAOptimizer();
  if(logger.UB != logger.LB) {
    // PHASE2_ColumnOptimization();
  }
  if (cur_best_coloring.checkColoring()) {
    cur_best_coloring.toString();
  }

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
    clock_t t = clock();
    MMT::status res = EADecision(logger.UB-1);
    switch (res) {
      case EA_TIME_OUT:
        return;
      default:
        logger.status = res;
        logger.lastItTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
    }
    logger.UB = cur_best_coloring.getNumColors();
    std::cout << "UB updated to " << logger.UB << " with status " << logger.status << "\n";
  }
}

MMT::status MMT::EADecision(int k) {
  // reset lastItNumOffsprings
  logger.lastItNumOffsprings = 0;

  // priority vector : for every vertex keeps track of the total number this
  //                   vertex is left uncolored in the current pool
  std::vector<int> priority(graph.n, pool_size);

  // pool : stores the current partial colorings
  //        init default pool with empty partial solutions
  std::vector<MMTPartialColoring> pool;

  // dist : distance matrix for each indv in the pool
  std::vector<std::vector<int> > dist(pool_size, std::vector<int>(pool_size, 0));

  // apply different initialization algorithms on the pool
  // 1/3 SEQ , 1/3 DSATUR , 1/3 TABU SEARCH

  // DSATUR Block
  int dsatur_block_size = pool_size/3;
  for (size_t i = 0; i < dsatur_block_size; i++) {
    MMTPartialColoring dsatur = MMTPartialColoring(k, &graph, L, T);
    if(dsatur.dsatur() || dsatur.tabuSearch()) {
      cur_best_coloring = dsatur;
      return INIT_DSATUR;
    }
    insertPool(dsatur, pool, priority);
  }

  // SEQ Block
  int seq_block_size = pool_size/3;
  for (size_t i = 0; i < seq_block_size; i++) {
    MMTPartialColoring dsatur = MMTPartialColoring(k, &graph, L, T);
    if(dsatur.greedy() || dsatur.tabuSearch()) {
      cur_best_coloring = dsatur;
      return INIT_GREEDY;
    }
    insertPool(dsatur, pool, priority);
  }

  // TABU SEARCH Block
  int tabusearch_block_size = pool_size - seq_block_size - dsatur_block_size;
  for (size_t i = 0; i < tabusearch_block_size; i++) {
    MMTPartialColoring tabusearch = MMTPartialColoring(k, &graph, L, T);
    if(tabusearch.tabuSearch()) {
      cur_best_coloring = tabusearch;
      return INIT_TABU;
    }
    insertPool(tabusearch, pool, priority);
  }

  // calculate set partition distances
  for (size_t i = 0; i < pool_size; i++) {
    for (size_t j = i+1; j < pool_size; j++) {
      int d_ij = pool[i].distanceTo(&pool[j], true);
      dist[i][j] = d_ij;
      dist[j][i] = d_ij;
    }
  }

  clock_t t = clock();
  size_t iter = 0;

  std::cout << "R = " << R << '\n';

  while (((float) clock() - t)/CLOCKS_PER_SEC < time_limit_sec) {
    /*  CROSSOVER  */
    size_t maxRejects = 20;
    MMTPartialColoring offspring(k, &graph, L, T);
    std::vector<int> distOffspringToPool(pool_size, 0);
    bool accept = true;
    bool acceptOffspring = false;
    for (size_t i = 0; true; i++) {
      auto parent_1 = std::next(std::begin(pool), (int) rand() % pool.size());
      auto parent_2 = std::next(std::begin(pool), (int) rand() % pool.size());
      while (parent_2 == parent_1) {
        parent_2 = std::next(std::begin(pool), (int) rand() % pool.size());
      }

      // generate offspring and if it is not already a solution improve by calling tabuSearch on it
      offspring = MMTPartialColoring(k, &graph, L, T);
      if(offspring.crossover(&(*parent_1), &(*parent_2)) || offspring.tabuSearch()) {
        cur_best_coloring = offspring;
        return EA;
      }

      /*
      if (!(iter % 5)) {
        std::cout << offspring.distanceTo(&(*parent_1), true) << ',' << offspring.distanceTo(&(*parent_2), true) << '\n';
      }
      */

      int maxFitness = 0;
      int elimIndv = -1;
      int neighborhoodSize = 0;
      for (size_t i = 0; i < pool_size; i++) {
        int temp = offspring.distanceTo(&pool[i], true);
        distOffspringToPool[i] = temp;
        if (temp < R) {
          neighborhoodSize++;
          if (neighborhoodSize > pool_size/2) {
            goto retry;
          }
        }
      }
      for (size_t i = 0; i < pool_size; i++) {
        if (pool[i].fitness > maxFitness) {
          maxFitness = pool[i].fitness;
          elimIndv = i;
        }
      }
      updatePool(offspring, &pool[elimIndv], pool, priority);
      updateDistance(dist, distOffspringToPool, elimIndv);
      break;

      retry:
      std::cout << ' ';
        // printPoolDistance(pool, false);

      // reject offspring if it is < R away from another indv in pool and fitness is worse
      /*
      accept = true;
      nodeid elimIndv = -1;
      for (size_t i = 0; i < pool_size; i++) {
        distOffspringToPool[i] = offspring.distanceTo(&pool[i], true);
        if (!acceptOffspring && distOffspringToPool[i] < R) {
          accept = false;
          if (offspring.fitness < pool[i].fitness) {
            acceptOffspring = true;
            elimIndv = i;
          }
        }
      }
      if (acceptOffspring) {
        //  DIRECT REPLACEMENT
        // std::cout << "DIRECT REPLACEMENT " << offspring.fitness << ',' << pool[elimIndv].fitness << '\n';

        updatePool(offspring, &pool[elimIndv], pool, priority);
        updateDistance(dist, distOffspringToPool, elimIndv);
        break;
      } else if (accept) {
        break;
      }*/
    }

    /*if (!accept && !acceptOffspring) {

      //  DIFFERENT MUTATION DUE TO TOO MANY REJECTIONS

      // there is a similar individual in the pool
      if ((float) rand()/RAND_MAX < pGreedy) {
        // drop offspring and generate new partial coloring with priorityGreedy()
        offspring = MMTPartialColoring(k, &graph, L, T);

        if(offspring.priorityGreedy(priority) || offspring.tabuSearch()) {
          cur_best_coloring = offspring;
          return EA;
        }
      }
    }

    if (!acceptOffspring) {

      //  STANDARD REPLACEMENT
      // std::cout << "STANDARD REPLACEMENT" << '\n';

      auto temp_pool = pool;
      std::sort(temp_pool.begin(), temp_pool.end());
      measure medianFitness = temp_pool[ (pool_size-1) / 2].fitness;
      int candidate_1 = (int) rand() % pool_size; // protect fittest
      int candidate_2 = -1;
      int minDist = std::numeric_limits<int>::max();
      for (size_t i = 0; i < pool_size; i++) {
        if (i != candidate_1 && dist[candidate_1][i] < minDist && pool[i].fitness <= medianFitness) {
          candidate_2 = i;
          minDist = dist[candidate_1][i];
        }
      }
      if (pool[candidate_1].fitness < pool[candidate_2].fitness) {
        updatePool(offspring, &pool[candidate_2], pool, priority);
        updateDistance(dist, distOffspringToPool, candidate_2);
      } else {
        updatePool(offspring, &pool[candidate_1], pool, priority);
        updateDistance(dist, distOffspringToPool, candidate_1);
      }
    }*/

    if (iter % 100 == 0) {
      printPoolDistance(pool, false);
      printPoolFitness(pool);
    }

    iter++;

    logger.totNumOffsprings++;
    logger.lastItNumOffsprings++;
  }
  return EA_TIME_OUT;
}

MMTPartialColoring* MMT::getColoring(){
  return &cur_best_coloring;
}

std::stringstream MMT::streamLogs(){
  std::stringstream logs;
  logs << logger.status << ',';
  logs << logger.totTimeInSec << ',' << logger.lastItTimeInSec << ',';
  logs << logger.totNumOffsprings << ',' << logger.lastItNumOffsprings<< ',';
  logs << logger.UB; /* << ',' << logger.LB << ',';
  logs << logger.colOpt;*/
  return logs;
}

void MMT::insertPool(MMTPartialColoring& new_individual, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority){
  // update priority
  for (const auto & uncol_v : new_individual.uncolored) {
    priority[uncol_v]--;
  }

  // update pool
  pool.push_back(new_individual);

  // update columns
  if(new_individual.evaluate() < measure_best_solution){
    measure_best_solution = new_individual.evaluate();
    //addStableSets(&new_individual);
  }
}

void MMT::updatePool(MMTPartialColoring& new_individual, MMTPartialColoring* old_individual, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority){
  // update
  for (const auto & uncol_v : old_individual->uncolored) priority[uncol_v]++;
  for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]--;

  // update pool
  *old_individual = new_individual;

  // update columns
  if(new_individual.evaluate() < measure_best_solution){
    measure_best_solution = new_individual.evaluate();
    //addStableSets(&new_individual);
  }
}

void MMT::updateDistance(std::vector<std::vector<int> >& dist, std::vector<int>& distOffspringToPool, int elimIndv){
  distOffspringToPool[elimIndv] = 0;
  dist[elimIndv] = distOffspringToPool;
  for (size_t i = 0; i < dist.size(); i++) {
    dist[i][elimIndv] = distOffspringToPool[i];
  }
}

void MMT::printPoolDistance(std::vector<MMTPartialColoring>& pool, bool expanded){
  assert(pool.size() != 0);
  std::cout << "pool distances : ";
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
    //std::cout << temp << " ; " << individual.uncolored.size() << "\t";
    sum += temp;
    best = temp < best ? temp : best;
  }
  std::cout << "Fitness | best = " << best << "; average = " << sum / pool.size() << '\n';

  /*
  for (auto& individual : pool) {
    auto indv_colors = individual.colors;
    for (size_t i = 0; i < graph.n; i++) {
      if (indv_colors[i] == 0) {
        std::cout << i << ' ';
      }
    }
    std::cout << '\n';
  }
  */
}


/*------------------------------------------------------------

  SECOND PHASE (INACTIVE)
  Working implementation of the 2nd phase of original MMT
  Due to solely non-improving results phase 2 is dropped here
  This emphasizes experimental results already documented in
  the original paper of MMT

---------------------------------------------------------------*/

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


// stable sets of each newly best partial coloring is added to columns
void MMT::addStableSets(MMTPartialColoring* new_best){
  new_best->lockColoring();
  for (auto& stable_set : new_best->color_classes) {
    columns.push(stable_set);
    if (columns.size() > numColOpt) columns.pop();
  }
}
