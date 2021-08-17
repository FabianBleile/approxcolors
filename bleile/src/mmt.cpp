#include "bleile/header/mmt.h"

MMT::MMT(Graph * graph, int L, int T, int timeLimit, int PS, bool set_bounds)
 : graph(graph), L(L), T(T), timeLimit(timeLimit), PS(PS),
 best_col(EvolPartialCol(graph->n, graph)), N(graph->n)

 {
   assert(PS >= 3);
   std::cout << "N = " << N << '\n';

   logger.UB = N+1;
   update_limit = 1000000 / N;
   delta_L = L * graph->dens;
   delta_PS = 1;

   if (set_bounds) adoptBounds();
}

void MMT::start(){
  clock_t t = clock();
  evolInit();
  evolOptimize();
  if(logger.UB != logger.LB) {
    // PHASE2_ColumnOptimization();
  }
  if (best_col.checkColoring()) {
    best_col.toString();
  }

  logger.totTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
}

void MMT::evolInit(){
  EvolPartialCol init = EvolPartialCol(logger.UB, &graph);
  if (init.dsatur() || init.tabuSearch(L, T)) {
    best_col = init;
    logger.UB = init.getNumColors();
  };
  std::cout << "Init UB to " << logger.UB << "\n";
}

void MMT::evolOptimize() {
  // document results for pop effect
  kLogData initKLogData = {logger.UB, 0, INIT_DSATUR};
  logger.kLogData.push_back({logger.UB, 0, INIT_DSATUR});
  // pool : stores the current partial colorings
  //        init default pool with empty partial solutions
  std::vector<EvolPartialCol> pool;

  clock_t T = clock();

  while (logger.UB > logger.LB) {
    clock_t t = clock();
    MMT::status res = evolDecision(logger.UB-1, pool);
    switch (res) {
      case EA_TIME_OUT:
        return;
      default:
        logger.status = res;
        logger.lastItTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
        std::cout << "time : " << (float) (clock() - T)/CLOCKS_PER_SEC << '\n';
    }
    logger.UB = best_col.getNumColors();
    std::cout << "UB updated to " << logger.UB << " with status " << logger.status << "\n";

    // document results for pop effect
    logger.kLogData.push_back({logger.UB, (float) (clock() - T)/CLOCKS_PER_SEC, logger.status});
  }
}

MMT::status MMT::evolDecision(int k, std::vector<EvolPartialCol>& pool) {

  // priority vector : for every vertex keeps track of the total number this
  //                   vertex is left uncolored in the current pool
  std::vector<int> priority(graph.n, PS);

  if (pool.size() < PS) {
    MMT::status res = initPool(k, pool, priority, false);
    if (res != UNSOLVED) {
      return res;
    }
    res = initPool(k, pool, priority, true);
    if (res != UNSOLVED) {
      return res;
    }
  } else {
    std::vector<EvolPartialCol> new_pool;
    for (auto& indv : pool) {
      indv.setK(k);
      if(indv.tabuSearch(L, T)) {
        best_col = indv;
        logger.lastItNumOffsprings = 0;
        return EA;
      }
      insertPool(indv, new_pool, priority);
    }
    pool = new_pool;
  }

  // reset lastItNumOffsprings
  size_t currentItNumOffsprings = 0;

  clock_t t = clock();
  R = N/10;
  int delta_R = R/5;
  std::cout << "R = " << R << '\n';
  float pGreedy = 0.1;

  while (((float) clock() - t)/CLOCKS_PER_SEC < timeLimit) {
    int poolDensityCounter = 0;

    for (size_t iter = 0; iter < update_limit; iter++) {

      EvolPartialCol offspring(k, &graph);
      std::vector<int> distOffspringToPool(PS, 0);

      int parent_1 = (int) rand() % PS;
      int parent_2 = (int) rand() % PS;
      while (parent_2 == parent_1) {
        parent_2 = (int) rand() % PS;
      }

      // generate offspring and if it is not already a solution improve by calling tabuSearch on it
      offspring = EvolPartialCol(k, &graph);
      if(offspring.crossover(pool[parent_1], pool[parent_2]) || offspring.tabuSearch(L, T)) {
        best_col = offspring;
        logger.lastItNumOffsprings = currentItNumOffsprings;
        return EA;
      }

      std::vector<int> sameIndvs;
      std::vector<int> closeIndvs;
      std::vector<int> nearIndvs;
      for (size_t i = 0; i < PS; i++) {
        int tempDist = offspring.distanceTo(&pool[i], true);
        if (tempDist < R) {
          nearIndvs.push_back(i);
          if (tempDist <= R/10) {
            closeIndvs.push_back(i);
            if (tempDist == 0) {
              sameIndvs.push_back(i);
            }
          }
        }
      }

      if (!closeIndvs.empty() || nearIndvs.size() > 3) {
        poolDensityCounter++;
      }

      if (sameIndvs.size() > 1) {
        offspring = EvolPartialCol(k, &graph);
        if(offspring.greedy() || offspring.tabuSearch(10*L, T)) {
          best_col = offspring;
          logger.lastItNumOffsprings = currentItNumOffsprings;
          return EA;
        }
        updatePool(offspring, &pool[sameIndvs[0]], pool, priority);
      } else if ((float) rand()/RAND_MAX < pGreedy && !closeIndvs.empty()) {
        int worst_col = 0;
        int worst_col_fitness;
        for (size_t i = 1; i < closeIndvs.size(); i++) {
          int temp_fitness = pool[closeIndvs[i]].evaluate();
          if (temp_fitness > pool[worst_col].evaluate()) {
            worst_col = closeIndvs[i];
            worst_col_fitness = temp_fitness;
          }
        }
        if (offspring.evaluate() < worst_col_fitness) {
          updatePool(offspring, &pool[worst_col], pool, priority);
        }
      } else if((float) rand()/RAND_MAX < pGreedy && nearIndvs.size() > 3) {
        // there is a near individual in the pool
        int temp_L = (float) poolDensityCounter > (float) 0.15*iter ?
          (int) 10*((float) poolDensityCounter/(iter+1))*L : L;

        // drop offspring and generate new partial coloring with priorityGreedy()
        offspring = EvolPartialCol(k, &graph);

        int res = (float) rand()/RAND_MAX < 0.3 ? offspring.priorityGreedy(priority) : offspring.greedy();

        if(res || offspring.tabuSearch(temp_L, T)) {
          best_col = offspring;
          logger.lastItNumOffsprings = currentItNumOffsprings;
          return EA;
        }

        int worst_col = 0;
        auto tempIndv = closeIndvs.empty() ? nearIndvs : closeIndvs;
        for (size_t i = 1; i < tempIndv.size(); i++) {
          if (pool[tempIndv[i]].evaluate() > pool[worst_col].evaluate()) {
            worst_col = tempIndv[i];
          }
        }
        updatePool(offspring, &pool[worst_col], pool, priority);
      } else {
        // delete worst parent and insert child to pool
        if (pool[parent_1].evaluate() <= pool[parent_2].evaluate()) {
          updatePool(offspring, &pool[parent_2], pool, priority);
        } else {
          updatePool(offspring, &pool[parent_1], pool, priority);
        }
      }

      if (poolDensityCounter > update_limit/5) {
        break;
      }

      logger.totNumOffsprings++;
      currentItNumOffsprings++;
    }
    std::cout << "update_limit = " << update_limit << " poolDensityCounter = " << poolDensityCounter << '\n';
    std::cout << "PS = " << PS << " N = " << N << '\n';
    if (poolDensityCounter <= update_limit/100 && PS < N/20) {
      for (size_t i = 0; i < delta_PS; i++) {
        EvolPartialCol newIndv = EvolPartialCol(k, &graph);

        if(newIndv.priorityGreedy(priority) || newIndv.tabuSearch(L,T)) {
          best_col = newIndv;
          logger.lastItNumOffsprings = currentItNumOffsprings;
          return INIT_DSATUR;
        }
        insertPool(newIndv, pool, priority);
        PS += delta_PS;
      }
    } else if (poolDensityCounter > update_limit/5) {
      R -= delta_R;
    }
    if (L < 250*N) {
      L += delta_L;
    }

    std::cout << "L = " << L << "; "
              << "PS = " << PS << "; "
              << "R = " << R
              << '\n';
  }
  return EA_TIME_OUT;
}

EvolPartialCol* MMT::getColoring(){
  return &best_col;
}

std::stringstream MMT::streamLogs(){
  std::stringstream logs;
  logs << logger.status << ',';
  logs << logger.UB << ',' << logger.LB << ',';
  logs << PS << ',' << L << ',';
  logs << logger.totTimeInSec << ',' << logger.lastItTimeInSec << ',';
  logs << logger.totNumOffsprings << ',' << logger.lastItNumOffsprings << '\n';
  return logs;
}

MMT::status MMT::initPool(int k, std::vector<EvolPartialCol>& pool, std::vector<int>& priority, bool diversify) {

  // clear pool from old inits
  pool.clear();

  // apply (different) initialization algorithms on the pool
  if (!diversify) {
    // DSATUR Block
    for (size_t i = 0; i < PS; i++) {
      EvolPartialCol dsatur_col = EvolPartialCol(k, &graph);
      if(dsatur_col.dsatur() || dsatur_col.tabuSearch(L,T)) {
        best_col = dsatur_col;
        logger.lastItNumOffsprings = 0;
        return MMT::INIT_DSATUR;
      }
      insertPool(dsatur_col, pool, priority);
    }
  } else {
    // GREEDY Block
    for (size_t i = 0; i < PS/3; i++) {
      EvolPartialCol greedy_col = EvolPartialCol(k, &graph);
      if(greedy_col.greedy() || greedy_col.tabuSearch(L,T)) {
        best_col = greedy_col;
        logger.lastItNumOffsprings = 0;
        return MMT::INIT_GREEDY;
      }
      insertPool(greedy_col, pool, priority);
    }
    // DSATUR Block
    for (size_t i = 0; i < PS/3; i++) {
      EvolPartialCol dsatur_col = EvolPartialCol(k, &graph);
      if(dsatur_col.dsatur() || dsatur_col.tabuSearch(L,T)) {
        best_col = dsatur_col;
        logger.lastItNumOffsprings = 0;
        return MMT::INIT_DSATUR;
      }
      insertPool(dsatur_col, pool, priority);
    }
    // TABU SEARCH Block
    int size_tabu_block = PS - pool.size();
    for (size_t i = 0; i < size_tabu_block; i++) {
      EvolPartialCol tabu_col = EvolPartialCol(k, &graph);
      if(tabu_col.dsatur() || tabu_col.tabuSearch(L,T)) {
        best_col = tabu_col;
        logger.lastItNumOffsprings = 0;
        return MMT::INIT_DSATUR;
      }
      insertPool(tabu_col, pool, priority);
    }
  }

  return UNSOLVED;
}

void MMT::insertPool(EvolPartialCol& new_individual, std::vector<EvolPartialCol>& pool, std::vector<int>& priority){
  // update priority
  for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

  // update pool
  pool.push_back(new_individual);
}

void MMT::updatePool(EvolPartialCol& new_individual, EvolPartialCol* old_individual, std::vector<EvolPartialCol>& pool, std::vector<int>& priority){
  // update
  for (const auto & uncol_v : old_individual->uncolored) priority[uncol_v]--;
  for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

  // update pool
  *old_individual = new_individual;
}

void MMT::removePool(int indexRemoveIndv, std::vector<EvolPartialCol>& pool, std::vector<int>& priority){
  // update
  for (const auto & uncol_v : pool[indexRemoveIndv].uncolored) priority[uncol_v]--;

  // remove from pool
  pool.erase(pool.begin() + indexRemoveIndv);
  PS--;
}

std::vector<int> MMT::getWorstIndvs(std::vector<EvolPartialCol>& pool, int returnSize){
  std::vector<int> sumDists(pool.size(),0);
  for (size_t i = 0; i < pool.size(); i++) {
    for (size_t j = i+1; j < pool.size(); j++) {
      int dist_ij = pool[i].distanceTo(&pool[j],true);
      sumDists[i] += dist_ij;
      sumDists[j] += dist_ij;
    }
  }
  vector<int> data_copy = sumDists;
  std::nth_element(data_copy.begin(), data_copy.begin() + returnSize, data_copy.end());
  vector<int> res;
  int kth_element = data_copy[returnSize - 1];
  for (int i = 0; i < sumDists.size(); i++)
    if (sumDists[i] <= kth_element)
      res.push_back(i);
  return res;
}

void MMT::adoptBounds(){
  std::ifstream bounds_file("bleile/bounds.txt");
  std::string s;
  while(getline(bounds_file, s)) {
    std::stringstream ss(s);
    std::string inst;
    int lb, ub, lazylb;
    ss >> inst >> lb >> ub >> lazylb;
    if (this->graph.instance == inst) {
      // logger.UB = ub;
      // logger.LB = lb;
      logger.LB = lazylb;
      std::cout << "Adopt bounds from file bounds.txt: LB = " << lb << ", UB = " << ub << '\n';
      break;
    }
  }
  bounds_file.close();
}

void MMT::printPoolDistance(std::vector<EvolPartialCol>& pool, bool expanded){
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

void MMT::printPoolFitness(std::vector<EvolPartialCol>& pool){
  int sum = 0;
  int best = std::numeric_limits<int>::max();
  for (auto& individual : pool) {
    int temp = individual.evaluate();
    //std::cout << temp << " ; " << individual.uncolored.size() << "\t";
    sum += temp;
    best = temp < best ? temp : best;
  }
  std::cout << "|\tbest = " << best << "; average = " << sum / pool.size() << '\n';
}

int COLORbleile(int ncount, int ecount, int *elist) {
  std::cout << "hello from c++" << '\n';

  return 0;

  // Graph g(ncount, ecount, elist);
  //
  // MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 28000, /*pool size*/ 10, true); // Graph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5
  //
  // mmt.start();
  //
  // return mmt.getColoring()->getNumColors();
}
