#include "bleile/header/mmt.h"

MMT::MMT(MMTGraph * graph, int L, int T, int timeLimit, int PS, bool setBounds)
 : graph(graph), L(L), T(T), timeLimit(timeLimit), PS(PS),
 cur_best_coloring(MMTPartialColoring(graph->n, graph, L, T)), N(graph->n)

 {
   assert(PS >= 3);
   std::cout << "N = " << N << '\n';

   logger.UB = N+1;
   R = N/5; // pool spacing
   updateLimit = 1000000 / N;
   deltaL = L * graph->dens;
   deltaPS = std::max((float) PS * graph->dens, (float) 1);

   if (setBounds) {
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
  if (init.dsatur() || init.tabuSearch()) {
    cur_best_coloring = init;
    logger.UB = init.getNumColors();
  };
  std::cout << "Init UB to " << logger.UB << "\n";
}

void MMT::PHASE1_EAOptimizer() {
  // document results for pop effect
  kLogData initKLogData = {logger.UB, 0, INIT_DSATUR};
  logger.kLogData.push_back({logger.UB, 0, INIT_DSATUR});
  // pool : stores the current partial colorings
  //        init default pool with empty partial solutions
  std::vector<MMTPartialColoring> pool;

  clock_t T = clock();

  while (logger.UB > logger.LB) {
    clock_t t = clock();
    MMT::status res = EADecision(logger.UB-1, pool);
    switch (res) {
      case EA_TIME_OUT:
        return;
      default:
        logger.status = res;
        logger.lastItTimeInSec = ((float) clock() - t)/CLOCKS_PER_SEC;
    }
    logger.UB = cur_best_coloring.getNumColors();
    std::cout << "UB updated to " << logger.UB << " with status " << logger.status << "\n";

    // document results for pop effect
    logger.kLogData.push_back({logger.UB, (float) (clock() - T)/CLOCKS_PER_SEC, logger.status});
  }
}

MMT::status MMT::EADecision(int k, std::vector<MMTPartialColoring>& pool) {

  // priority vector : for every vertex keeps track of the total number this
  //                   vertex is left uncolored in the current pool
  std::vector<int> priority(graph.n, PS);

  pool.clear();

  if (pool.size() < PS) {
    MMT::status res = initPool(k, pool, priority);
    if (res != UNSOLVED) {
      return res;
    }
  } else {
    std::vector<MMTPartialColoring> new_pool;
    for (auto& indv : pool) {
      indv.setK(logger.UB-1);
      if(indv.tabuSearch()) {
        cur_best_coloring = indv;
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
  float pGreedy = updatePGreedy(pool, R);

  while (((float) clock() - t)/CLOCKS_PER_SEC < timeLimit) {
    int poolDensityCounter = 0;

    for (size_t iter = 0; iter < updateLimit; iter++) {
      MMTPartialColoring offspring(k, &graph, L, T);
      std::vector<int> distOffspringToPool(PS, 0);

      int parent_1 = (int) rand() % PS;
      int parent_2 = (int) rand() % PS;
      while (parent_2 == parent_1) {
        parent_2 = (int) rand() % PS;
      }

      // generate offspring and if it is not already a solution improve by calling tabuSearch on it
      offspring = MMTPartialColoring(k, &graph, L, T);
      if(offspring.crossover(pool[parent_1], pool[parent_2]) || offspring.tabuSearch()) {
        cur_best_coloring = offspring;
        logger.lastItNumOffsprings = currentItNumOffsprings;
        return EA;
      }

      std::vector<int> nearIndvs;
      for (size_t i = 0; i < PS; i++) {
        if (offspring.distanceTo(&pool[i], false) < R) {
          nearIndvs.push_back(i);
        }
      }

      if (nearIndvs.size() > 2) {
        poolDensityCounter++;
        // there is a near individual in the pool
        if ((float) rand()/RAND_MAX < pGreedy) {
          // drop offspring and generate new partial coloring with priorityGreedy()
          offspring = MMTPartialColoring(k, &graph, 2*L, T);

          if(offspring.priorityGreedy(priority) || offspring.tabuSearch()) {
            cur_best_coloring = offspring;
            logger.lastItNumOffsprings = currentItNumOffsprings;
            return EA;
          }
        }

        int worstNearIndv = 0;
        for (size_t i = 1; i < nearIndvs.size(); i++) {
          if (pool[nearIndvs[i]].evaluate() > pool[worstNearIndv].evaluate()) {
            worstNearIndv = nearIndvs[i];
          }
        }
        updatePool(offspring, &pool[worstNearIndv], pool, priority);
      } else {
        // delete worst parent and insert child to pool
        if (pool[parent_1].evaluate() <= pool[parent_2].evaluate()) {
          updatePool(offspring, &pool[parent_2], pool, priority);
        } else {
          updatePool(offspring, &pool[parent_1], pool, priority);
        }
      }

      logger.totNumOffsprings++;
      currentItNumOffsprings++;
    }
    std::cout << "updateLimit = " << updateLimit << " poolDensityCounter = " << poolDensityCounter << '\n';
    std::cout << "PS = " << PS << " N = " << N << '\n';
    if (poolDensityCounter < updateLimit/5 && PS < N/20) {
      for (size_t i = 0; i < deltaPS; i++) {
        MMTPartialColoring newIndv = MMTPartialColoring(k, &graph, L, T);

        if(newIndv.priorityGreedy(priority) || newIndv.tabuSearch()) {
          cur_best_coloring = newIndv;
          logger.lastItNumOffsprings = currentItNumOffsprings;
          return INIT_DSATUR;
        }
        insertPool(newIndv, pool, priority);
        PS++;
      }
    } else if (poolDensityCounter > 4*updateLimit/5 && PS >=10) {
      int numRemoveIndvs = std::min(PS-5, deltaPS/2);
      if (numRemoveIndvs > 0) {
        std::vector<int> worstIndvs = getWorstIndvs(pool, deltaPS/2);
        std::sort(worstIndvs.begin(), worstIndvs.end());
        for (size_t i = 0; i < worstIndvs.size(); i++) {
          worstIndvs[i] -= i;
        }
        for (auto indx : worstIndvs) removePool(indx,pool,priority);

        PS = pool.size();
      }
    }
    if (L < 250*N) {
      L += deltaL;
    }
    pGreedy = updatePGreedy(pool, R);

    std::cout << "L = " << L << "; "
              << "PS = " << PS
              << '\n';
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
  logs << logger.totNumOffsprings << ',' << logger.lastItNumOffsprings << ',';
  logs << logger.UB << ',' << logger.LB << '\n';
  return logs;
}

MMT::status MMT::initPool(int k, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority) {

  // clear pool from old inits
  pool.clear();

  // apply (different) initialization algorithms on the pool
  // DSATUR Block
  for (size_t i = 0; i < PS; i++) {
    MMTPartialColoring dsatur = MMTPartialColoring(k, &graph, this->L, this->T);
    if(dsatur.dsatur() || dsatur.tabuSearch()) {
      cur_best_coloring = dsatur;
      logger.lastItNumOffsprings = 0;
      return MMT::INIT_DSATUR;
    }
    insertPool(dsatur, pool, priority);
  }

  return UNSOLVED;
}

void MMT::insertPool(MMTPartialColoring& new_individual, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority){
  // update priority
  for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

  // update pool
  pool.push_back(new_individual);
}

void MMT::updatePool(MMTPartialColoring& new_individual, MMTPartialColoring* old_individual, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority){
  // update
  for (const auto & uncol_v : old_individual->uncolored) priority[uncol_v]--;
  for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

  // update pool
  *old_individual = new_individual;
}

void MMT::removePool(int indexRemoveIndv, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority){
  // update
  for (const auto & uncol_v : pool[indexRemoveIndv].uncolored) priority[uncol_v]--;

  // remove from pool
  pool.erase(pool.begin() + indexRemoveIndv);
}

float MMT::updatePGreedy(std::vector<MMTPartialColoring>& pool, int R) {
  int sum = 0, size = pool.size();
  for (int i = 0; i < size; i++) {
    for (int j = i + 1; j < size; j++) {
      int dist;
      dist = pool[i].distanceTo(&pool[j], false);
      sum += dist;
    }
  }
  float avg = sum / (size*(size+1)/2 - size);
  return R/(3*avg);
}

std::vector<int> MMT::getWorstIndvs(std::vector<MMTPartialColoring>& pool, int returnSize){
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
  std::cout << "|\tbest = " << best << "; average = " << sum / pool.size() << '\n';
}

int COLORbleile(int ncount, int ecount, int *elist) {
  std::cout << "hello from c++" << '\n';
  return 0;
}
