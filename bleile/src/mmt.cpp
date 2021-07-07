#include "bleile/header/mmt.h"

MMT::MMT(MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size, bool setBounds)
 : graph(graph), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3),
 cur_best_coloring(MMTPartialColoring(graph->n, graph, L, T)), N(graph->n)

 {
   assert(pool_size >= 3);
   // compute lower bound
   // compute upper bound
   logger.UB = N+1;
   R = N / 10 + 1; // pool spacing

   if (setBounds) {
     std::ifstream bounds_file("bleile/bounds.txt");
     std::string s;
     while(getline(bounds_file, s)) {
       std::stringstream ss(s);
       std::string inst;
       int lb, ub;
       ss >> inst >> lb >> ub;
       if (this->graph.instance == inst) {
         logger.UB = ub;
         logger.LB = lb;
         std::cout << "Adopt bounds from file bounds.txt: LB = " << lb << ", UB = " << ub << '\n';
         break;
       }
     }
     bounds_file.close();
   }

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
  if (init.dsatur() || init.tabuSearch()) {
    cur_best_coloring = init;
    logger.UB = init.getNumColors();
  };
  std::cout << "Init UB to " << logger.UB << "\n";
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
  size_t currentItNumOffsprings = 0;

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
      logger.lastItNumOffsprings = currentItNumOffsprings;
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
      logger.lastItNumOffsprings = currentItNumOffsprings;
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
      logger.lastItNumOffsprings = currentItNumOffsprings;
      return INIT_TABU;
    }
    insertPool(tabusearch, pool, priority);
  }

  clock_t t = clock();
  size_t iter = 0;
  float pGreedy = updatePGreedy(pool, R);
  while (((float) clock() - t)/CLOCKS_PER_SEC < time_limit_sec) {
    if (!(iter % 2000)) {
      std::cout << "Anzahl Iterationen " << iter << '\t';
      printPoolFitness(pool);
      printPoolDistance(pool, false);
    }

    MMTPartialColoring offspring(k, &graph, L, T);
    std::vector<int> distOffspringToPool(pool_size, 0);

    int parent_1 = (int) rand() % pool_size;
    int parent_2 = (int) rand() % pool_size;
    while (parent_2 == parent_1) {
      parent_2 = (int) rand() % pool_size;
    }

    // generate offspring and if it is not already a solution improve by calling tabuSearch on it
    offspring = MMTPartialColoring(k, &graph, L, T);
    if(offspring.crossover(&pool[parent_1], &pool[parent_2]) || offspring.tabuSearch()) {
      cur_best_coloring = offspring;
      return;
    }

    std::vector<int> dist(pool_size);
    for (size_t i = 0; i < pool_size; i++) dist[i] = offspring.distanceTo(&pool[i], false);

    if (std::accumulate(dist.begin(), dist.end(), 0)/(2*pool_size) < R) {
      // there is a near individual in the pool
      if ((float) rand()/RAND_MAX < pGreedy) {
        // drop offspring and generate new partial coloring with priorityGreedy()
        offspring = MMTPartialColoring(k, &graph, 10*L, T);

        if(offspring.priorityGreedy(priority) || offspring.tabuSearch()) {
          cur_best_coloring = offspring;
          return;
        }
        std::cout << "pGreedy = " << pGreedy << '\n';
        pGreedy = updatePGreedy(pool, R);
      }
      updatePool(offspring, elimIndv, pool, priority);
      updateDistance(dist, distOffspringToPool, elimIndv);
    } else if (nearIndvs.size() == 1) {
      //std::cout << "1" << ' ';
      elimIndv = nearIndvs[0];
      if (offspring.fitness < pool[elimIndv].fitness) {
        updatePool(offspring, elimIndv, pool, priority);
        updateDistance(dist, distOffspringToPool, elimIndv);
      }
    } else {
      if (pool[parent_1].fitness < pool[parent_2].fitness) {
        elimIndv = parent_2;
      } else {
        elimIndv = parent_1;
      }
      updatePool(offspring, elimIndv, pool, priority);
      updateDistance(dist, distOffspringToPool, elimIndv);
      //std::cout << "0" << ' ';
    }

    // delete worst parent and insert child to pool
    //std::cout << offspring.evaluate() << " : " << parent_1->evaluate() << " : " << parent_2->evaluate() << " ";
    if (parent_1->evaluate() <= parent_2->evaluate()) {
      updatePool(offspring, &(*parent_2), pool, priority);
    } else {
      updatePool(offspring, &(*parent_1), pool, priority);
    }

    iter++;

    logger.totNumOffsprings++;
    currentItNumOffsprings++;
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

  return;
}

void MMT::updatePool(MMTPartialColoring& new_individual, int old_individual, std::vector<MMTPartialColoring>& pool, std::vector<int>& priority){
  // update
  for (const auto & uncol_v : pool[old_individual].uncolored) priority[uncol_v]++;
  for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]--;

  // update pool
  pool[old_individual] = new_individual;

  // update columns
  if(new_individual.evaluate() < measure_best_solution){
    measure_best_solution = new_individual.evaluate();
    //addStableSets(&new_individual);
  }

  return;
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
  return R/(2*avg);
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
}
