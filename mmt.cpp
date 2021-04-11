#include "mmt_graph.h"
#include "mmt_partial_coloring.h"

#include <time.h>
#include <vector>
#include <utility>
#include <queue>

class MMT {
public:
  MMT(int argc, char **av, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5) : graph(argc,av), L(L), T(T), time_limit_sec(time_limit_sec), pool_size((pool_size/3)*3), cur_best_coloring(MMTPartialColoring(graph.n, &graph, L, T)), pGreedy(pGreedy) {
    // compute lower bound
    // compute upper bound
    // UB = graph.n + 1;
    UB = 10;
    LB = 9;
  }

  void testing(){
    MMTPartialColoring A = MMTPartialColoring(7, &graph, L, T);
    MMTPartialColoring B = MMTPartialColoring(7, &graph, L, T);
    MMTPartialColoring C = MMTPartialColoring(7, &graph, L, T);

    A.tabuSearch();
    A.toString();
    B.tabuSearch();
    B.toString();
    C.crossover(&A, &B);
    C.lockColoring();
    C.toString();
  }

  void EAOptimizer() {
    clock_t t = clock();
    while (((float) clock() - t)/CLOCKS_PER_SEC < time_limit_sec && UB > LB) {
      cur_best_coloring = EADecision(UB);
      cur_best_coloring.toString();
      UB--;
    }
    std::cout << "solution found ? " << '\t';
    UB == LB ? std::cout << "YAAY in " << ((float) clock() - t)/CLOCKS_PER_SEC << " secs" << '\n' : std::cout << "no we timed out :o" << '\n';
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

          // copy priority vector and add little noise
          std::vector<int> temp_prio = priority;
          for (size_t i = 0; i < temp_prio.size() - 1; i++) {
            if ((float) rand()/RAND_MAX < priority_noise) {
              std::swap(temp_prio.at(i), temp_prio.at(i+1));
            }
          }

          if(offspring.priorityGreedy(temp_prio) || offspring.tabuSearch()) return offspring;
        }
      }

      // delete worst parent and insert child to pool
      if (pool[0].evaluate() <= pool[1].evaluate()) {
        updatePool(offspring, pool[1], pool, poolSimilarity, priority);
      } else {
        updatePool(offspring, pool[0], pool, poolSimilarity, priority);
      }
      iter++;
    }

    return cur_best_coloring;
  }

  MMTPartialColoring* getColoring(){
    return &cur_best_coloring;
  }

private:
  int L, T, time_limit_sec, pool_size;
  MMTGraph graph;
  int UB, LB;
  MMTPartialColoring cur_best_coloring;
  double pGreedy;
  const double priority_noise = 0.5;

  struct UInt32PairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const {
        assert(sizeof(std::size_t)>=8);  //Ensure that std::size_t, the type of the hash, is large enough
        //Shift first integer over to make room for the second integer. The two are
        //then packed side by side.
        return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
    }
  };

  void insertPool(const MMTPartialColoring& new_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update priority
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    pool.push_back(new_individual);
  }

  void updatePool(const MMTPartialColoring& new_individual, MMTPartialColoring& old_individual, std::vector<MMTPartialColoring>& pool, std::unordered_set<std::pair<int, measure>, UInt32PairHash>& poolSimilarity, std::vector<int>& priority){
    // update poolSimilarity
    poolSimilarity.erase(std::make_pair(old_individual.uncolored.size(), old_individual.evaluate()));
    poolSimilarity.insert(std::make_pair(new_individual.uncolored.size(), new_individual.evaluate()));

    // update
    for (const auto & uncol_v : old_individual.uncolored) priority[uncol_v]--;
    for (const auto & uncol_v : new_individual.uncolored) priority[uncol_v]++;

    // update pool
    old_individual = new_individual;
  }
};

int main(int argc, char **av) {

  MMT mmt(argc, av, /*L*/ 1000,/*T*/ 100, /*time limit*/ 20, /*pool size*/ 99, 0.1);

  mmt.EAOptimizer();
  // mmt.testing();

  return 0;
}
