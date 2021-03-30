extern "C" {
  #include "mmt_read.h"
}

#include "mmt_graph.h"
#include "mmt_partial_coloring.h"

int main(int argc, char **av) {
  int n = 0, m = 0;
  int *elist;
  read_graph(argc, av, &n, &m, &elist);

  MMTGraph g(n, m, &elist);

  g.toString();

  int GOAL = 6;

  MMTPartialColoring tbs(GOAL, &g, 3000, 100);
  tbs.tabuSearch();
  tbs.toString(GOAL+1);

  // MMTPartialColoring seq(GOAL, &g, 3000, 100);
  // seq.greedy();
  // seq.toString(GOAL+1);

  MMTPartialColoring dsatur(GOAL, &g, 3000, 100);
  dsatur.dsatur();
  dsatur.toString(GOAL+1);

  MMTPartialColoring crossover(GOAL, &g, 3000, 100);
  crossover.crossover(&tbs, &crossover);
  crossover.toString(GOAL+1);

  return 0;
}
