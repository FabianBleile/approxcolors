#include "mmt_graph.h"
#include "mmt_partial_coloring.h"

int main(int argc, char **av) {

  MMTGraph g(argc, av);

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
