#include "header/mmt.h"

void documentation(std::string instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ',' << i << '/' << imax <<',';
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

int main(int argc, char **av) {

  MMTGraph g(argc, av);

  /*
  MMTPartialColoring A(5, &g, 50, 5);
  A.greedy();
  A.toString(6);
  MMTPartialColoring B(5, &g, 50, 5);
  B.greedy();
  B.toString(6);

  std::cout << A.distanceTo(&B, true) << '\n';
  */

  MMT mmt(&g, /*L*/ 20000,/*T*/ 45, /*time limit*/ 120, /*pool size*/ 10, true); // MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

  mmt.start();

  documentation(g.instance, &mmt, 1, 4);

  return 0;
};
