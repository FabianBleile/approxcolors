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

  MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 150, /*pool size*/ 10, /*pgreedy*/0.05); // MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

  mmt.start();

  documentation(g.instance, &mmt, 1, 4);

  return 0;
};
