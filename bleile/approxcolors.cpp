#include "header/mmt.h"

void documentation(std::string instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ','; // << i << '/' << imax <<',';
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

int main(int argc, char **av) {

  MMTGraph g(argc, av);

  for (size_t i = 0; i < 1; i++) {
    MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 30, /*pool size*/ 10, true); // MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

    mmt.start();

    documentation(g.instance, &mmt, i, 4);
  }

  return 0;
};
