#include "header/mmt.h"

void documentation(char *instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ',' << i << '/' << imax;
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

int main(int argc, char **av) {

  MMTGraph g(argc, av);

  MMT mmt(&g, /*L*/ 1000,/*T*/ 45, /*time limit*/ 30, /*pool size*/ 5, /*pgreedy*/0.2);

  mmt.start();

  documentation(g.instance, &mmt, 1, 4);

  /*
  MMTPartialColoring D(30, &g, 10000000, 45);
  D.dsatur();
  clock_t t = clock();
  D.tabuSearch();
  std::cout << ((float) clock() - t)/CLOCKS_PER_SEC << '\n';
  D.toString(31);
  */

  return 0;
};
