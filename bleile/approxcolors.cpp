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

  MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 10, /*pool size*/ 20, /*pgreedy*/0.2);

  mmt.start();

  documentation(g.instance, &mmt, 1, 4);

  // MMTGraph g(argc, av);
  //
  // for (size_t i = 0; i < 1000000; i++) {
  //   MMTPartialColoring A(5, &g, 1, 0);
  //   A.greedy();
  //   MMTPartialColoring B(5, &g, 1, 0);
  //   B.greedy();
  //
  //   MMTPartialColoring C(5, &g, 1, 0);
  //   C.crossover(&A,&B);
  //
  //   if (A.distanceTo(&B, true) < A.distanceTo(&C, true) || A.distanceTo(&B, true) <  B.distanceTo(&C, true)) {
  //     A.toString();
  //     B.toString();
  //     C.toString();
  //     std::cout << A.distanceTo(&B, true) << '\n';
  //     std::cout << A.distanceTo(&C, true) << '\n';
  //     std::cout << B.distanceTo(&C, true) << '\n';
  //     std::cout << i << '\n';
  //     break;
  //   }
  // }

  // MMTPartialColoring D(5, &g, 1, 0);
  // D.crossover(&B,&C);
  // // D.toString();
  // std::cout << B.distanceTo(&D) << '\n';
  // std::cout << C.distanceTo(&D) << '\n';

  // MMTGraph g(argc, av);
  // PartialColoring pc(7, &g);
  // pc.greedy();
  // pc.toString();
  // pc.setK(5);
  // pc.toString();

  return 0;
};
