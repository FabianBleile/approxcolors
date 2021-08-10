#include "../header/mmt.h"
#include "../header/mmt_graph.h"
#include "../header/c_connector.h"
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif


int MMTbleile(int ncount, int ecount, int *elist) {
  MMTGraph g(ncount, ecount, elist);

  MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 28000, /*pool size*/ 10, true); // MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

  mmt.start();

  return 0;
}

#ifdef __cplusplus
}
#endif
