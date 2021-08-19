#include "../header/mmt.h"
#include <iostream>

extern "C" {

  #include "../header/c_connector.h"

  int MMTbleile(int ncount, int ecount, int *elist, int *ncolors,
                   COLORset **colorclasses, int L, int T, int time_limit, int lb) {
    return COLORbleile(ncount, ecount, elist, ncolors,
                     colorclasses, L, T, time_limit, lb);
  }

}
