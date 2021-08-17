#include "../header/mmt.h"
#include <iostream>

extern "C" {

  #include "../header/c_connector.h"

  int MMTbleile(int ncount, int ecount, int *elist) {
    return COLORbleile(ncount, ecount, elist);
  }

}
