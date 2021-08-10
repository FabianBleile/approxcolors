#include "bleile/header/c_connector.h"
#include "mmt_connector.h"

int COLORbleile(int ncount, int ecount, int *elist){
  return MMTbleile(ncount, ecount, elist);
}
