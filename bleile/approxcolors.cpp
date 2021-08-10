#include "header/mmt.h"

void documentation(std::string instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ','; // << i << '/' << imax <<',';
  doc << mmt->streamLogs().rdbuf();
  doc.close();
}

int main(int argc, char **av) {

  MMTGraph g(argc, av);

  std::vector<std::vector<MMT::kLogData> > v;
  std::vector<int> v_totTime;

  for (size_t i = 0; i < 1; i++) {
    MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 3000, /*pool size*/ 10, true); // MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

    mmt.start();

    documentation(g.instance, &mmt, i, 4);

    v.push_back(mmt.logger.kLogData);
    v_totTime.push_back(mmt.logger.totTimeInSec);
  }

  for (int i = 0; i < v.size(); i++) {
    std::cout << v[i][v[i].size()-1].k << ',' << v[i][v[i].size()-1].timeStampEnd - v[i][v[i].size()-2].timeStampEnd << ' ' << v_totTime[i] << '\n';

    for (size_t j = 0; j < v[i].size(); j++) {
      std::cout << v[i][j].k << ',' << v[i][j].timeStampEnd << '\n';
    }
  }

  return 0;
};
