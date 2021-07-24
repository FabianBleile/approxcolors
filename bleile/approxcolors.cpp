#include "header/mmt.h"

void documentation(std::string instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ','; // << i << '/' << imax <<',';
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

void cliqueAuswertung(MMTGraph& graph){
  std::ifstream bounds_file("bleile/bounds.txt");
  std::string s;
  int LB;
  while(getline(bounds_file, s)) {
    std::stringstream ss(s);
    std::string inst;
    int lb, ub, lazylb;
    ss >> inst >> lb >> ub >> lazylb;
    if (graph.instance == inst) {
      LB = lb;
      break;
    }
  }
  bounds_file.close();
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << graph.instance << ' '
            << graph.n << ' '
            << (*graph.getClique()).size() << ' '
            << (float) (*graph.getClique()).size()/graph.n << ' '
            << (float) (*graph.getClique()).size()/LB
            << '\n';
  doc.close();
}

int main(int argc, char **av) {

  MMTGraph g(argc, av);

  // cliqueAuswertung(g);

  std::vector<std::vector<MMT::kLogData> > v;

  for (size_t i = 0; i < 10; i++) {
    MMT mmt(&g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 120, /*pool size*/ 10, true); // MMTGraph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

    mmt.start();

    documentation(g.instance, &mmt, i, 4);

    v.push_back(mmt.logger.kLogData);
  }

  for (auto i : v) {
    std::cout << i[i.size()-1].k << ',' << i[i.size()-1].timeStampEnd << '\n';
  }

  return 0;
};
