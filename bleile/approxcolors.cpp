#include "header/mmt.h"

void documentation(std::string instance, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instance << ',' << i << '/' << imax <<',';
  doc << mmt->streamLogs().rdbuf() << '\n';
  doc.close();
}

void spaceAnalysis(){

}

int main(int argc, char **av) {

  MMTGraph g(argc, av);

  MMTPartialColoring explorer(11, &g, 1, 45);
  explorer.dsatur();
  std::vector<int> distDistrib = explorer.spaceAnalysis(250, 1000, 600);
  for (size_t i = 0; i < distDistrib.size(); i++) {
    std::cout << distDistrib[i] << ',';
  }
  std::cout << '\n';


  return 0;
};
