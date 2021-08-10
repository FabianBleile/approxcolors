#include "bleile/header/mmt_graph.h"

MMTGraph::MMTGraph(int argc, char **av) {
  // save graph instance
  char * file = av[1];
  std::stringstream filestream;
  filestream << file;
  while(std::getline(filestream, instance, '/')){
   // loop for last segment delimited by '/' there instance name sits
  }

  int *elist;
  // read_graph(argc, av, &n, &m, &elist);
  if(!writeToElist(file, &n, &m, &elist))
    std::cout << "Fehler beim Einlesen der Instanz." << '\n';

  initFromElist(n, m, elist);

  delete elist;
}

MMTGraph::MMTGraph(int ncount, int ecount, int *elist) {
  initFromElist(ncount, ecount, elist);
}

MMTGraph::MMTGraph(MMTGraph * input) {
  n = input->n;
  m = input->m;
  dens = input->dens;
  instance = input->instance;
  adjList = input->adjList;
}

void MMTGraph::initFromElist(int ncount, int ecount, int *elist){
  adjList = std::vector<std::vector<nodeid> >(n,std::vector<nodeid>());

  for (size_t i = 0; i < m; i++) {
    adjList[(elist)[2*i]].push_back((elist)[2*i + 1]);
    adjList[(elist)[2*i + 1]].push_back((elist)[2*i]);
  }

  dens = (float) m/(n*(n+1)/2);
}

bool MMTGraph::isAdj(const nodeid u, const nodeid v) const {
  if (!isValid(u) || !isValid(v) || u == v) {
    return false;
  }
  return std::find(adjList[u].begin(),adjList[u].end(),v) != adjList[u].end();
}

const std::vector<nodeid>* MMTGraph::getNeighbors(const nodeid u) const {
  assert(isValid(u));
  return &adjList[u];
}

int MMTGraph::getDegree(const nodeid u) const {
  assert(isValid(u));
  return adjList[u].size();
}

void MMTGraph::toString(int maxLines, bool real) const
{
  std::cout << "MMTGraph.toString() of " << this << '\n';
  assert(n != 0 && m != 0);
  std::cout << "n = " << n << " : m = " << m << '\n';
  for (auto u = 0; u < n; u++) {
    if(maxLines == 0) {
      std::cout << "\t..." << '\n';
      return;
    }
    for (auto &v : adjList[u]) {
      if (real || u < v) {
         std::cout << u << " " << v << '\t';
      }
    }
    std::cout << "("<< getDegree(u) <<")" << '\n';
    maxLines--;
  }
}

bool MMTGraph::isValid(const nodeid u) const {
    return u >= 0 && u < n;
}

bool MMTGraph::writeToElist(char *f, int *pncount, int *pecount, int **pelist) {

    std::ifstream graph_file(f);

    std::string line;
    size_t it = 0;
    for (size_t it = 0; getline(graph_file, line);) {
      std::istringstream iss(line);
      char ch;
      if (iss >> ch)
      {
          size_t from, to;
          std::string format;

          switch(ch) {
              case 'c':
                break;
              case 'p':
                  if ((*pncount)||(*pecount)) return false;
                  if (iss >> format >> (*pncount) >> (*pecount)) {
                      if ("edge" != format) return false;
                  }
                  *pelist = new int[2*(*pecount)];
                  break;
              case 'e':
                  if (iss >> from >> to){
                    (*pelist)[it++] = from-1;
                    (*pelist)[it++] = to-1;
                  }
                  break;
              default:
                  return false;
          }
      }
    }

    graph_file.close();
    return true;
}
