#include "../header/mmt_graph.h"

Graph::Graph(int argc, char **av) {
  // save graph instance
  char * file = av[1];
  std::stringstream filestream;
  filestream << file;

  while(std::getline(filestream, instance, '/')){
   // loop for last segment delimited by '/' there instance name sits
  }

  int *elist;

  if(!writeToElist(file, &n, &m, &elist))
    std::cout << "Fehler beim Einlesen der Instanz." << '\n';

  initFromElist(elist);

  delete elist;
}

Graph::Graph(int ncount, int ecount, int *elist) {
  assert(ncount != 0);
  assert(ecount != 0);
  assert(elist != NULL);
  n = ncount;
  m = ecount;
  initFromElist(elist);
}

Graph::Graph(Graph * input) {
  n = input->n;
  m = input->m;
  dens = input->dens;
  instance = input->instance;
  adjList = input->adjList;
}

void Graph::initFromElist(int *elist){
  adjList = std::vector<std::vector<nodeid> >(n,std::vector<nodeid>());

  for (size_t i = 0; i < m; i++) {
    adjList[(elist)[2*i]].push_back((elist)[2*i + 1]);
    adjList[(elist)[2*i + 1]].push_back((elist)[2*i]);
  }

  dens = (float) m/(n*(n+1)/2);
}

bool Graph::isAdj(const nodeid u, const nodeid v) const {
  if (!isValidNode(u) || !isValidNode(v) || u == v) {
    return false;
  }
  return std::find(adjList[u].begin(),adjList[u].end(),v) != adjList[u].end();
}

const std::vector<nodeid>* Graph::getNeighbors(const nodeid u) const {
  assert(isValidNode(u));
  return &adjList[u];
}

int Graph::getDegree(const nodeid u) const {
  assert(isValidNode(u));
  return adjList[u].size();
}

void Graph::toString(int maxLines, bool real) const
{
  std::cout << "Graph.toString() of " << this << '\n';
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

bool Graph::isValidNode(const nodeid u) const {
    return u >= 0 && u < n;
}

bool Graph::writeToElist(char *f, int *pncount, int *pecount, int **pelist) {

    std::ifstream graph_file(f);

    std::vector<std::vector<bool> > adjMat;

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
                  iss >> format >> (*pncount) >> (*pecount);
                  *pelist = new int[2*(*pecount)];
                  adjMat = std::vector<std::vector<bool> >((*pncount),std::vector<bool>((*pncount),false));
                  break;
              case 'e':
                  if (iss >> from >> to){
                    if (from < 0 || from > (*pncount) || to < 0 || to > (*pncount) ) {
                      break;
                    }
                    if (adjMat[from-1][to-1]) {
                      break;
                    } else {
                      adjMat[from-1][to-1] = true;
                      adjMat[to-1][from-1] = true;
                      (*pelist)[it++] = from-1;
                      (*pelist)[it++] = to-1;
                    }
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
