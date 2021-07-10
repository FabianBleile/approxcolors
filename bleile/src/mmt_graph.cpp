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
  read_graph(argc, av, &n, &m, &elist);

  // interpret every input graph as symmetric and if edges (i j) and (j i) are given handle it
  std::vector<std::unordered_set<nodeid> > adjListSets(n,std::unordered_set<nodeid>());
  for (size_t i = 0; i < m; i++) {
    adjListSets[(elist)[2*i]].insert((elist)[2*i + 1]);
    adjListSets[(elist)[2*i + 1]].insert((elist)[2*i]);
  }


  adjList = std::vector<std::vector<nodeid> >(n);
  for (size_t i = 0; i < n; i++) {
    adjList[i] = std::vector<nodeid>(adjListSets[i].begin(), adjListSets[i].end());
    std::sort(adjList[i].begin(), adjList[i].end());
  }
}

MMTGraph::MMTGraph(MMTGraph * input) {
  n = input->n;
  m = input->m;
  instance = input->instance;
  adjList = input->adjList;
}

bool MMTGraph::isAdj(const nodeid u, const nodeid v) const {
  assert(u != v && isValid(u) && isValid(v));
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
         std::cout << u << " " << v << " ; ";
      }
    }
    std::cout << "("<< getDegree(u) <<")" << '\n';
    maxLines--;
  }
}

bool MMTGraph::isValid(const nodeid u) const {
    return u >= 0 && u < n;
}
