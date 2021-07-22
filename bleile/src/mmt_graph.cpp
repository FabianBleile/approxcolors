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

  initFromElist(n, m, elist);

  cliqueTabuSearch(10000, 500);

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
  clique = input->clique;
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

const std::vector<nodeid>* MMTGraph::getClique() const {
  return &clique;
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

void MMTGraph::cliqueTabuSearch(int L, int T){
  // tabuList to store moves performed in recent history (iteration when move is valid again is stored)
  std::vector<size_t> tabuList(n, 0);
  // store best clique in graph member variable clique
  // for temp calculating use the following tclique
  int cliquefitness = 0;
  std::vector<bool> curclique(n, false);
  int bestcliquefitness = 0;
  std::vector<bool> bestclique = curclique;
  // num of conflicts in current clique 'curclique'
  std::vector<int> numfriends(n, 0);

  for (size_t it = 0; it < L; it++) {
    // find best next node
    nodeid curnode = -1;
    int curfitness = std::numeric_limits<int>::max();
    for (nodeid u = 0; u < n; u++) {
      if (!curclique[u]) {
        int tempfitness = (tabuList[u] > it)*n - numfriends[u];
        if (tempfitness < curfitness) {
          curnode = u;
          curfitness = tempfitness;
        }
      }
    }

    // remove conflicting vertices from clique
    for (nodeid w = 0; w < n; w++) {
      if (isAdj(curnode, w)) {
        numfriends[w]++;
      } else if (curclique[w]) {
        // remove from clique
        curclique[w] = 0;
        cliquefitness--;
        for (auto z : *getNeighbors(w)) {
          numfriends[z]--;
        }
      }
    }

    // add best vertex to clique
    tabuList[curnode] = it+T;
    curclique[curnode] = 1;
    cliquefitness++;

    if (cliquefitness > bestcliquefitness) {
      bestclique = curclique;
      bestcliquefitness = cliquefitness;
    }
  }

  for (size_t i = 0; i < n; i++) {
    if (bestclique[i]) {
      //std::cout << i << ' ';
      clique.push_back(i);
    }
  }

  std::cout << "clique found with k = " << clique.size() << '\n';

  return;
}
