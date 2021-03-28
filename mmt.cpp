extern "C" {
  #include "mmt_read.h"
}

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <numeric>

class MMTGraph {
public:
  const int n, m;

  MMTGraph(const int pncount, const int pecount, int **pelist) : n(pncount), m(pecount), adjList(pncount,std::unordered_set<int>()) {
    for (size_t i = 0; i < pecount; i++) {
      this->adjList[(*pelist)[2*i]].insert((*pelist)[2*i + 1]);
      this->adjList[(*pelist)[2*i + 1]].insert((*pelist)[2*i]);
    }
  }

  bool isAdj(const int u, const int v) const {
    assert(u != v && isValid(u) && isValid(v));
    return this->adjList[u].find(v) != this->adjList[u].end();
  }

  const std::unordered_set<int>* getNeighbors(const int u) const {
    assert(isValid(u));
    return &this->adjList[u];
  }

  int getDegree(const int u) const {
    assert(isValid(u));
    return adjList[u].size();
  }

  void toString(int maxLines = 5, bool real = true) const {
    std::cout << "MMTGraph.toString() of " << this << '\n';
    assert(this->n != 0 && this->m != 0);
    std::cout << "n = " << this->n << " : m = " << this->m << '\n';
    for (auto u = 0; u < this->n; u++) {
      if(maxLines == 0) {
        std::cout << "\t..." << '\n';
        return;
      }
      for (auto &v : this->adjList[u]) {
        if (real || u < v) {
           std::cout << u << " " << v << '\t';
        }
      }
      std::cout << "("<< getDegree(u) <<")" << '\n';
      maxLines--;
    }
  }
private:
  std::vector<std::unordered_set<int> > adjList;

  bool isValid(const int u) const {
    return u >= 0 && u < this->n;
  }
};


// maintains k stable sets and one more which stores all uncolored nodes
class MMTPartialColoring {
public:
  // empty constructor
  MMTPartialColoring(const int k, const MMTGraph * igraph) : k(k), colors(), graph(igraph) {
    assert(k!=0);
    assert(igraph != NULL);
    // generate random seq of nodes
    std::vector<int> rand_vect(graph->n);
    std::iota(rand_vect.begin(), rand_vect.end(), 0);
    std::random_shuffle(rand_vect.begin(), rand_vect.end());
    // insert nodes in (k+1)-st bucket
    for (const auto &u : rand_vect) this->colors.emplace(u, k);

  }
  // clear (k+1)-st bucket in a greedy way
  void greedy() {
    // SEQ
    // randomized insertion might not affect hashing in unordered_map ... idk(?)
    for (auto &kvp : this->colors) kvp.second = findMinAvailableColor(kvp.first);
  }
  // clear (k+1)-st bucket in a dsatur way
  void dsatur(){
    // compute initial degrees in G
    std::vector<int> v_dsatur;
    dsatur_selectMaxNodes([] (const MMTGraph* context, int u) -> int { return static_cast<const MMTGraph*>(context)->getDegree(u); }, v_dsatur, this->graph);
    std::vector<int> deg_dsatur(this->graph->n, 0);
    for (size_t i = 0; i < this->graph->n; i++) {
      // compute max degree
      for (size_t i = 0; i < this->graph->n; i++) {
        /* code */
      }
    }
  }

  void dsatur_selectMaxNodes(std::function<int(const MMTGraph*, const int&)> getDegree, std::vector<int>& v, const MMTGraph* context = null) const {
    v.clear();
    v = {0};
    for (size_t i = 1; i < this->graph->n; i++) {
      if(getDegree(context, i) > getDegree(context, v[0])) {
        v.clear();
        v.push_back(i);
      } else if (getDegree(context, i) == getDegree(context, v[0])) {
        v.push_back(i);
      }
    }
  }

  int findMinAvailableColor(int u) {
    const std::unordered_set<int> * u_neighbors = this->graph->getNeighbors(u);
    std::vector<bool> colorIsAvailable(k+1,true);
    for (const auto &v : *u_neighbors) {
      assert(isValidColor(this->colors[v]));
      colorIsAvailable[this->colors[v]] = false;
    }
    for (size_t i = 0; i < this->k; i++) { if (colorIsAvailable[i]) return i; }
    return this->k;
  }

  std::vector<std::unordered_set<int> >* getStableSets() {
    std::vector<std::unordered_set<int> > * stableSets = new std::vector<std::unordered_set<int> >(k+1, std::unordered_set<int>());
    for (const auto &kvp : this->colors) {
      assert(isValidColor(kvp.second));
      stableSets->operator[](kvp.second).emplace(kvp.first);
    }
    return stableSets;
  }

  void toString(int maxLines = 10) {
    std::cout << "MMTPartialColoring.toString() of " << this << '\n';
    std::vector<std::unordered_set<int> > * stableSets = getStableSets();
    for (const auto &set : *stableSets) {
      if(maxLines == 0) {
        std::cout << "\t..." << '\n';
        return;
      }
      for (const auto& u : set) {
        std::cout << u << ' ';
      }
      std::cout << '\n';
      maxLines--;
    }
  }

private:
  const int k;
  std::unordered_map<int,int> colors;
  const MMTGraph * graph;

  bool isValidColor(int value) const {
    return value >= 0 && value < this->k + 1;
  }
};

int main(int argc, char **av) {
  int n = 0, m = 0;
  int *elist;
  read_graph(argc, av, &n, &m, &elist);

  MMTGraph g(n, m, &elist);

  g.toString();

  MMTPartialColoring pc(5, &g);

  pc.dsatur();

  return 0;
}
