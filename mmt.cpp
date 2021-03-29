extern "C" {
  #include "mmt_read.h"
}

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <map>
#include <queue>

using nodeid = uint32_t;
using color = uint32_t;
using measure = int;

class MMTGraph {
public:
  const int n, m;

  MMTGraph(const int pncount, const int pecount, int **pelist) : n(pncount), m(pecount), adjList(pncount,std::unordered_set<nodeid>()) {
    for (size_t i = 0; i < pecount; i++) {
      adjList[(*pelist)[2*i]].insert((*pelist)[2*i + 1]);
      adjList[(*pelist)[2*i + 1]].insert((*pelist)[2*i]);
    }
  }

  bool isAdj(const nodeid u, const nodeid v) const {
    assert(u != v && isValid(u) && isValid(v));
    return adjList[u].find(v) != adjList[u].end();
  }

  const std::unordered_set<nodeid>* getNeighbors(const nodeid u) const {
    assert(isValid(u));
    return &adjList[u];
  }

  int getDegree(const nodeid u) const {
    assert(isValid(u));
    return adjList[u].size();
  }

  void toString(int maxLines = 5, bool real = true) const {
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
private:
  std::vector<std::unordered_set<nodeid> > adjList;

  bool isValid(const nodeid u) const {
    return u >= 0 && u < n;
  }
};


class MMTPartialColoring {
public:
  // empty constructor
  MMTPartialColoring(const int k, const MMTGraph * igraph, int L, int T) : k(k), colors(), graph(igraph), L(L), T(T) {
    assert(k!=0);
    assert(igraph != NULL);
    // generate ascending seq of nodes
    std::vector<nodeid> v(graph->n);
    std::iota(v.begin(), v.end(), 0);
    // insert nodes in (k+1)-st bucket
    for (const auto &u : v) setColor(u, k);
  }

  void tabuSearch(){
    assert(L>0 && T>0);
    std::queue<std::unordered_set<std::pair<nodeid, color>, UInt32PairHash>::iterator > tabuQueue;
    std::unordered_set<std::pair<nodeid, color>, UInt32PairHash> tabuList;

    for (size_t i = 0; i < L; i++) {
      if (uncolored.empty()) break;
      // randomly select an uncolored vertex v in V_(k+1)
      std::vector<nodeid> rand_uncolored_node;
      rand_uncolored_node.insert(rand_uncolored_node.begin(), uncolored.begin(), uncolored.end());
      std::shuffle(rand_uncolored_node.begin(), rand_uncolored_node.end(), std::default_random_engine());

      nodeid u = rand_uncolored_node[0];

      //explore neighborhood
      std::vector<int> costs(k+1, -graph->getDegree(u));
      const std::unordered_set<nodeid>* u_neighbors = graph->getNeighbors(u);
      for (const auto &v : *u_neighbors) costs[colors[v]] += graph->getDegree(v);
      color h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));
      assert(isValidColor(h));
      if (costs[h] >= 0 && tabuList.find(std::make_pair(u, h)) != tabuList.end()) {
        auto it = costs.begin();
        std::generate(costs.begin(), costs.end(), [&] () mutable {
          if (tabuList.find(std::make_pair(u, h)) != tabuList.end()) { // is tabu
            it++;
            return std::numeric_limits<int>::max();
          } else {
            return *it++;
          }
        });
        color temp_h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));
        if (costs[temp_h] < std::numeric_limits<int>::max()) h = temp_h;
      }
      tabuQueue.push(tabuList.insert(std::make_pair(u, h)).first);
      if (i >= T) {
        tabuList.erase(tabuQueue.front());
        tabuQueue.pop();
      }
      moveToColor(u, h);
    }
  }

  // clear (k+1)-st bucket in a greedy way
  void greedy() {
    // SEQ
    std::vector<nodeid> v(uncolored.begin(), uncolored.end());
    std::shuffle(v.begin(), v.end(), std::default_random_engine());
    for (const auto &u : v) setColor(u, findMinAvailableColor(u));
  }
  // clear (k+1)-st bucket in a dsatur way
  void dsatur(){
    // compute initial degrees in G
    std::vector<std::pair<int, int> > deg_dsatur(graph->n); // degU, degDsatur
    int u = 0;
    std::generate(deg_dsatur.begin(), deg_dsatur.end(), [&] () mutable { return std::make_pair(graph->getDegree(u++), 0); });
    for (nodeid i = 0; i < graph->n; i++) {
      // compute random node with max degree in G from all nodes with max saturation
      nodeid u = dsatur_selectMaxNode(deg_dsatur);
      // color max_node with the lowest available color
      setColor(u, findMinAvailableColor(u));
      // update degrees in G and C
      dsatur_updateSatDeg(u, deg_dsatur);
    }
  }

  int dsatur_selectMaxNode(std::vector<std::pair<int, int> >& degs) const {
    std::vector<nodeid> v = {0};
    for (size_t i = 1; i < graph->n; i++) {
      if(degs[v[0]].second < degs[i].second || (degs[v[0]].second == degs[i].second && degs[v[0]].first < degs[i].first)) {
        v.clear();
        v.push_back(i);
      } else if (degs[v[0]].second == degs[i].second && degs[v[0]].first == degs[i].first) {
        v.push_back(i);
      }
    }
    std::shuffle(v.begin(), v.end(), std::default_random_engine());
    return v[0];
  }

  void dsatur_updateSatDeg(nodeid u, std::vector<std::pair<int, int> >& degs){
    const std::unordered_set<nodeid> * u_neighbors = graph->getNeighbors(u);
    for (const auto &v : *u_neighbors) {
      if (degs[v].second != -1) {
        degs[v].first--; // remove from G (uncolored Graph)
        degs[v].second++; // add to C (colored Graph)
      }
    }
    degs[u].second = -1;
  }

  int findMinAvailableColor(nodeid u) {
    const std::unordered_set<nodeid> * u_neighbors = graph->getNeighbors(u);
    std::vector<bool> colorIsAvailable(k+1,true);
    for (const auto &v : *u_neighbors) {
      assert(isValidColor(colors[v]));
      colorIsAvailable[colors[v]] = false;
    }
    for (size_t i = 0; i < k; i++) { if (colorIsAvailable[i]) return i; }
    return k;
  }

  std::vector<std::unordered_set<nodeid> >* getStableSets() {
    std::vector<std::unordered_set<nodeid> > * stableSets = new std::vector<std::unordered_set<nodeid> >(k+1, std::unordered_set<nodeid>());
    for (const auto &kvp : colors) {
      assert(isValidColor(kvp.second));
      stableSets->operator[](kvp.second).emplace(kvp.first);
    }
    return stableSets;
  }

  void toString(int maxLines = 10) {
    std::cout << "MMTPartialColoring.toString() of " << this << '\n';
    std::vector<std::unordered_set<nodeid> > * stableSets = getStableSets();
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
  const int L, T;
  std::unordered_map<nodeid,color> colors;
  std::unordered_set<nodeid> uncolored;
  const MMTGraph * graph;

  measure evaluate() const {
    measure cost = 0;
    for (const auto & u : uncolored) {
      cost += graph->getDegree(u);
    }
    return cost;
  }

  bool isValidColor(color value) const {
    return value >= 0 && value < k + 1;
  }

  void setColor(nodeid u, color c){
    std::cout << "color " << u << " in color " << c << '\n';
    if(c == k) uncolored.insert(u);
    else uncolored.erase(u);
    colors[u] = c;
  }

  void moveToColor(nodeid u, color c) {
    const std::unordered_set<nodeid> * u_neighbors = graph->getNeighbors(u);
    for (const auto &v : *u_neighbors) if (colors[v] == c) colors[v] = k;
    uncolored.erase(u);
    colors[u] = c;
  }

  struct UInt32PairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const {
      assert(sizeof(std::size_t)>=8);  //Ensure that std::size_t, the type of the hash, is large enough
      //Shift first integer over to make room for the second integer. The two are
      //then packed side by side.
      return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
    }
  };
};

int main(int argc, char **av) {
  int n = 0, m = 0;
  int *elist;
  read_graph(argc, av, &n, &m, &elist);

  MMTGraph g(n, m, &elist);

  g.toString();

  MMTPartialColoring pc(5, &g, m, n);
  pc.tabuSearch();
  pc.toString();

  return 0;
}
