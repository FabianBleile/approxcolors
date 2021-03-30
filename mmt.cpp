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
#include <tuple>
#include <stdexcept>

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

  void crossover(const MMTPartialColoring* S1, const MMTPartialColoring* S2){
    assert(color_classes.size() == 0);
    assert(this->k == S1->k && this->k == S2->k);
    assert(this->graph == S1->graph && this->graph == S2->graph);

    // generate vectors with size of each color class for both parents
    std::vector<int> s1_c(k), s2_c(k);
    auto it = S1->color_classes.begin();
    std::generate(s1_c.begin(), s1_c.end(), [&] () mutable { return it != S1->color_classes.end() ? (*it++).size() : 0; });
    it = S2->color_classes.begin();
    std::generate(s2_c.begin(), s2_c.end(), [&] () mutable { return it != S2->color_classes.end() ? (*it++).size() : 0; });

    // init number of colored nodes from Parents not colored in child coloring
    int s1_n = std::accumulate(s1_c.begin(), s1_c.end(), 0);
    int s2_n = std::accumulate(s2_c.begin(), s2_c.end(), 0);

    // init currentColor
    color cur_color = 0;

    // iterate until all color classes of child have been populated
    // or there are no more colored nodes in parents which havn't already been colored in childs coloring
    while (cur_color < k && s1_n + s2_n > 0) {

      // SelectParent() :
      // A[0,1,2] containing pointer to graph coloring Si, vect si_c, int si_n for chosen parent
      // A[3,4,5] accordingly the other
      auto A = selectParent(S1, S2, &s1_c, &s2_c, &s1_n, &s2_n, cur_color);

      // select parent color with the greatest remaining size
      color h = std::distance(std::get<1>(A)->begin(), std::max_element(std::get<1>(A)->begin(), std::get<1>(A)->end()));

      // std::cout << "chosen color : " << h << " and h vect " << (std::get<0>(A)->color_classes[h]).size() << '\n';

      for (const auto &u : std::get<0>(A)->color_classes[h]) {

        // insert returns <it to element, bool if was inserted>
        // if inserted update si_c and si_n
        if (this->colors[u] == this->k) {
          this->colors[u] = cur_color;
          (*std::get<1>(A))[h]--;     // lower color class size of parent
          std::get<2>(A)--;           // lower total sum of nodes color in parent coloring but not in child coloring
          try {
            color u_col_non_parent = std::get<3>(A)->colors.at(u);
            if (u_col_non_parent < k) {
              (*std::get<4>(A))[u_col_non_parent]--;     // lower color class size of non parent
              std::get<5>(A)--;                                   // lower total sum of nodes color in non parent coloring but not in child coloring
            }
          } catch (std::out_of_range e) {
            assert(true == false);    // just don't jump in here please :D
          }
        }
      }

      // color class h of parent should be empty by now
      assert((*std::get<1>(A))[h] == 0);

      // move to next color
      cur_color++;
    }
    toString();
    tabuSearch();
  }

  std::tuple<const MMTPartialColoring*, std::vector<int>*, int*, const MMTPartialColoring*, std::vector<int>*, int* > selectParent(const MMTPartialColoring* s1, const MMTPartialColoring* s2, std::vector<int>* s1_c, std::vector<int>* s2_c, int* s1_n, int* s2_n, int cur_color){
    if ( ( *s1_n > 0 && *s2_n > 0 && cur_color % 2 ) || *s1_n > 0) {
      return std::make_tuple(s1, s1_c, s1_n, s2, s2_c, s2_n);
    } else {
      return std::make_tuple(s2, s2_c, s2_n, s1, s1_c, s1_n);
    }
  }

  void tabuSearch(){
    assert(L>0 && T>0);
    assert(color_classes.size() == 0);
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

    lockColoring();
  }

  // clear (k+1)-st bucket in a greedy way
  void greedy() {
    assert(color_classes.size() == 0);
    // SEQ
    std::vector<nodeid> v(uncolored.begin(), uncolored.end());
    std::shuffle(v.begin(), v.end(), std::default_random_engine());
    for (const auto &u : v) setColor(u, findMinAvailableColor(u));

    tabuSearch();
  }
  // clear (k+1)-st bucket in a dsatur way
  void dsatur(){
    assert(color_classes.size() == 0);
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

    tabuSearch();
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

  void toString(int maxLines = 10) const {
    std::cout << "MMTPartialColoring.toString() of " << this << '\n';

    for (size_t i = 0; i < k+1; i++) {
      if(maxLines == 0) {
        std::cout << "\t..." << '\n';
        return;
      }
      for (const auto &kvp : colors) if (kvp.second == i) std::cout << kvp.first << ' ';
      std::cout << '\n';
      maxLines--;
    }
  }

private:
  const int k;
  const int L, T;
  std::unordered_map<nodeid,color> colors;
  std::unordered_set<nodeid> uncolored;
  std::vector<std::unordered_set<nodeid> > color_classes;
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

  void lockColoring(){
    color_classes = std::vector<std::unordered_set<nodeid> >(k, std::unordered_set<nodeid>());
    for (const auto& kvp : colors) {
      if (kvp.second < k) color_classes[kvp.second].insert(kvp.first);
    }
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

  int GOAL = 6;

  MMTPartialColoring tbs(GOAL, &g, 3000, 100);
  tbs.tabuSearch();
  tbs.toString(GOAL+1);

  // MMTPartialColoring seq(GOAL, &g, 3000, 100);
  // seq.greedy();
  // seq.toString(GOAL+1);

  MMTPartialColoring dsatur(GOAL, &g, 3000, 100);
  dsatur.dsatur();
  dsatur.toString(GOAL+1);

  MMTPartialColoring crossover(GOAL, &g, 3000, 100);
  crossover.crossover(&tbs, &crossover);
  crossover.toString(GOAL+1);

  return 0;
}
