#include "mmt_partial_coloring.h"

// empty constructor
MMTPartialColoring::MMTPartialColoring(const int k, MMTGraph * igraph, int L, int T) : k(k), colors(), graph(igraph), L(L), T(T) {
  assert(k!=0);
  assert(igraph != NULL);
  // generate ascending seq of nodes
  std::vector<nodeid> v(graph->n);
  std::iota(v.begin(), v.end(), 0);
  // insert nodes in (k+1)-st bucket
  for (const auto &u : v) setColor(u, k);
}

bool MMTPartialColoring::crossover(MMTPartialColoring* S1, MMTPartialColoring* S2){
  assert(color_classes.size() == 0);
  assert(this->k == S1->k && this->k == S2->k);
  assert(this->graph == S1->graph && this->graph == S2->graph);

  S1->lockColoring();
  S2->lockColoring();

  // generate vectors with size of each color class for both parents
  std::vector<int> s1_c(k), s2_c(k);
  auto it = S1->color_classes.begin();
  std::generate(s1_c.begin(), s1_c.end(), [&] () mutable { return it != S1->color_classes.end() ? (*it++).size() : 0; });
  it = S2->color_classes.begin();
  std::generate(s2_c.begin(), s2_c.end(), [&] () mutable { return it != S2->color_classes.end() ? (*it++).size() : 0; });

  // init number of colored nodes from Parents not colored in child coloring
  int* s1_n = new int(std::accumulate(s1_c.begin(), s1_c.end(), 0));
  int* s2_n = new int(std::accumulate(s2_c.begin(), s2_c.end(), 0));

  // init currentColor
  color cur_color = 0;

  // iterate until all color classes of child have been populated
  // or there are no more colored nodes in parents which havn't already been colored in childs coloring
  while (cur_color < k && *s1_n + *s2_n > 0) {

    // SelectParent() :
    // A[0,1,2] containing pointer to graph coloring Si, vect si_c, int si_n for chosen parent
    // A[3,4,5] accordingly the other
    auto A = selectParent(S1, S2, &s1_c, &s2_c, s1_n, s2_n, cur_color);

    // select parent color with the greatest remaining size
    color h = std::distance(std::get<1>(A)->begin(), std::max_element(std::get<1>(A)->begin(), std::get<1>(A)->end()));

    // std::cout << "chosen color : " << h << " and h vect " << (std::get<0>(A)->color_classes[h]).size() << '\n';

    for (const auto &u : std::get<0>(A)->color_classes[h]) {

      // insert returns <it to element, bool if was inserted>
      // if inserted update si_c and si_n
      if (this->colors[u] == this->k) {
        this->colors[u] = cur_color;
        (*std::get<1>(A))[h]--;     // lower color class size of parent
        #pragma GCC diagnostic ignored "-Wunused-value"
        (*std::get<2>(A))--;           // lower total sum of nodes color in parent coloring but not in child coloring
        try {
          color u_col_non_parent = std::get<3>(A)->colors.at(u);
          if (u_col_non_parent < k) {
            (*std::get<4>(A))[u_col_non_parent]--;     // lower color class size of non parent
            (*std::get<5>(A))--;                                   // lower total sum of nodes color in non parent coloring but not in child coloring
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
  return evaluate() == 0;
}

std::tuple<const MMTPartialColoring*, std::vector<int>*, int*, const MMTPartialColoring*, std::vector<int>*, int* > MMTPartialColoring::selectParent(const MMTPartialColoring* s1, const MMTPartialColoring* s2, std::vector<int>* s1_c, std::vector<int>* s2_c, int* s1_n, int* s2_n, int cur_color){

  if ( ( *s1_n > 0 && *s2_n > 0 && !(cur_color % 2) ) || *s2_n == 0) {
    return std::make_tuple(s1, s1_c, s1_n, s2, s2_c, s2_n);
  } else {
    return std::make_tuple(s2, s2_c, s2_n, s1, s1_c, s1_n);
  }
}

// Mutation
// set a pair (node, color) tabu for T steps
bool MMTPartialColoring::tabuSearch(){
  // Umleitung aktiv
  // return tabuSearchSimplified();

  assert(L>0 && T>0);
  assert(color_classes.size() == 0);
  std::queue<std::pair<nodeid, color> > tabuQueue;
  std::unordered_set<std::pair<nodeid, color>, UInt32PairHash> tabuList;

  std::pair<MMTPartialColoring, measure > local_best = std::make_pair(*this, evaluate());

  for (size_t i = 0; i < L; i++) {
    // solution discovered ? if yes break and return
    if (uncolored.empty()) break;

    // choose u from random vect
    assert(uncolored.size() != 0);
    auto random_it = std::next(std::begin(uncolored), (int) rand() % uncolored.size());
    nodeid u = *random_it;

    // color u with first available color if possible
    color h = findMinAvailableColor(u);

    // in case it is the new solution has better fitness for any function I considered by now
    // so we can just add u to that color and continue tabu search with a next random uncolored node
    if (h != k)
    {
      setColor(u, h);
    }
    // h = k iff there is no free color available for u
    // assert for every neighbor j is f(S_j) >= f(S*)
    else
    {
      //explore neighborhood for every color 0 to k-1
      std::vector<int> costs(k, -graph->getDegree(u));
      for (const auto &v : *graph->getNeighbors(u)) if(colors[v] != k) costs[colors[v]] += graph->getDegree(v);
      h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));

      // If (node, color) is tabu then search for best not tabu color
      if (tabuList.find(std::make_pair(u, h)) != tabuList.end()) {
        auto it = costs.begin() - 1;
        std::generate(costs.begin(), costs.end(), [&] () mutable {
          it++;
          if (tabuList.find(std::make_pair(u, h)) != tabuList.end()) { // is tabu
            return std::numeric_limits<int>::max();
          } else {
            return *it;
          }
        });
        color temp_h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));

        // is any color not tabu ? if yes update target color
        // else proceed with previously calculated tabued color
        if (costs[temp_h] < std::numeric_limits<int>::max()) {
          h = temp_h;
        }
      }

      assert(h < k); // at this point we want to add u to a color (color h is the set of all uncolored vertices)
      moveToColor(u, h);
    }

    // std::cout << "(node, color) : " << u << " , " << h << '\n';

    tabuList.insert(std::make_pair(u, h));
    tabuQueue.push(std::make_pair(u, h));

    if (i >= T) {
      tabuList.erase(tabuQueue.front());
      tabuQueue.pop();
    }
  }

  return greedy();
}

// Mutation MMT
// set a single node tabu for T steps independant of the color.
// minimize over |delta(V_{k+1})|
bool MMTPartialColoring::tabuSearchSimplified(){
  assert(L>0 && T>0);
  assert(color_classes.size() == 0);
  std::queue<nodeid > tabuQueue;
  std::unordered_set<nodeid> tabuList;

  for (size_t i = 0; i < L; i++) {
    // solution discovered ? if yes break and return
    if (uncolored.size() == 0) break;

    // collect non tabu nodes
    nodeid u = -1;
    std::vector<nodeid> nonTabuNodes;
    for (auto id : uncolored) {
      if (tabuList.find(id) != tabuList.end())
        nonTabuNodes.push_back(id);
    }

    // choose non tabu node at random or random uncolored node
    if (!nonTabuNodes.empty()) {
      u = *std::next(std::begin(nonTabuNodes), (int) rand() % nonTabuNodes.size());
    } else {
      u = *std::next(std::begin(uncolored), (int) rand() % uncolored.size());
    }

    //explore neighborhood for every color 0 to k-1
    std::vector<int> costs(k, 0);
    for (const auto &v : *graph->getNeighbors(u)) if(colors[v] != k) costs[colors[v]] += graph->getDegree(v);
    color h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));

    moveToColor(u, h);

    tabuList.insert(u);
    tabuQueue.push(u);

    if (i >= T) {
      tabuList.erase(tabuQueue.front());
      tabuQueue.pop();
    }
  }

  return greedy();
}

// clear (k+1)-st bucket in a greedy way
bool MMTPartialColoring::greedy() {
  assert(color_classes.size() == 0);
  // SEQ
  std::vector<nodeid> v(uncolored.begin(), uncolored.end());
  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
  for (const auto &u : v) setColor(u, findMinAvailableColor(u));

  return evaluate() == 0;
}

// clear (k+1)-st bucket in a greedy way
bool MMTPartialColoring::priorityGreedy(const std::vector<int>& priority_v) {
  assert(color_classes.size() == 0);
  assert(priority_v.size() == graph->n);
  // generate ascending seq of nodes
  std::vector<nodeid> nodes_v(graph->n);
  std::iota(nodes_v.begin(), nodes_v.end(), 0);

  std::sort(nodes_v.begin(), nodes_v.end(), [&](const nodeid & left, const nodeid & right) -> bool {
    return priority_v.at(left) < priority_v.at(right);
  });

  // copy priority vector and add little noise
  for (size_t i = 0; i < nodes_v.size() - 1; i++) {
    if ((float) rand()/RAND_MAX < 0.5) {
      std::swap(nodes_v.at(i), nodes_v.at(i+1));
    }
  }

  for (const auto &u : nodes_v) setColor(u, findMinAvailableColor(u));

  return evaluate() == 0;
}


// clear (k+1)-st bucket in a dsatur way
bool MMTPartialColoring::dsatur(){
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
  return evaluate() == 0;
}

int MMTPartialColoring::dsatur_selectMaxNode(std::vector<std::pair<int, int> >& degs) const {
  std::vector<nodeid> v = {0};
  for (size_t i = 1; i < graph->n; i++) {
    if(degs[v[0]].second < degs[i].second || (degs[v[0]].second == degs[i].second && degs[v[0]].first < degs[i].first)) {
      v.clear();
      v.push_back(i);
    } else if (degs[v[0]].second == degs[i].second && degs[v[0]].first == degs[i].first) {
      v.push_back(i);
    }
  }
  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
  return v[0];
}

void MMTPartialColoring::dsatur_updateSatDeg(nodeid u, std::vector<std::pair<int, int> >& degs){
  const std::unordered_set<nodeid> * u_neighbors = graph->getNeighbors(u);
  for (const auto &v : *u_neighbors) {
    if (degs[v].second != -1) {
      degs[v].first--; // remove from G (uncolored Graph)
      degs[v].second++; // add to C (colored Graph)
    }
  }
  degs[u].second = -1;
}

int MMTPartialColoring::findMinAvailableColor(nodeid u) {
  const std::unordered_set<nodeid> * u_neighbors = graph->getNeighbors(u);
  std::vector<bool> colorIsAvailable(k+1,true);
  for (const auto &v : *u_neighbors) {
    assert(isValidColor(colors[v]));
    colorIsAvailable[colors[v]] = false;
  }
  for (size_t i = 0; i < k; i++) { if (colorIsAvailable[i]) return i; }
  return k;
}

void MMTPartialColoring::toString(int maxLines) const {
  std::cout << "MMTPartialColoring.toString() of " << this << '\n';

  if (this->uncolored.empty()) {
    std::cout << "successfully colored with " << k << " colors." << '\n';
  }

  for (size_t i = 0; i < k+1; i++) {
    if(maxLines == 0) {
      std::cout << "\t..." << '\n';
      return;
    }
    std::cout << "( color " << i << " )" << ':';
    for (const auto &kvp : colors) if (kvp.second == i) std::cout << kvp.first << ' ';
    std::cout << '\n';
    maxLines--;
  }
}

measure MMTPartialColoring::evaluate() const {
  return uncolored.size();
}

int MMTPartialColoring::getNumColors() const {
  color maxcolor = 0;
  for (const auto &kvp : colors) {
    if (kvp.second > maxcolor) {
      maxcolor = kvp.second;
    }
  }
  return maxcolor + 1;
}

bool MMTPartialColoring::isValidColor(color value) const {
  return value >= 0 && value < k + 1;
}

void MMTPartialColoring::setColor(nodeid u, color c){
  if(c == k) uncolored.insert(u);
  else uncolored.erase(u);
  colors[u] = c;
}

void MMTPartialColoring::moveToColor(nodeid u, color c) {
  const std::unordered_set<nodeid> * u_neighbors = graph->getNeighbors(u);
  for (const auto &v : *u_neighbors) {
    if (colors[v] == c) setColor(v, k);
  }
  setColor(u,c);
}

void MMTPartialColoring::lockColoring(){
  if (color_classes.size() != 0) return;
  color_classes = std::vector<std::unordered_set<nodeid> >(k, std::unordered_set<nodeid>());
  for (const auto& kvp : colors) {
    if (kvp.second < k) color_classes[kvp.second].insert(kvp.first);
  }
}

std::size_t MMTPartialColoring::UInt32PairHash::operator()(const std::pair<uint32_t, uint32_t> &p) const {
    assert(sizeof(std::size_t)>=8);  //Ensure that std::size_t, the type of the hash, is large enough
    //Shift first integer over to make room for the second integer. The two are
    //then packed side by side.
    return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
}

int MMTPartialColoring::distanceTo(MMTPartialColoring* S, bool exact) {
  assert(k == S->k);
  this->lockColoring();
  std::vector<std::vector<double> > matIntersec(k+1,std::vector<double>(k+1,graph->n));
  for (auto& kvp : colors) {
    matIntersec[kvp.second][S->colors[kvp.first]]--;
  }

  return exact ? exactDistance(matIntersec) : approxDistance(matIntersec);
}

// distance implementation proposed by D.C. Porumbel, J.-K. Hao, and P. Kuntz
int MMTPartialColoring::approxDistance(std::vector<std::vector<double> >& matIntersec){
  int max_cost = 0;
  for (auto& intersec : matIntersec) {
    max_cost += (int) *std::min_element(intersec.begin(), intersec.end());
  }
  return std::max(0, max_cost - k*graph->n);
}

int MMTPartialColoring::exactDistance(std::vector<std::vector<double> >& matIntersec){
  HungarianAlgorithm HungAlgo;
	vector<int> assignment;

	double r = HungAlgo.Solve(matIntersec, assignment);

	return r - k*graph->n;
}
