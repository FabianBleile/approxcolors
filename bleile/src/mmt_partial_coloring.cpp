#include "bleile/header/mmt_partial_coloring.h"

/*

    *******************************************************
    *******************************************************

                    PartialCol

    *******************************************************
    *******************************************************

*/

PartialCol::PartialCol(const int k, Graph * graph) : k(k), graph(graph), colors(graph->n) {
  assert(k!=0);
  assert(graph != NULL);

  // generate ascending seq of nodes
  std::vector<nodeid> v(graph->n);
  std::iota(v.begin(), v.end(), 0);
  // insert nodes in (k+1)-st bucket
  for (const auto &u : v) setColor(u, k);
}

/*
    Calculate the distance between two PartialCols.
    Either exact or lose
*/
int PartialCol::distanceTo(PartialCol* S, bool exact) {
  assert(k == S->k);

  // populate intersection matrix
  std::vector<std::vector<double> > mat_intersec(k+1,std::vector<double>(k+1,graph->n));
  int num_uncolored = 0;
  for (size_t i = 0; i < colors.size(); i++) {
    mat_intersec[colors[i]][S->colors[i]]--;
    if (S->colors[i] == k) num_uncolored++;
  }
  num_uncolored += k*graph->n - std::accumulate(mat_intersec[k].begin(), mat_intersec[k].end()-1, 0);
  // remove uncolored nodes from distance calculation
  mat_intersec.pop_back();
  for (auto& vIntersec : mat_intersec) {
    vIntersec.pop_back();
  }

  return exact ? exactDistance(mat_intersec, num_uncolored) : approxDistance(mat_intersec, num_uncolored);
}

/*
    SEQ : Try to assign a color to all uncolored vertices
    Either by given order or random order
*/
bool PartialCol::greedy(const vector<nodeid>& v) {
  // SEQ
  if (v.empty()) {
    std::vector<nodeid> temp(uncolored.begin(), uncolored.end());
    // obtain a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(temp.begin(), temp.end(), std::default_random_engine(seed));
    for (const auto &u : temp) setColor(u, getMinAvailableColor(u));
  } else if (v.size() == graph->n) {
    for (const auto &u : v) setColor(u, getMinAvailableColor(u));
  }

  return evaluate() == 0;
}

/*
    DSatur : Try to assign a color to all uncolored vertices
*/
bool PartialCol::dsatur(){
  // shuffle uncolored nodes
  std::vector<nodeid> v(uncolored.begin(), uncolored.end());
  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));

  // compute initial degrees in G
  std::vector<int > satdegree(graph->n, 0);
  std::vector<int > freedegree(graph->n);
  for (nodeid u = 0; u < graph->n; u++) {
    freedegree[u] = graph->getDegree(u);
  }
  for (nodeid i = 0; i < graph->n; i++) {
    // compute random node with max degree in G from all nodes with max saturation
    nodeid u = dsatur_selectMaxNode(v, satdegree, freedegree);
    // color max_node with the lowest available color
    setColor(u, getMinAvailableColor(u));
    // update degrees in G and C
    dsatur_updateDeg(u, satdegree, freedegree);

  }
  return evaluate() == 0;
}

/*
    Return Node with max saturation degree for DSatur
*/
int PartialCol::dsatur_selectMaxNode(const std::vector<nodeid>& shuffled_nodes, std::vector<int>& satdegree, std::vector<int>& freedegree) const {
  int maxsatdegree_node = 0;
  for (auto i : shuffled_nodes) {
    if( (satdegree[i] > satdegree[maxsatdegree_node]) ||
          ((satdegree[maxsatdegree_node] == satdegree[i]) &&
          (freedegree[i] > freedegree[maxsatdegree_node])) ) {

      maxsatdegree_node = i;

    }
  }
  // obtain a time-based seed:
  return maxsatdegree_node;
}

/*
    Update saturation and free degrees for DSatur
*/
void PartialCol::dsatur_updateDeg(nodeid u, std::vector<int>& satdegree, std::vector<int>& freedegree){
  // set satdegree of chosen node to -1 in order to not select again
  satdegree[u] = -1;
  // get neighbors of u and update sat and freedegree of neighbors
  color u_color = colors[u];
  const std::vector<nodeid> * u_neighbors = graph->getNeighbors(u);
  for (const auto &v : *u_neighbors) {
    if (satdegree[v] != -1) {
      freedegree[v]--;
      bool new_color = true;
      const std::vector<nodeid> * v_neighbors = graph->getNeighbors(v);
      for (const auto &w : *v_neighbors) {
        if (colors[w] == u_color && w != u) {
          new_color = false;
          break;
        }
      }
      if (new_color) {
        satdegree[v]++;
      }
    }
  }
}

/*
    Given Node u; return minimum color with no conflicting nodes
*/
int PartialCol::getMinAvailableColor(nodeid u) {
  const std::vector<nodeid> * u_neighbors = graph->getNeighbors(u);
  std::vector<bool> colorIsAvailable(k+1,true);
  for (const auto &v : *u_neighbors) {
    colorIsAvailable[colors[v]] = false;
  }
  for (size_t i = 0; i < k; i++) { if (colorIsAvailable[i]) return i; }
  return k;
}

/*
    Evaluate Partial Coloring by the sum over the degrees of uncolored vertices
*/
int PartialCol::evaluate() {
  int cost = 0;
  for (const auto & u : uncolored) {
    int deg = graph->getDegree(u);
    cost += deg < k ? 0 : deg;
  }
  fitness = cost;
  return cost;
  // return uncolored.size();
}

/*
    Migrate PratialCol from k color classes to k_new classes
*/
/*
void EvolPartialCol::setK(int k_new){
  if (k_new > this->k) {
    std::cout << "requesting to increase k - this functionality is not implemented yet" << '\n';
  } else if (k_new < this->k){
    this->k = k_new;
    uncolored.clear();
    for (size_t i = 0; i < colors.size(); i++) {
      if (colors[i] >= k_new) {
        setColor(i, k_new);
      }
    }
  }
}
*/


void EvolPartialCol::setK(int k_new){
  if (k_new > this->k) {
    std::cout << "requesting to increase k - this functionality is not implemented yet" << '\n';
  } else if (k_new < this->k){
    buildColorClasses();
    std::vector<int> cc_cost(k, 0);
    for (size_t u = 0; u < colors.size(); u++) {
      if (colors[u] < k) {
        cc_cost[colors[u]] += graph->getDegree(u);
      }
    }
    this->k = k_new;
    for (size_t u = 0; u < colors.size(); u++) {
      setColor(u, k_new);
    }
    for (size_t i = 0; i < k_new; i++) {
      color best_color = std::distance(cc_cost.begin(), std::max_element(cc_cost.begin(), cc_cost.end()));
      for (auto u : this->color_classes[best_color]) {
        setColor(u, i);
      }
      cc_cost[best_color] = 0;
    }
  }
}

/*
    Set color of given node u to c
    Attention: there is no checking for conflicts at this point
    Require: move u to color c is not raising conflicts
*/
void PartialCol::setColor(nodeid u, color c){
  if(c == k) uncolored.insert(u);
  else uncolored.erase(u);
  colors[u] = c;
}

/*
    Set color of given node u to c while moving conflicting nodes to uncolored
*/
void PartialCol::moveToColor(nodeid u, color c) {
  assert(0 <= c && c <= k);
  const std::vector<nodeid> * u_neighbors = graph->getNeighbors(u);
  for (const auto &v : *u_neighbors) {
    if (colors[v] == c) setColor(v, k);
  }
  setColor(u,c);
}

/*
    count used colors
*/
int PartialCol::getNumColors() const {
  int numColors = 0;
  std::vector<bool> colorUsed(k, false);
  for (size_t i = 0; i < colors.size(); i++) {
    if (colors[i] < k && !colorUsed[colors[i]]) {
      colorUsed[colors[i]] = true;
      numColors++;
    }
  }
  return numColors;
}

/*
    Outout Coloring
    - ordered by color clases
    - sorted by numerical
*/
void PartialCol::toString(int maxLines) const {
  std::cout << "EvolPartialCol.toString() of " << this << '\n';

  if (this->uncolored.empty()) {
    std::cout << "successfully colored with " << getNumColors() << " colors." << '\n';
  }

  for (size_t i = 0; i < k+1; i++) {
    if(maxLines == 0) {
      std::cout << "\t..." << '\n';
      return;
    }
    std::cout << "Color " << i << ": ";
    for (size_t j = 0; j < colors.size(); j++) if (colors[j] == i) std::cout << j << ' ';
    std::cout << '\n';
    maxLines--;
  }
}

/*
    distance implementation proposed by D.C. Porumbel, J.-K. Hao, and P. Kuntz
*/
int PartialCol::approxDistance(std::vector<std::vector<double> >& mat_intersec, int num_uncolored){
  int max_cost = 0;
  for (auto& intersec : mat_intersec) {
    max_cost += (int) *std::min_element(intersec.begin(), intersec.end());
  }
  return std::max(0, graph->n - num_uncolored - (k*graph->n - max_cost));
}

/*
    HungarianAlgorithm to calculate set-theoretic partition distance
*/
int PartialCol::exactDistance(std::vector<std::vector<double> >& mat_intersec, int num_uncolored){
  HungarianAlgorithm HungAlgo;
	vector<int> assignment;

	double r = HungAlgo.Solve(mat_intersec, assignment);

	return graph->n - num_uncolored - (k*graph->n - r);
}

/*
    compare fitness of two PartialCols
*/
bool PartialCol::operator<(const PartialCol& S) const
{
  return (this->fitness < S.fitness);
}




/*

    *******************************************************
    *******************************************************

                    EvolPartialCol

    *******************************************************
    *******************************************************

*/




/*
    Empty constructor taking k, graph and tabuSearch parameters L,T
*/
EvolPartialCol::EvolPartialCol(const int k, Graph * graph) : PartialCol(k, graph) {
  // everything handled by the super constructor
}

/*
    GPX crossover between two PartialCols
*/
bool EvolPartialCol::crossover(EvolPartialCol& S1, EvolPartialCol& S2){
  assert(color_classes.size() == 0);
  assert(this->k == S1.k && this->k == S2.k);
  assert(this->graph == S1.graph && this->graph == S2.graph);

  S1.buildColorClasses();
  S2.buildColorClasses();

  // generate vectors with size of each color class for both parents
  std::vector<int> s1_c(k), s2_c(k);
  auto it = S1.color_classes.begin();
  std::generate(s1_c.begin(), s1_c.end(), [&] () mutable { return it != S1.color_classes.end() ? (*it++).size() : 0; });
  it = S2.color_classes.begin();
  std::generate(s2_c.begin(), s2_c.end(), [&] () mutable { return it != S2.color_classes.end() ? (*it++).size() : 0; });

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
    auto A = selectParent(&S1, &S2, &s1_c, &s2_c, &s1_n, &s2_n, cur_color);

    // select parent color with the greatest remaining size
    color h = std::distance(std::get<1>(A)->begin(), std::max_element(std::get<1>(A)->begin(), std::get<1>(A)->end()));

    // std::cout << "chosen color : " << h << " and h vect " << (std::get<0>(A)->color_classes[h]).size() << '\n';

    for (const auto &u : std::get<0>(A)->color_classes[h]) {

      // insert returns <it to element, bool if was inserted>
      // if inserted update si_c and si_n
      if (this->colors[u] == this->k) {
        setColor(u,cur_color);
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

/*
    GPX crossover: select parent two get color class from
*/
std::tuple<const EvolPartialCol*, std::vector<int>*, int*, const EvolPartialCol*, std::vector<int>*, int* > EvolPartialCol::selectParent(const EvolPartialCol* s1, const EvolPartialCol* s2, std::vector<int>* s1_c, std::vector<int>* s2_c, int* s1_n, int* s2_n, int cur_color){
  if ( ( *s1_n > 0 && *s2_n > 0 && !(cur_color % 2) ) || *s2_n == 0) {
    return std::make_tuple(s1, s1_c, s1_n, s2, s2_c, s2_n);
  } else {
    return std::make_tuple(s2, s2_c, s2_n, s1, s1_c, s1_n);
  }
}

/*
    Tabu Search performed for L iterations with a tenure of T
*/
bool EvolPartialCol::tabuSearch(int L, int T){

  assert(L>0 && T>0);
  if (color_classes.size() != 0) {
    color_classes.clear();
  }

  // tabuList to store moves performed in recent history (iteration when move is valid again is stored)
  std::vector<std::vector<int>> tabuList(graph->n, std::vector<int>(k, 0));

  for (size_t it = 0; it < L; it++) {
    // solution discovered ? if yes break and return
    if (uncolored.empty()) break;
    // choose u from random vect
    auto random_it = std::next(std::begin(uncolored), (int) rand() % uncolored.size());
    nodeid u = *random_it;

    // init and populate cost vect, explore neighborhood for every color 0 to k-1
    // high enough constant to not use tabued moves
    nodeid h = getOptimalColor(u, tabuList, it);
    // perform move (u,h)
    moveToColor(u, h);
    // add move to TabuList
    tabuList[u][h] = it + T + (int) rand() % 10;
  }

  return evaluate() == 0;
}

/*
    Return best color concerning tabu-list and fitness
*/
int EvolPartialCol::getOptimalColor(nodeid u, std::vector<std::vector<int>>& tabuList, int it){
  assert(colors[u] == k);
  color h = -1;
  int K = graph->n * graph->n;
  std::vector<int> costs(k, -graph->getDegree(u));
  for (const auto &v : *graph->getNeighbors(u)) {
    if(colors[v] != k) {
      costs[colors[v]] += graph->getDegree(v);
    }
  }
  for (color c = 0; c < k; c++) {
    if (tabuList[u][c] > it) {
      costs[c] += K;
    }
  }
  return std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));
}

/*
    Priority Greedy : Given priority oerder; color nodes in a greedy way
*/
bool EvolPartialCol::priorityGreedy(const std::vector<int>& priority_v) {
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
    if ((float) rand()/RAND_MAX < 0.1) {
      std::swap(nodes_v.at(i), nodes_v.at(i+1));
    }
  }

  for (const auto &u : nodes_v) setColor(u, getMinAvailableColor(u));

  return evaluate() == 0;
}

/*
    populate color classes. coloring is stored in color vector
    transition that to color clases
*/
void EvolPartialCol::buildColorClasses(){
  if (color_classes.size() != 0) return;
  color_classes = std::vector<std::unordered_set<nodeid> >(k, std::unordered_set<nodeid>());
  for (size_t i = 0; i < colors.size(); i++) {
    if (colors[i] < k) color_classes[colors[i]].insert(i);
  }
}

/*
    check for conflicts in coloring
*/
bool EvolPartialCol::checkColoring(){
  buildColorClasses();
  for (auto cclass : color_classes) {
    for (auto u : cclass) {
      for (auto w : cclass) {
        if (u != w && graph->isAdj(u, w)) {
          std::cout << "ERROR - THIS COLORING IS NOT LEGAL!!!" << '\n';
          return false;
        }
      }
    }
  }
  std::cout << "coloring has no conflicts" << '\n';
  return true;
}
