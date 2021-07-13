#include "bleile/header/mmt_partial_coloring.h"

/*

    *******************************************************
    *******************************************************

                    PartialColoring

    *******************************************************
    *******************************************************

*/

PartialColoring::PartialColoring(const int k, MMTGraph * graph) : k(k), graph(graph), colors(graph->n) {
  assert(k!=0);
  assert(graph != NULL);

  // generate ascending seq of nodes
  std::vector<nodeid> v(graph->n);
  std::iota(v.begin(), v.end(), 0);
  // insert nodes in (k+1)-st bucket
  for (const auto &u : v) setColor(u, k);
}

int PartialColoring::distanceTo(PartialColoring* S, bool exact) {
  assert(k == S->k);

  // populate intersection matrix
  std::vector<std::vector<double> > matIntersec(k+1,std::vector<double>(k+1,graph->n));
  int num_uncolored = 0;
  for (size_t i = 0; i < colors.size(); i++) {
    matIntersec[colors[i]][S->colors[i]]--;
    if (S->colors[i] == k) num_uncolored++;
  }
  num_uncolored += k*graph->n - std::accumulate(matIntersec[k].begin(), matIntersec[k].end()-1, 0);
  // remove uncolored nodes from distance calculation
  matIntersec.pop_back();
  for (auto& vIntersec : matIntersec) {
    vIntersec.pop_back();
  }

  return exact ? exactDistance(matIntersec, num_uncolored) : approxDistance(matIntersec, num_uncolored);
}

// clear (k+1)-st bucket in a greedy way
bool PartialColoring::greedy(const vector<nodeid>& v) {
  // SEQ
  if (v.empty()) {
    std::vector<nodeid> temp(uncolored.begin(), uncolored.end());
    // obtain a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(temp.begin(), temp.end(), std::default_random_engine(seed));
    for (const auto &u : temp) setColor(u, findMinAvailableColor(u));
  } else if (v.size() == graph->n) {
    for (const auto &u : v) setColor(u, findMinAvailableColor(u));
  }

  return evaluate() == 0;
}

// clear (k+1)-st bucket in a dsatur way
bool PartialColoring::dsatur(){
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
    setColor(u, findMinAvailableColor(u));
    // update degrees in G and C
    dsatur_updateDeg(u, satdegree, freedegree);

  }
  return evaluate() == 0;
}

int PartialColoring::dsatur_selectMaxNode(const std::vector<nodeid>& shuffled_nodes, std::vector<int>& satdegree, std::vector<int>& freedegree) const {
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

void PartialColoring::dsatur_updateDeg(nodeid u, std::vector<int>& satdegree, std::vector<int>& freedegree){
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

int PartialColoring::findMinAvailableColor(nodeid u) {
  const std::vector<nodeid> * u_neighbors = graph->getNeighbors(u);
  std::vector<bool> colorIsAvailable(k+1,true);
  for (const auto &v : *u_neighbors) {
    colorIsAvailable[colors[v]] = false;
  }
  for (size_t i = 0; i < k; i++) { if (colorIsAvailable[i]) return i; }
  return k;
}

measure PartialColoring::evaluate() {
  measure cost = 0;
  for (const auto & u : uncolored) {
    int deg = graph->getDegree(u);
    cost += deg < k ? 0 : deg;
  }
  fitness = cost;
  return cost;
  // return uncolored.size();
}

void PartialColoring::setK(int k){
  if (k > this->k) {
    std::cout << "Diese Funktionalität muss bei Bedarf noch implementiert werden" << '\n';
  } else if (k < this->k){
    this->k = k;
    uncolored.clear();
    for (size_t i = 0; i < colors.size(); i++) {
      if (colors[i] >= k) {
        setColor(i, k);
      }
    }
  }
}

// requires move to not result in any conflicts
void PartialColoring::setColor(nodeid u, color c){
  if(c == k) uncolored.insert(u);
  else uncolored.erase(u);
  colors[u] = c;
}

void PartialColoring::moveToColor(nodeid u, color c) {
  const std::vector<nodeid> * u_neighbors = graph->getNeighbors(u);
  for (const auto &v : *u_neighbors) {
    if (colors[v] == c) setColor(v, k);
  }
  setColor(u,c);
}

int PartialColoring::getNumColors() const {
  color maxcolor = 0;
  for (size_t i = 0; i < colors.size(); i++) {
    if (colors[i] > maxcolor) {
      maxcolor = colors[i];
    }
  }
  return maxcolor + 1;
}

void PartialColoring::toString(int maxLines) const {
  std::cout << "MMTPartialColoring.toString() of " << this << '\n';

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

// distance implementation proposed by D.C. Porumbel, J.-K. Hao, and P. Kuntz
int PartialColoring::approxDistance(std::vector<std::vector<double> >& matIntersec, int num_uncolored){
  int max_cost = 0;
  for (auto& intersec : matIntersec) {
    max_cost += (int) *std::min_element(intersec.begin(), intersec.end());
  }
  return std::max(0, graph->n - num_uncolored - (k*graph->n - max_cost));
}

int PartialColoring::exactDistance(std::vector<std::vector<double> >& matIntersec, int num_uncolored){
  HungarianAlgorithm HungAlgo;
	vector<int> assignment;

	double r = HungAlgo.Solve(matIntersec, assignment);

  // for (size_t i = 0; i < assignment.size(); i++) {
  //   std::cout << assignment[i] << ',' << graph->n - matIntersec[i][assignment[i]] << '\t';
  // }
  // std::cout << '\n';

	return graph->n - num_uncolored - (k*graph->n - r);
}

bool PartialColoring::operator<(const PartialColoring& S) const
{
  return (this->fitness < S.fitness);
}




/*

    *******************************************************
    *******************************************************

                    MMTPartialColoring

    *******************************************************
    *******************************************************

*/




// empty constructor
MMTPartialColoring::MMTPartialColoring(const int k, MMTGraph * graph, int L, int T) : PartialColoring(k, graph), L(L), T(T) {
  // everything handled by the super constructor
}

bool MMTPartialColoring::crossover(MMTPartialColoring& S1, MMTPartialColoring& S2){
  assert(color_classes.size() == 0);
  assert(this->k == S1.k && this->k == S2.k);
  assert(this->graph == S1.graph && this->graph == S2.graph);

  S1.lockColoring();
  S2.lockColoring();

  // generate vectors with size of each color class for both parents
  std::vector<int> S1_cc(k), S2_cc(k);
  auto it = S1.color_classes.begin();
  std::generate(S1_cc.begin(), S1_cc.end(),[&] () mutable {return it != S1.color_classes.end() ? (*it++).size() : 0; } );
  it = S2.color_classes.begin();
  std::generate(S2_cc.begin(), S2_cc.end(),[&] () mutable {return it != S2.color_classes.end() ? (*it++).size() : 0; } );

  // init currentColor
  color cur_color = 0;

  // iterate until all color classes of child have been populated
  // or there are no more colored nodes in parents which havn't already been colored in childs coloring
  bool alternate = true;
  MMTPartialColoring* act_parent = &S1;
  std::vector<int>* act_parent_cc = &S1_cc;
  MMTPartialColoring* pas_parent = &S2;
  std::vector<int>* pas_parent_cc = &S2_cc;

  for (size_t cur_color = 0; cur_color < k; cur_color++) {
    // select parent color with the greatest remaining size
    color h = std::distance(act_parent_cc->begin(),std::max_element(act_parent_cc->begin(),act_parent_cc->end()));

    for (const auto &u : act_parent->color_classes[h]) {
      if (colors[u] == k) {
        setColor(u, cur_color);
        if (pas_parent->colors[u] < k) {
          pas_parent_cc->at(h)--;
        }
      }
    }
    act_parent_cc->at(h) = 0;

    // SelectParent() :
    if (alternate) {
      auto temp_parent = act_parent;
      act_parent = pas_parent;
      pas_parent = temp_parent;
      auto temp_parent_cc = act_parent_cc;
      act_parent_cc = pas_parent_cc;
      pas_parent_cc = temp_parent_cc;
    }
  }

  return evaluate() == 0;
}

// Mutation
// set a pair (node, color) tabu for T steps
bool MMTPartialColoring::tabuSearch(){

  /*
    Laufzeitüberlegung aktuell:
      - zufälligen Knoten wählen O(n)
      - costVect aufbauen für init O(k), populate O(delta), calculate h 2*O(k)
      - Fall Schritt tabu: O(delta)
    => O(n + delta + 3*k)
  */

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
    color h = -1;
    int K = graph->n * graph->n;
    std::vector<int> costs(k, 0);
    for (const auto &v : *graph->getNeighbors(u)) {
      if(colors[v] != k) {
        costs[colors[v]] += graph->getDegree(v);
      }
    }
    for (color c = 0; c < k; c++) {
      if (costs[c] == 0) {
        h = c;
        goto color_found;
      } else if (tabuList[u][c] >= it) {
        costs[c] += K;
      }
    }
    h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));

  color_found:
    // perform move (u,h)
    moveToColor(u, h);
    // add move to TabuList
    tabuList[u][h] = it + T;
  }

  // for (size_t i = 0; i < isChosen.size(); i++) {
  //   if(isChosen[i] == false) {
  //     std::cout << i << ' ';
  //   }
  // }
  // std::cout << '\n';

  return evaluate() == 0;
}

// Mutation MMT
// set a single node tabu for T steps independant of the color.
// minimize over |delta(V_{k+1})|
bool MMTPartialColoring::tabuSearchSimplified(){
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
    color h = -1;
    int K = graph->n * graph->n;
    std::vector<int> costs(k, 0);
    for (const auto &v : *graph->getNeighbors(u)) {
      if(colors[v] != k) {
        costs[colors[v]]++;
      }
    }
    for (color c = 0; c < k; c++) {
      if (costs[c] == 0) {
        h = c;
        goto color_found;
      } else if (tabuList[u][c] >= it) {
        costs[c] += K;
      }
    }
    h = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));

  color_found:
    // perform move (u,h)
    moveToColor(u, h);
    // add move to TabuList
    tabuList[u][h] = it + T;
  }

  return greedy();
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

void MMTPartialColoring::lockColoring(){
  if (color_classes.size() != 0) return;
  color_classes = std::vector<std::unordered_set<nodeid> >(k, std::unordered_set<nodeid>());
  for (size_t i = 0; i < colors.size(); i++) {
    if (colors[i] < k) color_classes[colors[i]].insert(i);
  }
}

bool MMTPartialColoring::checkColoring(){
  lockColoring();
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
  std::cout << "CALCULATED COLORING IS LEGAL" << '\n';
  return true;
}
