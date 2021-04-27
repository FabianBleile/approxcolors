
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include "hungarian.h"
#include "mmt_partial_coloring.h"

class PartialColoringCluster {
public:
  PartialColoringSOM(int num_center, int max_dist, const int k, MMTGraph * graph) : num_center(num_center), max_dist(max_dist){
    // init centers
    for (size_t i = 0; i < num_center; i++) {
      PartialColoring partCol(k, graph);
      partCol.greedy();
      centers.push_back(partCol);
    }

    for (size_t i = 0; i < num_center; i++) {
      float sum_dist = 0;
      for (size_t j = 0; j < num_center; j++) {
        if (i == j) continue;
        sum_dist += centers[i].distanceTo(&centers[j], false);
      }
      next_center_dists.push_back(sum_dist/(num_center));
    }
  }

  void feed(PartialColoring* input){

    int next_center_dist;
    int nearest_center = getNearestNode(input, next_center_dist);

    // replace old center if new input fits better
    if (next_center_dist > next_center_dists[nearest_center])
      replaceNode(input,nearest_center);
  }

  void test(PartialColoring* input) {
    int next_center_dist;
    int nearest_center = getNearestNode(input);
    logger.winner[nearest_center]++;
  }

private:
    int num_center, max_dist;
    std::vector<PartialColoring> centers;
    std::vector<float> next_center_dists;
    LogData logger;

    int getNearestCenter(PartialColoring* input, int& next_center_dist = 0){
      int min_approx_dist = max_dist;
      int min_exact_dist = max_dist;
      int nearest_center;
      for (int i = 0; i < n; i++) {
        int cur_approx_dist = distances.push_back(input->distanceTo(&nodes[i], false));
        if (cur_approx_dist <= min_approx_dist) {
          int cur_exact_dist = distances.push_back(input->distanceTo(&nodes[i], true));
          if (cur_exact_dist < min_exact_dist) {
            // new nearest node found
            next_center_dist = min_approx_dist;
            min_approx_dist = cur_approx_dist;
            min_exact_dist = cur_exact_dist;
            nearest_center = i;
          }
        }
      }

      return nearest_center;
    }

    void replaceNode(PartialColoring& new_center, int old_center){
      nodes[old_center] = new_center;
      return;
    }

    struct LogData {
      int n, max_dist;
      std::vector<size_t> winner;
    }
}
