
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include "hungarian.h"
#include "mmt_partial_coloring.h"

class ColorSelfOrganizingMap {
public:
  ColorSelfOrganizingMap(int num_center, int max_dist) : num_center(num_center), max_dist(max_dist){
    //
  }

  void feed(PartialColoring* input){

    float cur_avg_dist;
    int nearest_center = getNearestNode(input, cur_avg_dist);

    // replace old center if new input fits better
    if (cur_avg_dist < avg_dists[nearest_center])
      replaceNode(input,nearest_center);
  }

  void test(PartialColoring* input) {
    float cur_avg_dist;
    int nearest_center = getNearestNode(input, cur_avg_dist);
    logger.winner[nearest_center]++;
  }

private:
    int num_center, max_dist;
    std::vector<PartialColoring> centers;
    std::vector<float> avg_dists;
    LogData logger;

    int getNearestCenter(PartialColoring* input, float& avg_dist){
      int min_approx_dist = max_dist;
      int min_exact_dist = max_dist;
      int nearest_center;
      for (int i = 0; i < n; i++) {
        int cur_approx_dist = distances.push_back(input->distanceTo(&nodes[i], false));
        avg_dist += cur_approx_dist;
        if (cur_approx_dist <= min_approx_dist) {
          int cur_exact_dist = distances.push_back(input->distanceTo(&nodes[i], true));
          if (cur_exact_dist < min_exact_dist) {
            // new nearest node found
            min_approx_dist = cur_approx_dist;
            min_exact_dist = cur_exact_dist;
            nearest_center = i;
          }
        }
      }
      // remove distance to nearest neighbor from avg_dist;
      avg_dist -= min_approx_dist;
      // calculate avg
      avg_dist /= (n-1);

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
