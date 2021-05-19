
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <numeric>
#include "bleile/utils/hungarian.h"
#include "bleile/header/mmt_partial_coloring.h"

class PartialColoringCluster {
public:
  PartialColoringCluster(int num_center, int max_dist, const int k, MMTGraph * graph) : num_center(num_center), max_dist(max_dist), k(k) {
    // init centers
    for (size_t i = 0; i < num_center; i++) {
      PartialColoring partCol(k, graph);
      partCol.greedy();
      centers.push_back(partCol);

      logger.winner.push_back(0);
    }

    for (size_t i = 0; i < max_dist; i++) {
      logger.nearest_hist.push_back(0);
      logger.tot_dist_hist.push_back(0);
    }

    for (size_t i = 0; i < num_center; i++) {
      float sum_dist = 0;
      for (size_t j = 0; j < num_center; j++) {
        if (i == j) continue;
        sum_dist += centers[i].distanceTo(&centers[j], false);
      }
      next_center_dists.push_back(sum_dist/(num_center));
    }

    std::cout << "FINISHED CONSTRUCTOR wiht center size " << centers.size() << '\n';
  }

  void feed(PartialColoring& input){
    int nearest_center_dist;
    int next_center_dist;
    int nearest_center = getNearestCenter(input, nearest_center_dist, next_center_dist);

    // replace old center if new input fits better
    if (next_center_dist > next_center_dists[nearest_center])
      replaceNode(input,nearest_center);
  }

  void test(PartialColoring& input) {
    int nearest_center_dist;
    int next_center_dist;
    int nearest_center = getNearestCenter(input, nearest_center_dist, next_center_dist, true);
    if (nearest_center_dist > logger.max_nearest_dist) {
      logger.max_nearest_dist = nearest_center_dist;
    }
    logger.winner[nearest_center]++;
    logger.nearest_hist[nearest_center_dist]++;
  }

  /*
  void readFromFile(char* filename){

  }
  */

  void writeToFile(){
    char filename[ ] = "test.clu";
    std::ofstream doc;
    doc.open (filename, std::fstream::app);
    doc << "c Clustering for " << filename << " create by Fabian Bleile" << '\n';
    doc << "c Provides a search space spanning cluster for k = " << k << '\n';
    doc << "c # of tested partial colorings = " <<
      std::accumulate(logger.winner.begin(),
                      logger.winner.end(),
                      decltype(logger.winner)::value_type(0)
                    )
        << '\n';
    doc << "c #center   k    max_dist=#nodes of graph  max_nearest_dist center_updated\n";
    doc << "p " << num_center << ' ' << k << ' ' << max_dist << ' ' << logger.max_nearest_dist << ' ' << logger.center_updated << '\n';
    doc << "w ";
    for (auto& hits : logger.winner) {
      doc << hits << ' ';
    }
    doc << '\n';
    doc << "w ";
    for (auto& hits : logger.nearest_hist) {
      doc << hits << ' ';
    }
    doc << '\n';
    doc << "w ";
    for (auto& hits : logger.tot_dist_hist) {
      doc << hits << ' ';
    }
    doc << '\n';
    // for (auto& center : centers) {
    //   doc << "i ";
    //   for (size_t i = 0; i < k+1; i++) {
    //     for (const auto& kvp : center.colors) {
    //       if (kvp.second == i) doc << kvp.first << ' ';
    //     }
    //   }
    //   doc << '\n';
    // }
    doc.close();
    // for (auto& center : centers) {
    //   center.toString();
    // }
  }

  void toString() {
    std::cout << "cluster distances : " << '\n';
    int sum = 0;
    for (int i = 0; i < num_center; i++) {
      for (int j = 0; j < num_center; j++) {
        int approx = centers[i].distanceTo(&centers[j], false);
        int exact = centers[i].distanceTo(&centers[j],true);
        // std::cout << exact << ' ';
        sum += exact;
      }
      // std::cout << '\n';
    }
    std::cout << "avg = " << sum / (num_center*(num_center-1)) << '\n';
    std::cout << "max nearest_center_dist = " << logger.max_nearest_dist << '\n';
  }

private:
  struct LogData {
    std::vector<size_t> winner;
    std::vector<size_t> nearest_hist;
    std::vector<size_t> tot_dist_hist;
    int max_nearest_dist = 0;
    size_t center_updated = 0;
  };

  int num_center, max_dist;
  std::vector<PartialColoring> centers;
  std::vector<float> next_center_dists;
  const int k;
  LogData logger;

  int getNearestCenter(PartialColoring& input, int& nearest_center_dist, int& next_center_dist, bool update = false){
    // int min_approx_dist = max_dist;
    int nearest_approx_dist = max_dist;
    nearest_center_dist = max_dist;
    std::vector<int> nearest_center;
    for (int i = 0; i < num_center; i++) {
      int cur_center_dist = input.distanceTo(&centers[i], true);
      if (update) logger.tot_dist_hist[cur_center_dist]++;
      if (cur_center_dist < nearest_center_dist) {
        // new nearest node found
        next_center_dist = nearest_center_dist;
        // min_approx_dist = cur_approx_dist;
        nearest_center_dist = cur_center_dist;
        nearest_center.clear();
        nearest_center.push_back(i);
      } else if (cur_center_dist == nearest_center_dist) {
        next_center_dist = nearest_center_dist;
        nearest_center.push_back(i);
      }
    }

    return *std::next(std::begin(nearest_center), (int) rand() % nearest_center.size());
  }

  void replaceNode(PartialColoring& new_center, int old_center){
    logger.center_updated++;
    centers[old_center] = new_center;
    return;
  }
};
