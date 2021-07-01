
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <numeric>
#include "bleile/utils/hungarian.h"
#include "bleile/header/mmt_partial_coloring.h"

class PartialColoringCluster {
public:
  MMTGraph * graph;

  PartialColoringCluster(const char * filename, MMTGraph * g) {
    graph = g;

    // Read from the text file
    std::ifstream clusterfile(filename);

    std::string strline;
    while (getline (clusterfile, strline)) {
      std::stringstream ss;
      ss << strline;
      char strlinetype;
      ss >> strlinetype;
      switch (strlinetype) {
        case 'c':{
          break;
        }
        case 'p':{
          // num_center n
          std::string strtemp;
          ss >> strtemp;
          num_center = stoi(strtemp);
          ss >> strtemp;
          k = stoi(strtemp);
          ss >> strtemp;
          n = stoi(strtemp);
          break;
        }
        case 'i':{
          std::vector<nodeid> coloring;
          std::string strnode;
          while (ss >> strnode) {
            coloring.push_back(stoi(strnode));
          }
          PartialColoring tempCol(k, graph);
          tempCol.greedy(coloring);
          centers.push_back(tempCol);
          break;
        }
        default:{
          break;
        }
      }
    }
    // Close the file
    clusterfile.close();
    for (size_t i = 0; i < num_center; i++) {
      logger.start_center.push_back(0);
      logger.end_center.push_back(0);
    }
    for (size_t i = 0; i < n; i++) {
      logger.nearest_hist.push_back(0);
      logger.tot_dist_hist.push_back(0);
    }
  }

  PartialColoringCluster(int num_center, int n, int k, MMTGraph * g) : num_center(num_center), n(n), k(k) {
    graph = g;

    // init centers
    for (size_t i = 0; i < num_center; i++) {
      PartialColoring col(k, graph);
      col.greedy();
      centers.push_back(col);

      logger.start_center.push_back(0);
      logger.end_center.push_back(0);
    }

    for (size_t i = 0; i < n; i++) {
      logger.nearest_hist.push_back(0);
      logger.tot_dist_hist.push_back(0);
    }

    for (size_t i = 0; i < num_center; i++) {
      int next_center_dist = n;
      for (size_t j = 0; j < num_center; j++) {
        if (i == j) continue;
        int cur_center_dist = centers[i].distanceTo(&centers[j], true);
        next_center_dist = std::min(cur_center_dist, next_center_dist); ;
      }
      next_center_dists.push_back(next_center_dist);
    }
  }

  void initFeeding(int timeLimit, int ratioLimit){
    std::vector<bool> ratio(ratioLimit, 1);
    int sumRatio = ratioLimit;

    clock_t t = clock();
    for (size_t it = 0; ((float) clock() - t)/CLOCKS_PER_SEC < timeLimit && sumRatio != 0 ; it++) {
      it %= ratioLimit;
      PartialColoring col(k, graph);
      col.greedy();
      bool setToNewCenter = feed(col);
      sumRatio += setToNewCenter - ratio[it];
      ratio[it] = setToNewCenter;
      if (it == 0) {
        std::cout << "# updates in the last " << ratioLimit << " iters : " << sumRatio << " in " << ((float) clock() - t)/CLOCKS_PER_SEC << " secs" <<'\n';
      }
    }
    std::cout << "time " << ((float) clock() - t)/CLOCKS_PER_SEC << '\n';
    std::cout << "# updates in the last " << ratioLimit << " iters : " << sumRatio << '\n';
  }

  void initTesting(int timeLimit){
    clock_t t = clock();
    for (size_t it = 0; ((float) clock() - t)/CLOCKS_PER_SEC < timeLimit; it++) {
      PartialColoring col(k, graph);
      col.greedy();
      test(col, true);
    }
  }

  bool feed(PartialColoring& input){
    int nearest_center_dist;
    int next_center_dist;
    int nearest_center = getNearestCenter(input, nearest_center_dist, next_center_dist);

    if (next_center_dist > next_center_dists[nearest_center]) {
      replaceNode(input,nearest_center);
      return true;
    }
    return false;
  }

  void test(PartialColoring& input, bool start) {
    int nearest_center_dist;
    int next_center_dist;
    int nearest_center = getNearestCenter(input, nearest_center_dist, next_center_dist, true);
    if (nearest_center_dist > logger.max_nearest_dist) {
      logger.max_nearest_dist = nearest_center_dist;
    }
    if (start) {
      logger.start_center[nearest_center]++;
    } else {
      logger.end_center[nearest_center]++;
    }
    logger.nearest_hist[nearest_center_dist]++;
  }

  void writeAnalysisToFile(const char * filename){
    std::ofstream doc;
    doc.open (filename, std::fstream::app);
    doc << '\n';
    doc << "w ";
    for (auto& hits : logger.end_center) {
      doc << hits << ' ';
    }
    doc << '\n';
    doc << "w ";
    for (auto& hits : logger.nearest_hist) {
      if (hits != 0) {
        doc << hits << ' ';
      } else {
        doc << '.';
      }
    }
    doc << '\n';
    doc << "w ";
    for (auto& hits : logger.tot_dist_hist) {
      if (hits != 0) {
        doc << hits << ' ';
      } else {
        doc << '.';
      }
    }
    doc << '\n';
    doc.close();
  }

  void writeToFile(const char * filename){
    // char filename[ ] = "test.clu";
    std::ofstream doc;
    doc.open (filename, std::fstream::app);
    doc << "c Clustering for " << filename << " create by Fabian Bleile" << '\n';
    doc << "c Provides a search space spanning cluster for k = " << k << '\n';
    doc << "c # of tested partial colorings = " <<
      std::accumulate(logger.start_center.begin(),
                      logger.start_center.end(),
                      decltype(logger.start_center)::value_type(0)
                    )
        << '\n';
    doc << "c #center   k    n  max_nearest_dist center_updated\n";
    doc << "p " << num_center << ' ' << k << ' ' << n << ' ' << logger.max_nearest_dist << ' ' << logger.center_updated << '\n';
    for (auto& center : centers) {
      doc << "i ";
      for (size_t i = 0; i < k+1; i++) {
        for (int j = 0; j < center.colors.size(); j++) {
          if (center.colors[j] == i) {
            doc << j << ' ';
          }
        }
      }
      doc << '\n';
    }
    doc.close();
    // for (auto& center : centers) {
    //   center.toString();
    // }
  }

  void toString(bool extended) {
    std::cout << "cluster distances : " << '\n';
    int sum = 0;
    for (int i = 0; i < num_center; i++) {
      for (int j = 0; j < num_center; j++) {
        int approx = centers[i].distanceTo(&centers[j], false);
        int exact = centers[i].distanceTo(&centers[j],true);
        if (extended) std::cout << exact << ' ';
        sum += exact;
      }
      if (extended) std::cout << '\n';
    }
    std::cout << "avg = " << sum / (num_center*(num_center-1)) << '\n';
    std::cout << "max nearest_center_dist = " << logger.max_nearest_dist << '\n';
  }
  int k;

private:
  struct LogData {
    std::vector<size_t> start_center;
    std::vector<size_t> end_center;
    std::vector<size_t> nearest_hist;
    std::vector<size_t> tot_dist_hist;
    int max_nearest_dist = 0;
    size_t center_updated = 0;
  };

  int num_center, n;
  std::vector<PartialColoring> centers;
  std::vector<float> next_center_dists;
  LogData logger;

  int getNearestCenter(PartialColoring& input, int& nearest_center_dist, int& next_center_dist, bool update = false){
    // int min_approx_dist = n;
    int nearest_approx_dist = n;
    nearest_center_dist = n;
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
