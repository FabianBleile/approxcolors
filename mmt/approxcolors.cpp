/**
    This file is part of exactcolors.

    exactcolors is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    exactcolors is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
    the implementation the hybrid evolutionary algorithm MMT by Malaguti,
    Monaci and Toth arises from Fabian Bleile as part of the bachelor's thesis
      Randomisierte Algorithmen für Graphenfärbung, University of Bonn, 2021.

    In addition to the MMT mainly two features are added
    1. a new distance measure on feasible partial colorings is introduced and
      algorithmically used to guarantee a diverse pool of colorings. the distance
      measure is based on the set-theoretic partition distance over color classes
      with 'special' treatment of the set of uncolored vertices.
    2. recycling of feasible partial (k+1) colorings, after having discovered a
      complete and feasible (k+1) coloring, for better runtime performance concerning the
      subsequent solving of the k-GCP.

    approxcolors.cpp is the entry to call the mmt isolated from exactcolors
*/

#include "header/mmt.h"

void documentation(std::string instanceName, MMT* mmt, int i, int imax){
  char filename[ ] = "mmt_documentation.txt";
  std::ofstream doc;
  doc.open (filename, std::fstream::app);
  doc << instanceName << ','; // << i << '/' << imax <<',';
  doc << mmt->streamLogs().rdbuf();
  doc.close();
}

int main(int argc, char **av) {

  Graph g(argc, av);

  std::vector<std::vector<MMT::kLogData> > v;
  std::vector<int> v_totTime;

  for (size_t i = 0; i < 1; i++) {
    MMT mmt(g, /*L*/ 10000,/*T*/ 45, /*time limit*/ 120, /*pool size*/ 10, /* lower bound */ 2); // Graph * graph, int L, int T, int time_limit_sec, int pool_size = 99, double pGreedy = 0.5

    mmt.start();

    documentation(g.instanceName, &mmt, i, 4);

    v.push_back(mmt.logger.kLogData);
    v_totTime.push_back(mmt.logger.totTimeInSec);
  }

  for (int i = 0; i < v.size(); i++) {
    std::cout << v[i][v[i].size()-1].k << ',' << v[i][v[i].size()-1].timeStampEnd - v[i][v[i].size()-2].timeStampEnd << ' ' << v_totTime[i] << '\n';

    for (size_t j = 0; j < v[i].size(); j++) {
      std::cout << v[i][j].k << ',' << v[i][j].timeStampEnd << '\n';
    }
  }

  return 0;
};
