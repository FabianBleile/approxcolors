/*

  MMTgraph.cpp : adjacence matrix graph implementation

*/

#include <iostream>
#include <vector>

class MMTgraph {
  public:
    MMTgraph(int *pncount, int *pecount, int **pelist) {

      std::vector<int> v(5,1);

      for (size_t i = 0; i < *pecount; i++) {
        std::cout << pelist[i] << '\n';
      }
    }
  private:
    std::vector<std::vector<bool> > adjMat;
    int numVertices;
};
