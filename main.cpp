#include <iostream>

#include "Graph.h"
#include "Algorithm.h"

int main() {
  Algorithm a("CollegeMsgData.txt");
  a.DFS();
  a.kruskals();
	a.pageRank();
}