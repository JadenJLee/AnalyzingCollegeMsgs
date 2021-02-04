#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "Graph.h"
#include "Edge.h"

using std::string;

class Algorithm {

public:

  
  /**
   * Initializes g_ for the algorithms.
   *
   * @param filename path of data file
   */
  Algorithm(const string & filename);

  /**
   * Executes a depth-first search traversal and prints to a txt file.
   */
  void DFS();
  void DFSRecur(Vertex v, vector<bool> &visited);

  /**
   * Helper function for DFS
   */
  void DFSUtil(Vertex v, vector<bool> &visited);

  /**
   * Executes the PageRank algorithm on our graph. Returns a vector containing
   * the "rank" of the vertices (students in our dataset).
   */
  std::vector<std::string> pageRank();

  /**
   * Executes Kruskal's Algorithm on our graph. Returns a vector containing
   * a minimum spanning tree (?).
   */
  std::vector<Edge> kruskals();


  /* HELPER FUNCTIONS TO READ IN FILE (will be used in Algorithm constructor to initialize g_) */

  // take in filename and turn into vector, where each entry is a file line
  std::vector<std::string> file_to_vector(const std::string & filename);

  // take in string and split by white space
  std::vector<std::string> file_line_to_vector(std::string s);

private:
  Graph g_;
  std::string filename_;
  double alpha = 0.85;
  //double tolerance = 0.0001;
  double tolerance = 1;

};