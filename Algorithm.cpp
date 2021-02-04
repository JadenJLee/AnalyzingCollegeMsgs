#include "Algorithm.h"
#include "Graph.h"
#include "Edge.h"

#include <stdlib.h>
#include <fstream>

using std::vector;
using std::string;
using std::unordered_map;
using std::memcpy;
using std::replace;

// Turn file into vector of strings (each entry is a line of file)
vector<string> Algorithm::file_to_vector(const string & filename) {
	std::ifstream in(filename);
  vector<string> vec;
  // Check if object is valid
  if(!in)
  {
      std::cerr << "Cannot open the File : "<< filename <<std::endl;
      return vec;
  }

  string str;
  // Read the next line from File until it reaches the end.
  while (std::getline(in, str))
  {
      // Line contains string of length > 0 then save it in vector
      if(str.size() > 0)
          vec.push_back(str);
  }
  //Close The File
  in.close();

  return vec;
} 

// Take a string and splits it into a vector (by whitespace)
vector<string> Algorithm::file_line_to_vector(string s) {
  vector<string> result;
  std::istringstream iss(s);
  for(string s; iss >> s; ) 
    result.push_back(s); 

  return result;
}

Algorithm::Algorithm(const string & filename) : g_(true, false) {
  filename_ = filename;
  // Step 1: split file into file lines
  vector<string> file_vector = file_to_vector(filename);

  // Step 2: iterate through lines and split each into vector w/ 3 entries
  for (string line : file_vector) {
    vector<string> line_vec = file_line_to_vector(line);

    // Step 3: check if vertices in [0] and [1] positions exist
    // if doesn't exist: create new vertex in graph
    if (!g_.vertexExists(line_vec[0])) {
      g_.insertVertex(line_vec[0]);
      g_.addIndex(line_vec[0]);
    }
    if (!g_.vertexExists(line_vec[0])) {
      g_.insertVertex(line_vec[1]);
      g_.addIndex(line_vec[1]);
    }

    // Step 4: add edge between [0] and [1] and generate a random number for edge weight
    g_.insertEdge(line_vec[0], line_vec[1]);
    g_.setEdgeWeight(line_vec[0], line_vec[1], (std::rand() % 20) + 1); // random edge weight from 1 to 20
    g_.setEdgeLabel(line_vec[0], line_vec[1], line_vec[0] + " " + line_vec[1]);
  }
}

void Algorithm::DFSRecur(Vertex v, vector<bool> &visited) {
  
  visited[g_.getIndex(v)] = true; //marks this Vertex v as visited 
  std::cout << v << " "; //prints out vertex to console 
  auto adj = g_.getAdjacent(v);  //vector of adjacent vertices to given Vertex v  
  std::vector<Vertex>::iterator it;   //iterator to iterate thru adjacent vector list 
  for (it = adj.begin(); it != adj.end(); ++it) {
    if (!visited[g_.getIndex(*it)]) { //*it = Vertex; // checks if Vertex has not been visited yet 
      DFSRecur(*it, visited); //Recurs thru graph 
    }
  }
  
}

void Algorithm::DFS() {
  // traverse and print path of g_
  auto vertices = g_.getVertices(); //builds vector of all vertices in graph 
  auto Vsize = vertices.size();  //captures size of vertices vector 
  
  vector<bool> visited(Vsize); //initialize visited vector to size of vertices vector 
  auto v = g_.getStartingVertex(); //a vertex to begin the DFS traversal with 
  DFSRecur(v, visited); //DFS recursive function to recur thru the graph 

}

// return a vector containing ranks of people (each entry of vector is a person) from most popular to least
vector<string> Algorithm::pageRank() {
  vector<string> vertexList(g_.getVertices().size());
  // create map
  unordered_map<Vertex, unordered_map<Vertex, Edge>> graphAdjacencyList = g_.getAdjacencyList();

  // create adjacency matrix
  std::vector<std::vector<double>> adjacencyMatrix(g_.getVertices().size(), vector<double> (g_.getVertices().size(),0));
  vector<double> rank(g_.getVertices().size());
  for(auto it = graphAdjacencyList.begin(); it != graphAdjacencyList.end(); ++it) {
    vertexList[g_.getIndex(it->first)] = it->first;
  }
  // normalize
  double normalizationVal = 1 / adjacencyMatrix.size();

  // Step 1 Build transition matrix using list
  for (auto it = graphAdjacencyList.begin(); it != graphAdjacencyList.end(); ++it) {
    for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      int col = g_.getIndex(it->first);
      int row = g_.getIndex(it2->first);
      double val = it2->second.getWeight();
      adjacencyMatrix[row][col] = val;
    }
  }
  // Step 2: Normalizing transition Matrix to become Markov Matrix (Each column sums to 1 and where each matrix entry represents a probability)
  for(size_t x = 0; x < adjacencyMatrix.size(); ++x) {
    double colSum = 0;
    for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
      colSum += adjacencyMatrix[x][y];
    }
    if(colSum > 0) {
      for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
         adjacencyMatrix[x][y] = adjacencyMatrix[x][y] / colSum;
      }
    } else {
      //If there is a whole column full of 0s (Also called dangling nodes), then the column is normalized by having each entry equal to 1 / number of vertices
      for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
         adjacencyMatrix[x][y] = normalizationVal;
      }
    }
  }

  // Step 3: Recompute transition matrix to become page rank matrix
  // Using formula pageRankMatrix = alphaValue * transitionMatrix  + matrix filled with ones * (1 - alphaValue) / size of matrix
  for(size_t x = 0; x < adjacencyMatrix.size(); ++x) {
    for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
      double val = adjacencyMatrix[x][y];
      val = (alpha * val) + ((1.0 - alpha) / adjacencyMatrix.size());
      adjacencyMatrix[x][y] = val;
    }
  }

  // Creating a Current Power Matrix
  std::vector<std::vector<double>> currentPowerMatrix(g_.getVertices().size(), vector<double> (g_.getVertices().size(),0));
  for(size_t x = 0; x < adjacencyMatrix.size(); ++x) {
    for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
      currentPowerMatrix[x][y] = (x == y) ? 1.0 : 0.0;
    }
  }
  // Step 4: Until the norm squared difference of current power matrix is lower than tolerance, take matrix multiplication 
  bool hasConverged = false;
  int numSteps = 0;
  while(!hasConverged) {
    //Matrix Multiplication
    std::vector<std::vector<double>> productMatrix(g_.getVertices().size(), vector<double> (g_.getVertices().size(),0));
    for(size_t x = 0; x < adjacencyMatrix.size(); ++x) {
      for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
        double sum = 0.0;
        for(size_t z = 0; z < adjacencyMatrix.size(); ++z) {
          sum += currentPowerMatrix[x][z] * adjacencyMatrix[z][y];
        }
        productMatrix[x][y] = sum;
      }
    }

    currentPowerMatrix = productMatrix;
    double diff = 0.0;
    double diff_normSquared = 0.0;
    for(size_t y = 0; y < adjacencyMatrix.size(); ++y) {
      for(size_t x = 1; x < adjacencyMatrix.size(); ++x) {
        diff = (currentPowerMatrix[x][y] - currentPowerMatrix[0][y]);
        diff_normSquared += diff * diff;
      }
    }
    if (diff_normSquared < tolerance) {
      hasConverged = true;
    } else {
      ++numSteps;
    }
  }

  // Creating rank vector
  for(size_t x = 0; x < adjacencyMatrix.size(); ++x) {
    rank[x] = currentPowerMatrix[0][x];
  }

  //Sorts vertices by rank using bubble sort algorithm
  size_t i = 0;
  size_t j = 1;
  bool canSort = true;
  while(canSort) {
    if(rank[i] < rank[j]) {
      double temp_i, temp_j;
      temp_i = rank[i];
      temp_j = rank[j];
      rank[i] = temp_j;
      rank[j] = temp_i;
      replace(vertexList.begin(), vertexList.end(), vertexList[i], vertexList[j]);
      i = 0;
      j = i;
    } else {
      ++i;
      ++j;
    }
    if(j >= adjacencyMatrix.size()) {
      canSort = false;
    }
  }

  // save results in file
  std::ofstream myfile;
  myfile.open("results/pagerank_results.txt");
  myfile << "The following are the ranked vertices:\n";
  for (size_t i = 0; i < vertexList.size(); i++) {
      myfile << vertexList[i] << "\n";
  }
  myfile.close();
  return vertexList;
}

// To represent Disjoint Sets 
// taken from geeksforgeeks
struct DisjointSets { 
    int *parent, *rank; 
    int n; 
    // unordered_map<Vertex, unordered_map<Vertex, Edge>> *parent;
    // unordered_map<Vertex, unordered_map<Vertex, Edge>> *rank;
  
    // Constructor. 
    DisjointSets(int n) { 
        // Allocate memory 
        this->n = n; 
        parent = new int[n+1]; 
        rank = new int[n+1]; 
  
        // Initially, all vertices are in 
        // different sets and have rank 0. 
        for (int i = 0; i <= n; i++) 
        { 
            rank[i] = 0; 
  
            //every element is parent of itself 
            parent[i] = i; 
        } 
    } 
  
    // Find the parent of a node 'u' 
    // Path Compression 
    int find(int elem) { 
        /* Make the parent of the nodes in the path 
           from u--> parent[u] point to parent[u] */
        if (elem != parent[elem]) 
            parent[elem] = find(parent[elem]); 
        return parent[elem]; 
    } 
  
    // Union by rank 
    void setUnion(int x, int y) { 
        x = find(x), y = find(y); 
  
        /* Make tree with smaller height 
           a subtree of the other tree  */
        if (rank[x] > rank[y]) 
            parent[y] = x; 
        else // If rnk[x] <= rnk[y] 
            parent[x] = y; 
  
        if (rank[x] == rank[y]) 
            rank[y]++; 
    } 
}; 

// 1. sort all edges in non-decreasing order of their weight
// 2. Pick the smallest edge, check if it forms a cycle
// 3. if no cycle, include it, if there is a cycle discard it.
// 4. repeat steps 2 and 3 until there are V-1 edges in the tree
vector<Edge> Algorithm::kruskals() {
  vector<Edge> output;

  vector<Edge> edges = g_.getEdges();
  vector<Vertex> vertices = g_.getVertices();

  // step 1. sort all edges in non-decreasing order of their weight
  std::sort(edges.begin(), edges.end());

  // create a disjoint sets
  DisjointSets disjointset(vertices.size());

  // step 2. pick the smallest edge
  for (auto it = edges.begin(); it != edges.end(); ++it) {
    Vertex nodeU = it->getSource();
    Vertex nodeV = it->getDest();

    int indexU = g_.getIndex(nodeU);
    int indexV = g_.getIndex(nodeV);

    int disjointsetU = disjointset.find(indexU);
    int disjointsetV = disjointset.find(indexV);

    // check if it forms a cycle
    if (disjointsetU != disjointsetV) {
      // push the nodeU and nodeV into output vector
      // idk how to do this
      
      disjointset.setUnion(disjointsetU, disjointsetV);
      output.push_back(*it); // push the edge back
    }
  }

  // save result in file
  std::ofstream myfile;
  myfile.open("results/kruskals_results.txt");
  myfile << "The following are the edges in the constructed MST:\n";
  int minimumCost = 0;
  for (size_t i = 0; i < output.size(); i++) {
      myfile << output[i].getSource() << " -- " << output[i].getDest()
            << " == " << output[i].getWeight() << "\n";
      minimumCost = minimumCost + output[i].getWeight();
  }
  myfile << "Minimum Cost of Spanning Tree: " << minimumCost << "\n";
  myfile.close();

  // return result of kruskal's algorithm
  return output;
}