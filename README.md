# CS 225 Final Project
**Authors**: Alice Huang, Jaden Lee, Sumedh Vaidyanathan, Rohan Sreerama

## About This Project
We created a graph data structure to assist in running algorithms on our chosen dataset (a college messaging network). The project is explained more in our GOALS document.

## Build
- To compile, type `make`
- To run the program, type `./finalproj`
- To compile the test suite, type `make test`
- To run the test suite, type `./test`

## Results
- To access results of running depth first search, kruskal, or page rank algorithm, navigate to the results folder. 
- To replicate the process we used to generate these results:
  - Run the program.
  - The kruskals_results.txt file should generate on its own.
  - The dfs_results.txt file contains the result that is printed to the console (which we copy-pasted into the txt file).
  - The pagerank_results.txt file should generate on its own. Note that this result takes quite a long time to finish running. If you just want to see the results of the first two txt files faster, just comment out the line of code that generates our pagerank results.
  - To speed up pagerank algorithm, make sure to change the value of the tolerance variable located in Algorithm.h to a larger value. Note: changing tolerance may result in an inaccurate popularity rank output
