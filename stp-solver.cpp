#include <iostream>
#include <chrono>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;
using namespace std::chrono;

/**
 * Outputs a text describing how to use this program.
 */
void printHelp() {
	std::cout << "You need to provide a STP file containing the "
	  "Steiner-tree problem. The algorithm will then find the optimal "
          "Solution and print the total cost." << endl;
}

/**
 * Yields a optimal steiner tree for the given STP instance.
 * Applies Shore, Foulds and Gibbons' branch and bound algorithm.
 *
 * \param graph
 *	the Graph containing both steiner and terminal nodes
 * \param weight
 *	the weight of each edge
 * \param terminals
 *	the set of terminal nodes, as opposed to steiner nodes
 * \param isTerminal
 *	mapping for each node whether it is a terminal or a steiner node
 * \param tree
 *	will hold the resulting tree
 *
 * \return
 *	the total cost of the steiner tree
 */
double calcSolution(Graph &graph, EdgeArray<double> weight, List<node> &terminals, NodeArray<bool> &isTerminal,  Graph &tree) {
	return 0;
}

/**
 * Checks the provided command line arguments for validity.
 * Executes the branch and bound algorithm and writes the result to SVG and
 * command line.
 */
int main(int argc, char* argv[]) {
	if(argc != 2) {
		printHelp();
	}
	else {
		EdgeWeightedGraph<double> graph;
		List<node> terminals;
		NodeArray<bool> isTerminal(graph);
		if(GraphIO::readSTP(graph, terminals, isTerminal, argv[1])) {
			cout << "starting algorithm" << endl;
			high_resolution_clock::time_point startTime = high_resolution_clock::now();
			double totalCost = 0; // TODOT
			duration<double> timeSpan = duration_cast<duration<double>>(high_resolution_clock::now() - startTime);
			cout << "Calculation finished after " << timeSpan.count() << " seconds." << endl;
			cout << "The optimal solution costs " << totalCost << "." << endl;
		}
		else {
			std::cerr << "unable to read file: " << argv[1] << endl;
		}
	}
	return 0;
}
