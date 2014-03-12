#include <iostream>
#include <chrono>
#include <limits>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;
using namespace std::chrono;

/**
 * Used for comparing edges by their respective weight.
 */
class EdgeComparer
{
private:
	EdgeArray<double> m_weights;
public:
	EdgeComparer(EdgeArray<double> weights) : m_weights(weights) {}	
	int compare(const edge &e, const edge &f) const {
		// TODO: check for infinity?
		return m_weights[e] - m_weights[f];
	}
	OGDF_AUGMENT_COMPARER(edge)
};

/**
 * Outputs a text describing how to use this program.
 */
void printHelp()
{
	std::cout << "You need to provide a STP file containing the "
	  "Steiner-tree problem. The algorithm will then find the optimal "
          "Solution and print the total cost." << endl;
}

double bnbInternal(
  Graph &graph, 
  EdgeArray<double> &weights, 
  List<node> &terminals, 
  Graph &tree, 
  double upperBound,
  EdgeComparer &comp,
  double prevCost = 0)
{
	cout << terminals.size() << "  |  " << prevCost << endl;

	// TODO: Construct tree	
	double result = std::numeric_limits<double>::max();
	// TODO: compare bounds

	if(terminals.size() < 2) {
		// all terminals are connected
		result = prevCost;
	} 
	else {	
		edge branchingEdge = NULL;
		double maxPenalty = 0;

		// calculate penalties for nodes
		forall_listiterators(node, it, terminals) {
			List<edge> edges;
			graph.adjEdges(*it, edges);
			if(edges.size() > 0) {
				edges.quicksort(comp);
				edge e1 = *(edges.get(0));
				
				if(edges.size() > 1) {
					edge e2 = *(edges.get(1));

					double tmp = weights[e2] - weights[e1];
					if(tmp >= maxPenalty) {
						maxPenalty = tmp;
						branchingEdge = e1;
					}
				} else {
					if(!edges.empty()) {
						branchingEdge = e1;
					}
					break;
				}
			} else {
				break;
			}
		}
		
		// branching edge has been found or there is no feasible solution	
		if(branchingEdge != NULL && weights[branchingEdge] < std::numeric_limits<double>::max()) {
			// "remove" branching edge by assigning infinite weight
			double branchingEdgeWeight = weights[branchingEdge];
			weights[branchingEdge] = std::numeric_limits<double>::max();

			// first branch: Inclusion of the edge
			// remove source node of edge and calculate new edge weights
			List<edge> edges, hiddenEdges;
			node nodeToRemove = branchingEdge->source();
			node targetNode = branchingEdge->target();
		
			edge e;
			forall_adj_edges(e, nodeToRemove) {
				if(e->target() != nodeToRemove) {
					graph.reverseEdge(e);
				}

				edge f = graph.searchEdge(e->source(), targetNode);
				if(f == NULL) {
					f = graph.searchEdge(targetNode, e->source());
				}
				if(f != NULL) {
					if(weights[f] < weights[e]) {
						hiddenEdges.pushFront(e);
						graph.hideEdge(e);
					}
					else {
						hiddenEdges.pushFront(f);
						graph.hideEdge(f);
					}
				}
				graph.moveTarget(e, targetNode);
			}
			// nodeToRemove is isolated at this point
			// thus no need to actually remove it
			// (easier to keep track of CopyGraph mapping)
			
			// remove node from terminals tooo
			ListIterator<node> it  = terminals.search(nodeToRemove),
			                   it2 = terminals.search(targetNode);
			bool remNodeIsTerminal = it.valid();
			if(remNodeIsTerminal) {
				terminals.del(it);
			}
			bool targetNodeIsTerminal = it2.valid();
			if(!targetNodeIsTerminal) {
				terminals.pushFront(targetNode);
			}

			// calculate result on modified graph
			result = bnbInternal(graph, weights, terminals, tree, upperBound, comp, branchingEdgeWeight + prevCost);
			
			// restore previous graph	
			if(remNodeIsTerminal) {
				terminals.pushFront(nodeToRemove);
			}
			if(!targetNodeIsTerminal) {
				terminals.del(terminals.search(targetNode));
			}
		       
			forall_listiterators(edge, it, edges) {
				graph.moveTarget(*it, nodeToRemove);
			}

			forall_listiterators(edge, it, hiddenEdges) {
				graph.restoreEdge(*it);
			}

			// sencond branch: Exclusion of the edge
			double exEdgeResult = bnbInternal(graph, weights, terminals, tree, upperBound, comp, prevCost);

			// decide which branch returned best result
			if(exEdgeResult < result) {
				result = exEdgeResult;
			}

			// finally: restore the branching edge
			weights[branchingEdge] = branchingEdgeWeight;
		}
	}
	return result;
}

/**
 * Yields an optimal steiner tree for the given STP instance.
 * Applies Shore, Foulds and Gibbons' branch and bound algorithm.
 *
 * \param graph
 *	the Graph containing both steiner and terminal nodes
 * \param weights
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
double calcSolution(
  const Graph &graph, 
  const EdgeArray<double> &weights, 
  const List<node> &terminals, 
//  const NodeArray<bool> &isTerminal, 
  Graph &tree)
{
	tree.clear();
	
	GraphCopy copyGraph(graph);
	EdgeArray<double> copyWeights(copyGraph);
	List<node> copyTerminals;

	edge e;
	forall_edges(e, copyGraph) {
		copyWeights[e] = weights[copyGraph.original(e)];
	}

	forall_listiterators(node, it, terminals) {
		node v = copyGraph.copy(*it);
		copyTerminals.pushFront(v);
	}
	
	EdgeComparer comp(copyWeights);
	
	return bnbInternal(
	  copyGraph, 
	  copyWeights, 
	  copyTerminals, 
	  tree, 
	  std::numeric_limits<double>::max(), 
	  comp); 
}

/**
 * Checks the provided command line arguments for validity.
 * Executes the branch and bound algorithm and writes the result to SVG and
 * command line.
 */
int main(int argc, char* argv[])
{
	if(argc != 2) {
		printHelp();
	}
	else {
		EdgeWeightedGraph<double> graph;
		List<node> terminals;
		NodeArray<bool> isTerminal(graph);
		if(GraphIO::readSTP(graph, terminals, isTerminal, argv[1])) {
			Graph tree;
			EdgeArray<double> weights(graph);

			// EdgeWeightedGraph will be removed soon from OGDF
			edge e;
			forall_edges(e, graph) {
				weights[e] = graph.weight(e);
			}	

			cout << "starting algorithm.." << endl;

			high_resolution_clock::time_point startTime = high_resolution_clock::now();
			double totalCost = calcSolution(graph, weights, terminals, tree);
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
