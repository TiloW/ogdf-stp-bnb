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
class AdjEdgeComparer
{
private:
	EdgeArray<double> m_weights;
public:
	AdjEdgeComparer(EdgeArray<double> weights) : m_weights(weights) {}	
	int compare(const adjEntry &adj1, const adjEntry &adj2) const {
		return m_weights[adj1->theEdge()] - m_weights[adj2->theEdge()];
	}
	OGDF_AUGMENT_COMPARER(adjEntry)
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
  EdgeArray<double> weights, 
  List<node> &terminals, 
  Graph &tree, 
  double upperBound,
  AdjEdgeComparer &comp,
  double prevCost = 0)
{
	// TODO: Construct tree	
	double result = std::numeric_limits<double>::max();
	// TODO: compare bounds

	cout << terminals.size() << endl;

	if(terminals.size() < 2) {
		result = prevCost;
	} 
	else {
		edge branchingEdge = NULL;
		double maxPenalty = 0;

		cout << "  calculating penalties" << endl;
		// calculate penalties for nodes
		forall_listiterators(node, it, terminals) {
			List<adjEntry> adjEntries;
			cout << "  getting entries.." << endl;
			graph.adjEntries(*it, adjEntries);
			cout << "  got entries!" << endl;

			if(adjEntries.size() > 0) {
				adjEntries.quicksort(comp);
				edge e1 = (*(adjEntries.get(0)))->theEdge();
				
				if(adjEntries.size() > 1) {
					edge e2 = (*(adjEntries.get(1)))->theEdge();

					double tmp = weights[e2] - weights[e1];
					if(tmp >= maxPenalty) {
						maxPenalty = tmp;
						branchingEdge = e1;
					}
				} else {
					if(!adjEntries.empty()) {
						branchingEdge = e1;
					}
					break;
				}
			} else {
				break;
			}
		}

		cout << "  before branch" << endl;

		// branching edge has been found or there is no feasible solution	
		if(branchingEdge != NULL && weights[branchingEdge] < std::numeric_limits<double>::max()) {
			// first branch: Include the edge
		
			cout << "  branching on edge: " << branchingEdge << endl;
			graph.hideEdge(branchingEdge);
	
			// remove source node of edge and calc new edge weights
			List<edge> edges, hiddenEdges;
			node nodeToRemove = branchingEdge->source();
			node targetNode = branchingEdge->target();
			cout << "  retrieving adjacent edges.." << endl;
			graph.adjEdges(nodeToRemove, edges);
			
			cout << "  starting edge iteration" << endl;
			forall_listiterators(edge, it, edges) {
				cout << "  investigating edge " << *it << endl;
				edge e = *it;	
				if(e->target() != nodeToRemove) {
					graph.reverseEdge(e);
				}

				edge f = graph.searchEdge(e->source(), targetNode);
				if(f != NULL) {
					if(weights[f] < weights[e]) {
						hiddenEdges.pushFront(e);
						graph.hideEdge(e);
					}
					else {
						hiddenEdges.pushFront(e);
						graph.hideEdge(f);
					}
				}
				graph.moveTarget(e, targetNode);
			}
			graph.delNode(nodeToRemove);

			cout << "  node removed from graph" << endl;
			
			// remove node from terminals too
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

			cout << "  node removed from terminal set" << endl;
			
			cout << "  entering inclusion branch" << endl;
			// calculate result on modified graph
			result = bnbInternal(graph, weights, terminals, tree, upperBound, comp, weights[branchingEdge] + prevCost);
			
			// restore previous graph	
			node restoredNode = graph.newNode();
			if(remNodeIsTerminal) {
				terminals.pushFront(restoredNode);
			}
			if(!targetNodeIsTerminal) {
				terminals.del(terminals.search(targetNode));
			}
		       
			cout << "  moving edges to old positions" << endl; 	
			forall_listiterators(edge, it, edges) {
				graph.moveTarget(*it, restoredNode);
			}

			cout << "  revealing edges" << endl;
			forall_listiterators(edge, it, hiddenEdges) {
				cout << "  restoring edge: " << *it << endl;
				graph.restoreEdge(*it);
				cout << "  edge " << *it << " restored." << endl;
			}

			cout << "  entering exlusion branch" << endl;
			// sencond branch: Exclude the edge
			double exEdgeResult = bnbInternal(graph, weights, terminals, tree, upperBound, comp, prevCost);

			// finally: restore the branching edge
			cout << "  restoring branchingEdge";
			graph.restoreEdge(branchingEdge);
			cout << "  rstep finished!" << endl;

			if(exEdgeResult < result) {
				result = exEdgeResult;
			}
		}
	}
	return result;
}

/**
 * Yields a optimal steiner tree for the given STP instance.
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
  const EdgeArray<double> weights, 
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
	
	AdjEdgeComparer comp(copyWeights);
	
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
