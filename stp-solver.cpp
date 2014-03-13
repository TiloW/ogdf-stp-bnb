#include <iostream>
#include <chrono>
#include <limits>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/fileformats/GraphIO.h>

// temporary includes
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <sstream>
#include <string>

using namespace ogdf;
using namespace std::chrono;

/**
 * Used for comparing edges by their respective weight.
 */
class EdgeComparer
{
private:
	EdgeArray<double> m_origWeights;
	EdgeArray<edge> m_mapping;
public:
	EdgeComparer(const EdgeArray<edge> &mapping, const EdgeArray<double> &origWeights) : 
	  m_origWeights(origWeights), 
	  m_mapping(mapping) 
	{}
	int compare(const edge &e, const edge &f) const {
		// TODO: check for infinity?
		return m_origWeights[m_mapping[e]] - m_origWeights[m_mapping[f]];
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

void writeSVG(const Graph &graph, const List<node> &terminals, const std::string &name)
{
	GraphAttributes attr(graph,
	  GraphAttributes::edgeGraphics | 
	  GraphAttributes::nodeGraphics | 
	  GraphAttributes::nodeLabel | 
	  GraphAttributes::nodeStyle);

	node v;
	forall_nodes(v, graph) {
		attr.width(v) = attr.height(v) = 20;
		std::stringstream ss;
		ss << v;
		attr.label(v) = ss.str();
		if(terminals.search(v).valid()) {
			attr.fillColor(v) = Color::Red;
		}

	}

	FMMMLayout layout;
	layout.useHighLevelOptions(true);
	layout.unitEdgeLength(40);
	layout.newInitialPlacement(true);
	layout.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	layout.call(attr);
	GraphIO::drawSVG(attr, name);
}

edge moveTarget(Graph &graph, EdgeArray<edge> mapping, edge e, node newTarget) {
	edge origEdge = mapping[e];
	OGDF_ASSERT(origEdge != NULL);
	graph.delEdge(e);
	edge f = graph.newEdge(e->source(), newTarget);
	OGDF_ASSERT(f != NULL);
	mapping[f] = origEdge;

	return f;
}

edge moveSource(Graph &graph, EdgeArray<edge> mapping, edge e, node newSource) {
	edge origEdge = mapping[e];
	OGDF_ASSERT(origEdge != NULL);
	graph.delEdge(e);
	edge f = graph.newEdge(newSource, e->target());
	OGDF_ASSERT(f != NULL);
	mapping[f] = origEdge;

	return f;
}

double bnbInternal(
  Graph &graph, 
  EdgeArray<edge> mapping,
  const EdgeArray<double> &origWeights,
  List<node> &terminals,
  Graph &tree, 
  double upperBound,
  const EdgeComparer &comp,
  double prevCost = 0,
  int depth = 0)
{
	std::stringstream ss;
	ss << "svg/" << depth << ".svg";
	//writeSVG(graph, terminals, ss.str());

	for(int i = 0; i < depth; i++) cout << " ";
	cout << terminals.size() << "  |  " << prevCost << endl;

	OGDF_ASSERT(isLoopFree(graph));
	OGDF_ASSERT(isParallelFreeUndirected(graph));

	// validate mapping - only for debug
	edge f;
	forall_edges(f, graph) {
		OGDF_ASSERT(mapping[f] != NULL);
	}
	
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

					double tmp = origWeights[mapping[e2]] - origWeights[mapping[e1]];
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

		OGDF_ASSERT(graph.consistencyCheck());
		OGDF_ASSERT(isParallelFreeUndirected(graph));		
		OGDF_ASSERT(branchingEdge == NULL || mapping[branchingEdge] != NULL);
	
		// branching edge has been found or there is no feasible solution	
		if(branchingEdge != NULL && origWeights[mapping[branchingEdge]] < std::numeric_limits<double>::max()) {
			// remove branching edge
			graph.hideEdge(branchingEdge);

			for(int i = 0; i < depth; i++) cout << " ";
			cout << "branching edge " << branchingEdge << endl;

			// first branch: Inclusion of the edge
			// remove source node of edge and calculate new edge weights
			List<edge> edges, hiddenEdges;
			node nodeToRemove = branchingEdge->source();
			node targetNode = branchingEdge->target();
	
			OGDF_ASSERT(targetNode != NULL);
			OGDF_ASSERT(nodeToRemove != NULL);
			OGDF_ASSERT(graph.searchEdge(targetNode, nodeToRemove) == NULL);
			OGDF_ASSERT(graph.searchEdge(nodeToRemove, targetNode) == NULL);
			graph.adjEdges(nodeToRemove, edges);	
			forall_listiterators(edge, it, edges) {
				OGDF_ASSERT(branchingEdge != *it);	
				OGDF_ASSERT((*it)->opposite(nodeToRemove) != targetNode);
			}

			List<edge> movedEdges;
			edge e;
			forall_adj_edges(e, nodeToRemove) {			
				for(int i = 0; i <= depth; i++) cout << " ";
				cout << "moving edge " << e << endl;
	
				OGDF_ASSERT(e != branchingEdge);
				OGDF_ASSERT(e->target() == nodeToRemove || e->source() == nodeToRemove);
				OGDF_ASSERT(e->opposite(nodeToRemove) != targetNode);

				node w = e->opposite(nodeToRemove);
				edge f = graph.searchEdge(w, targetNode);
				if(f == NULL) {
					f = graph.searchEdge(targetNode, w);
				}
				if(f != NULL) {
					if(origWeights[mapping[f]] < origWeights[mapping[e]]) {
						hiddenEdges.pushFront(e);
						graph.hideEdge(e);
						OGDF_ASSERT(graph.consistencyCheck());
					}
					else {
						hiddenEdges.pushFront(f);
						graph.hideEdge(f);
						OGDF_ASSERT(graph.consistencyCheck());
					}
				}
				if(e->target() == nodeToRemove) {
					OGDF_ASSERT(e->target() == nodeToRemove);
					OGDF_ASSERT(e->source() != targetNode);
					//graph.moveTarget(e, targetNode);
					movedEdges.pushFront(moveTarget(graph, mapping, e, targetNode));
					OGDF_ASSERT(graph.consistencyCheck());
				}
				else {
					OGDF_ASSERT(e->source() == nodeToRemove)
					OGDF_ASSERT(e->target() != targetNode)
					//graph.moveSource(e, targetNode);
					movedEdges.pushFront(moveSource(graph, mapping, e, targetNode));
					OGDF_ASSERT(graph.consistencyCheck());
				}
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
			result = bnbInternal(
			  graph, 
			  mapping,
			  origWeights, 
			  terminals, 
			  tree, 
			  upperBound, 
			  comp, 
			  origWeights[mapping[branchingEdge]] + prevCost, 
			  depth+1);
			
			// restore previous graph	
			if(remNodeIsTerminal) {
				terminals.pushFront(nodeToRemove);
			}
			if(!targetNodeIsTerminal) {
				terminals.del(terminals.search(targetNode));
			}
		       
			forall_listiterators(edge, it, movedEdges) {
				edge e = *it;
				OGDF_ASSERT(e->source() != nodeToRemove && e->target() != nodeToRemove);
				
				if(e->target() == targetNode) {
					moveTarget(graph, mapping, e, nodeToRemove);
				}
				else {
					OGDF_ASSERT(e->source() == targetNode);
					moveSource(graph, mapping, e, nodeToRemove);
				}
			}

			forall_listiterators(edge, it, hiddenEdges) {
				graph.restoreEdge(*it);
			}

			// sencond branch: Exclusion of the edge
			double exEdgeResult = bnbInternal(
			  graph,
			  mapping, 
			  origWeights, 
			  terminals, 
			  tree, 
			  upperBound, 
			  comp, 
			  prevCost, 
			  depth+1);

			// decide which branch returned best result
			if(exEdgeResult < result) {
				result = exEdgeResult;
			}

			// finally: restore the branching edge
			graph.restoreEdge(branchingEdge);
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
	
	GraphCopy copiedGraph(graph);
	EdgeArray<edge> mapping(copiedGraph);
	List<node> copyTerminals;

	edge e;
	forall_edges(e, copiedGraph) {
		OGDF_ASSERT(copiedGraph.original(e));
		mapping[e] = copiedGraph.original(e);
	}


	forall_listiterators(node, it, terminals) {
		node v = copiedGraph.copy(*it);
		copyTerminals.pushFront(v);
	}
	
	EdgeComparer comp(mapping, weights);
	
	return bnbInternal(
	  copiedGraph, 
	  mapping,
	  weights, 
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
