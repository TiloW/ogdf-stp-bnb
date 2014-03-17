#include <iostream>
#include <chrono>
#include <limits>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/fileformats/GraphIO.h>

#define MAX_WEIGHT std::numeric_limits<double>::max()

// temporary includes
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <sstream>
#include <string>

using namespace ogdf;
using namespace std::chrono;

/**
 * Prints a simple text describing how to use this program.
 */
void printHelp()
{
	std::cout << "You need to provide a STP file containing the "
	  "Steiner-tree problem. The algorithm will then find the optimal "
          "Solution and print the total cost." << endl;
}

/**
 * Prints the given Steiner tree problem.
 *
 * Only used for debugging.
 */
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

/**
 * Used to validate the current mapping of edges to orignal edges
 *
 * Only used for debugging.
 */
bool validateMapping(const Graph &graph, const EdgeArray<edge> &mapping, const EdgeArray<double> &weights)
{
	edge e;
	forall_edges(e, graph) {
		OGDF_ASSERT(mapping[e] != NULL);
		OGDF_ASSERT(weights[mapping[e]]);
	}
	return  true;
}

/**
 * Calculates the optimal Steinter tree recursivly.
 * Should not be called directly but by calcSolution.
 * Each edge is either included or excluded, which gives rise to up to two new branches in each step.
 * 
 * \param graph
 *	the graph to be examined, though not const, no node will be modified
 *	edges may be replaced in the process while the resulting graph will be similiar to the input
 * \param mapping
 *	maps each edge of the graph to an edge with weight stored in origWeights
 * \param origWeights
 *	the weight of the original edges
 * \param terminals
 *	a list of terminals, the recursion will stop when the number of terminals reaches one.
 * \param isTerminal
 *	maps nodes to whether they are a terminal or not
 * \param tree
 *	the resulting Steinertree with edges in the original graph, not yet implemented
 * \param upperBound
 *	the currently best known upper bound for an optimal solution
 * \param prevCost
 *	the cost accumulated in previous recursion steps (previously included edges)
 * \param depth
 *	the number of recursive descents, only used for debugging
 * 
 * \return
 *	the total weight of the optimal solution (including prevCost)
 *	note: This might be higher then the actual solution if no solution
 *	satisfying the upper bound can be found.
 */
double bnbInternal(
  Graph &graph, 
  EdgeArray<edge> &mapping,
  const EdgeArray<double> &origWeights,
  List<node> &terminals,
  NodeArray<bool> &isTerminal,
  Graph &tree, 
  double upperBound,
  double prevCost = 0,
  int depth = 0)
{
	double result = MAX_WEIGHT;
	
	if(prevCost < upperBound) {
		// TODO: Construct tree	
		if(terminals.size() < 2) {
			// all terminals are connected
			result = prevCost;
		} 
		else { 	
			edge branchingEdge = NULL;
			double maxPenalty = -1;

			// calculate penalties for nodes
			double sumOfMinWeights = 0; // b
			double sumOfMinTermWeights = 0; // c
			double absoluteMinTermWeight = MAX_WEIGHT;
			for(ListConstIterator<node> it = terminals.begin(); sumOfMinWeights < MAX_WEIGHT && it.valid(); ++it) {
				double minTermWeight = MAX_WEIGHT,
				       minWeight = MAX_WEIGHT,
				       secondMinWeight = MAX_WEIGHT;
				edge minEdge = NULL;

				// investigate all edges of each terminal
				// calculate lower boundary and find branching edge
				List<edge> adjEdges;
				graph.adjEdges(*it, adjEdges);
				for(ListConstIterator<edge> itEdge = adjEdges.begin(); itEdge.valid(); ++itEdge) {
					edge e = *itEdge;
					if(origWeights[mapping[e]] < minWeight) {
						secondMinWeight = minWeight;
						minWeight = origWeights[mapping[e]];
						minEdge = e;
					}
					else {
						if(origWeights[mapping[e]] < secondMinWeight) {
							secondMinWeight = origWeights[mapping[e]];
						}
					}
					
					if(isTerminal[e->opposite(*it)] && origWeights[mapping[e]] < minTermWeight) {
						minTermWeight = origWeights[mapping[e]];
						if(minTermWeight < absoluteMinTermWeight) {
							absoluteMinTermWeight = minTermWeight;
						}
					}
				}

				if(sumOfMinTermWeights < MAX_WEIGHT && minTermWeight < MAX_WEIGHT) {
					sumOfMinTermWeights += minTermWeight;
				}
				else {
					sumOfMinTermWeights = MAX_WEIGHT;
				}
				OGDF_ASSERT(absoluteMinTermWeight <= sumOfMinTermWeights);

				// is terminal isolated or has only one edge?
				// if so we can break here
				if(minWeight == MAX_WEIGHT || 
				   secondMinWeight == MAX_WEIGHT) {
					branchingEdge = minEdge;
					if(minWeight == MAX_WEIGHT) {
						sumOfMinWeights = MAX_WEIGHT;
					}
				} else {
					sumOfMinWeights += minWeight;
					// update branching edge if need be
					double penalty = secondMinWeight - minWeight;
					if(penalty > maxPenalty) {
						maxPenalty = penalty;
						branchingEdge = minEdge;
					}
				}
			}

			// compare bounds for this graph
			if(branchingEdge != NULL) {
				double maxCost = upperBound - prevCost;
				if(sumOfMinTermWeights < MAX_WEIGHT) {
					sumOfMinTermWeights -= absoluteMinTermWeight;
				}
				if(maxCost <= sumOfMinWeights && maxCost <= sumOfMinTermWeights) {
					branchingEdge = NULL;
				}
			}

			OGDF_ASSERT(branchingEdge == NULL || mapping[branchingEdge] != NULL);
			
			// branching edge has been found or there is no feasible solution
			if(branchingEdge != NULL && origWeights[mapping[branchingEdge]] < MAX_WEIGHT) {	
				// remove branching edge
				node nodeToRemove = branchingEdge->source();
				node targetNode = branchingEdge->target();
				edge origBranchingEdge = mapping[branchingEdge];
				graph.delEdge(branchingEdge);
				branchingEdge->~EdgeElement();

				// first branch: Inclusion of the edge
				// remove source node of edge and calculate new edge weights	
				OGDF_ASSERT(targetNode != NULL);
				OGDF_ASSERT(nodeToRemove != NULL);
				OGDF_ASSERT(graph.searchEdge(targetNode, nodeToRemove) == NULL);
				OGDF_ASSERT(graph.searchEdge(nodeToRemove, targetNode) == NULL);
				
				List<node> delEdges, movedEdges;
				List<edge> origDelEdges;
				ListConstIterator<edge> itNext;
				
				List<edge> edges;
				graph.adjEdges(nodeToRemove, edges);
				for(ListConstIterator<edge> it = edges.begin(); it.valid(); it = itNext) {
					edge e = *it;
					itNext = it.succ();

					OGDF_ASSERT(e != branchingEdge);
					OGDF_ASSERT(e->target() == nodeToRemove || e->source() == nodeToRemove);
					OGDF_ASSERT(e->opposite(nodeToRemove) != targetNode);

					node w = e->opposite(nodeToRemove);
					edge f = graph.searchEdge(w, targetNode);
					if(f == NULL) {
						f = graph.searchEdge(targetNode, w);
					}
					bool deletedEdgeE = false;
					if(f != NULL) {
						if(origWeights[mapping[f]] < origWeights[mapping[e]]) {
							delEdges.pushFront(e->target());
							delEdges.pushFront(e->source());
							origDelEdges.pushFront(mapping[e]);
							graph.delEdge(e);
							e->~EdgeElement();

							deletedEdgeE = true;
						}
						else {
							delEdges.pushFront(f->target());
							delEdges.pushFront(f->source());
							origDelEdges.pushFront(mapping[f]);
							graph.delEdge(f);
							e->~EdgeElement();
						}
					}
					if(!deletedEdgeE) {	
						if(e->target() == nodeToRemove) {
							OGDF_ASSERT(e->source() != targetNode);
							movedEdges.pushFront(e->source());
							graph.moveTarget(e, targetNode);
						}
						else {
							OGDF_ASSERT(e->source() == nodeToRemove)
							OGDF_ASSERT(e->target() != targetNode)
							movedEdges.pushFront(e->target());
							graph.moveSource(e, targetNode);
						}
					}
				}
				// nodeToRemove is isolated at this point
				// thus no need to actually remove it
				// (easier to keep track of CopyGraph mapping)
				
				// remove node from terminals too
				bool remNodeIsTerminal = isTerminal[nodeToRemove],
                                     targetNodeIsTerminal = isTerminal[targetNode];
				if(remNodeIsTerminal) {
					terminals.del(terminals.search(nodeToRemove));
					isTerminal[nodeToRemove] = false;
				}
				if(!targetNodeIsTerminal) {
					terminals.pushFront(targetNode);
					isTerminal[targetNode] = true;
				}

				// calculate result on modified graph
				result = bnbInternal(
				  graph, 
				  mapping,
				  origWeights, 
				  terminals, 
				  isTerminal,
				  tree, 
				  upperBound, 
				  origWeights[mapping[branchingEdge]] + prevCost, 
				  depth+1);
				
				// update upper bound according to inclusion branch
				if(result < upperBound) {
					upperBound = result;
				}
				
				// restore previous graph	
				if(remNodeIsTerminal) {
					terminals.pushFront(nodeToRemove);
					isTerminal[nodeToRemove] = true;
				}
				if(!targetNodeIsTerminal) {
					terminals.del(terminals.search(targetNode));
					isTerminal[targetNode] = false;
				}
			      
				// restore moved edges 
				while(!movedEdges.empty()) {
					node v = movedEdges.popFrontRet();

					edge e = graph.searchEdge(v, targetNode);
					if(e == NULL) {
						e = graph.searchEdge(targetNode, v);
					}
					OGDF_ASSERT(e != NULL);
					OGDF_ASSERT(e->opposite(targetNode) != nodeToRemove);
		
					if(e->source() == v) {
						graph.moveTarget(e, nodeToRemove);
					}
					else {
						graph.moveSource(e, nodeToRemove);
					}
				}

				// restored deleted edges
				while(!delEdges.empty()) {
					node source = delEdges.popFrontRet();
					node target = delEdges.popFrontRet();
					edge e = graph.newEdge(source, target);
					OGDF_ASSERT(!origDelEdges.empty());
					mapping[e] = origDelEdges.popFrontRet();
				}
				OGDF_ASSERT(origDelEdges.empty());
				
				// sencond branch: Exclusion of the edge
				double exEdgeResult = bnbInternal(
				  graph,
				  mapping, 
				  origWeights, 
				  terminals, 
				  isTerminal,
				  tree, 
				  upperBound,
				  prevCost, 
				  depth+1);

				// decide which branch returned best result
				if(exEdgeResult < result) {
					result = exEdgeResult;
				}

				// finally: restore the branching edge
				edge f = graph.newEdge(nodeToRemove, targetNode);
				mapping[f] = origBranchingEdge;
			}
		}
		//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
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
  Graph &tree)
{
	tree.clear();
	
	Graph copiedGraph;
	NodeArray<node> origNodes(copiedGraph);
	NodeArray<node> copiedNodes(graph);
	EdgeArray<edge> origEdges(copiedGraph);
	List<node> copiedTerminals;
	NodeArray<bool> isTerminal(copiedGraph, false);

	node v;
	forall_nodes(v, graph) {
		node copiedNode = copiedGraph.newNode();
		origNodes[copiedNode] = v;
		copiedNodes[v] = copiedNode;
	}

	edge e;
	forall_edges(e, graph) {
		origEdges[copiedGraph.newEdge(
		  copiedNodes[e->source()],
		  copiedNodes[e->target()]
		)] = e;
	}

	forall_listiterators(node, it, terminals) {
		copiedTerminals.pushFront(copiedNodes[*it]);
		isTerminal[copiedNodes[*it]] = true;
	}

	return bnbInternal(
	  copiedGraph, 
	  origEdges,
	  weights, 
	  copiedTerminals, 
	  isTerminal,
	  tree,
	  MAX_WEIGHT); 
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
