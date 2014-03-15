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
 * Used for comparing edges by their respective weight.
 *
class EdgeComparer
{
private:
	const EdgeArray<double> *m_pOrigWeights;
	const EdgeArray<edge> *m_pMapping;
public:
	EdgeComparer(const EdgeArray<edge> *mapping, const EdgeArray<double> *origWeights)
	{
		m_pMapping = mapping;
		m_pOrigWeights = origWeights;
	}

	int compare(const edge &e, const edge &f) const {
		// TODO: check for infinity?
		OGDF_ASSERT((*m_pMapping)[e] != NULL);
		OGDF_ASSERT((*m_pMapping)[f] != NULL);
		return (*m_pOrigWeights)[(*m_pMapping)[e]] - (*m_pOrigWeights)[(*m_pMapping)[f]];
	}
	OGDF_AUGMENT_COMPARER(edge)
};
*/

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

bool validateMapping(const Graph &graph, const EdgeArray<edge> &mapping, const EdgeArray<double> &weights)
{
	edge e;
	forall_edges(e, graph) {
		OGDF_ASSERT(mapping[e] != NULL);
		OGDF_ASSERT(weights[mapping[e]]);
	}
	return  true;
}

double bnbInternal(
  Graph &graph, 
  EdgeArray<edge> &mapping,
  const EdgeArray<double> &origWeights,
  List<node> &terminals,
  Graph &tree, 
  double upperBound,
  double prevCost = 0,
  int depth = 0)
{
	double result = MAX_WEIGHT;
	
	// TODO: compare bounds
	if(prevCost < upperBound) {
		/*
		std::stringstream ss;
		ss << "svg/" << depth << ".svg";
		writeSVG(graph, terminals, ss.str());


		for(int i = 0; i < depth; i++) cout << " ";
		cout << terminals.size() << "  |  " << prevCost << endl;

		OGDF_ASSERT(isLoopFree(graph));
		OGDF_ASSERT(isParallelFreeUndirected(graph));
		OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
		*/
		
		// TODO: Construct tree	
		if(terminals.size() < 2) {
			// all terminals are connected
			result = prevCost;
		} 
		else { 	
			edge branchingEdge = NULL;
			double maxPenalty = -1;

			// calculate penalties for nodes
			bool continueSearch = true;
			for(ListConstIterator<node> it = terminals.begin(); continueSearch && it.valid(); ++it) {
				double minWeight = MAX_WEIGHT,
				       secondMinWeight = MAX_WEIGHT;
				edge minEdge = NULL;

				// investigate all edges of each terminal
				List<edge> adjEdges;
				graph.adjEdges(*it, adjEdges);
				for(ListConstIterator<edge> itEdge = adjEdges.begin(); continueSearch && itEdge.valid(); ++itEdge) {
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
				}

				// is terminal isolated or has only one edge?
				// if so we can break here
				if(minWeight == MAX_WEIGHT || 
				   secondMinWeight == MAX_WEIGHT) {
					branchingEdge = minEdge;
					continueSearch = false;
				} else {
					// update branching edge if need be
					double penalty = secondMinWeight - minWeight;
					if(penalty > maxPenalty) {
						maxPenalty = penalty;
						branchingEdge = minEdge;
					}
				}
			}
			OGDF_ASSERT(branchingEdge == NULL || mapping[branchingEdge] != NULL);
			/*
			OGDF_ASSERT(graph.consistencyCheck());
			OGDF_ASSERT(isParallelFreeUndirected(graph));		
			OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
			*/
			// branching edge has been found or there is no feasible solution
			if(branchingEdge != NULL && origWeights[mapping[branchingEdge]] < MAX_WEIGHT) {
				//for(int i = 0; i < depth; i++) cout << " ";
				//cout << "branching edge " << branchingEdge << endl;
				
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
				
				List<edge> edges;
				graph.adjEdges(nodeToRemove, edges);
				/*
				forall_listiterators(edge, it, edges) {
					OGDF_ASSERT(branchingEdge != *it);
					OGDF_ASSERT((*it)->opposite(nodeToRemove) != targetNode);
				}
				for(int i = 0; i <= depth; i++) cout << " ";
				cout << "found " << edges.size() << " edges that need to be moved" << endl;
				*/
				List<node> delEdges, movedEdges;
				List<edge> origDelEdges, origMovedEdges;
				ListConstIterator<edge> itNext;
				for(ListConstIterator<edge> it = edges.begin(); it.valid(); it = itNext) {
					edge e = *it;
					itNext = it.succ();

					//for(int i = 0; i <= depth; i++) cout << " ";
					//cout << "checking edge " << e << endl;
		
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
						//for(int i = 0; i <= depth; i++) cout << " ";
						//cout << "> found conflicting edge: " << f << ".. deleting edge with higher weight" << endl;
						if(origWeights[mapping[f]] < origWeights[mapping[e]]) {
							delEdges.pushFront(e->target());
							delEdges.pushFront(e->source());
							origDelEdges.pushFront(mapping[e]);
							graph.delEdge(e);
							e->~EdgeElement();
							//OGDF_ASSERT(graph.consistencyCheck());

							deletedEdgeE = true;
						}
						else {
							delEdges.pushFront(f->target());
							delEdges.pushFront(f->source());
							origDelEdges.pushFront(mapping[f]);
							graph.delEdge(f);
							e->~EdgeElement();
							//OGDF_ASSERT(graph.consistencyCheck());
						}
					}
					if(!deletedEdgeE) {
						//for(int i = 0; i <= depth; i++) cout << " ";
						//cout << "> edge " << e << " was not deleted and will thus be moved" << endl;
						
						origMovedEdges.pushFront(mapping[e]);
						if(e->target() == nodeToRemove) {
							OGDF_ASSERT(e->source() != targetNode);
							movedEdges.pushFront(e->source());
							graph.moveTarget(e, targetNode);
							//OGDF_ASSERT(graph.consistencyCheck());
							//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
						}
						else {
							OGDF_ASSERT(e->source() == nodeToRemove)
							OGDF_ASSERT(e->target() != targetNode)
							movedEdges.pushFront(e->target());
							graph.moveSource(e, targetNode);
							//OGDF_ASSERT(graph.consistencyCheck());
							//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
						}
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
				OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
				result = bnbInternal(
				  graph, 
				  mapping,
				  origWeights, 
				  terminals, 
				  tree, 
				  upperBound, 
				  origWeights[mapping[branchingEdge]] + prevCost, 
				  depth+1);
				OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
				
				// update upper bound according to inclusion branch
				if(result < upperBound) {
					upperBound = result;
					cout << result << endl;
				}
				
				// restore previous graph	
				if(remNodeIsTerminal) {
					terminals.pushFront(nodeToRemove);
				}
				if(!targetNodeIsTerminal) {
					terminals.del(terminals.search(targetNode));
				}
			      
				// restore moved edges 
				forall_listiterators(node, it, movedEdges) {
					node v = *it;

					edge e = graph.searchEdge(v, targetNode);
					if(e == NULL) {
						e = graph.searchEdge(targetNode, v);
					}
					OGDF_ASSERT(e != NULL);
					OGDF_ASSERT(e->opposite(targetNode) != nodeToRemove);
		
					//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
					//edge f;
					//cout << v << endl;
					if(e->source() == v) {
						//f = 
						graph.moveTarget(e, nodeToRemove);
						//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
					}
					else {
						//f = 
						graph.moveSource(e, nodeToRemove);
						//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
					}
					//for(int i = 0; i <= depth; i++) cout << " ";
					//cout << "moved edge " << f << " to former position" << endl;
				}

				// restored deleted edges
				while(!delEdges.empty()) {
					node source = delEdges.popFrontRet();
					node target = delEdges.popFrontRet();
					edge e = graph.newEdge(source, target);
					OGDF_ASSERT(!origDelEdges.empty());
					mapping[e] = origDelEdges.popFrontRet();
					//for(int i = 0; i <= depth; i++) cout << " ";
					//cout << "restored formerly deleted edge " << e << endl;
				}
				OGDF_ASSERT(origDelEdges.empty());
				//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
				
				// sencond branch: Exclusion of the edge
				double exEdgeResult = bnbInternal(
				  graph,
				  mapping, 
				  origWeights, 
				  terminals, 
				  tree, 
				  upperBound,
				  prevCost, 
				  depth+1);

				// decide which branch returned best result
				if(exEdgeResult < result) {
					result = exEdgeResult;
					if(result < upperBound) {
						cout << result << endl;
					}
				}

				// finally: restore the branching edge
				edge f = graph.newEdge(nodeToRemove, targetNode);
				mapping[f] = origBranchingEdge;
				//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
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
  Graph &tree)
{
	tree.clear();
	
	Graph copiedGraph;
	NodeArray<node> origNodes(copiedGraph);
	NodeArray<node> copiedNodes(graph);
	EdgeArray<edge> origEdges(copiedGraph);
	List<node> copiedTerminals;

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
	}
	
	//EdgeComparer comp(&origEdges, &weights);

	//writeSVG(graph, terminals, "original.svg");
	//writeSVG(copiedGraph, copiedTerminals, "copy.svg");

	return bnbInternal(
	  copiedGraph, 
	  origEdges,
	  weights, 
	  copiedTerminals, 
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
