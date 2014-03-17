#include "STPSolver.h"

STPSolver::STPSolver(
  const Graph &graph,
  const EdgeArray<double> &weights, 
  const List<node> &terminals) : 
    m_originalGraph(graph),
    m_originalWeights(weights),
    m_originalTerminals(terminals)
{
	m_upperBound = MAX_WEIGHT;
	m_graph = Graph();
	m_mapping = EdgeArray<edge>(m_graph);
	m_isTerminal = NodeArray<bool>(m_graph, false);
	
	NodeArray<node> copiedNodes(m_originalGraph);

	node v;
	forall_nodes(v, m_originalGraph) {
		node copiedNode = m_graph.newNode();
		copiedNodes[v] = copiedNode;
	}

	edge e;
	forall_edges(e, m_originalGraph) {
		m_mapping[m_graph.newEdge(
		  copiedNodes[e->source()],
		  copiedNodes[e->target()]
		)] = e;
	}

	forall_listiterators(node, it, m_originalTerminals) {
		m_terminals.pushFront(copiedNodes[*it]);
		m_isTerminal[copiedNodes[*it]] = true;
	} 
}

void STPSolver::writeSVG(const std::string &filename)
{
	GraphAttributes attr(m_graph,
	  GraphAttributes::edgeGraphics |
	  GraphAttributes::nodeGraphics |
	  GraphAttributes::nodeLabel |
	  GraphAttributes::nodeStyle);

	node v;
	forall_nodes(v, m_graph) {
		attr.width(v) = attr.height(v) = 20;
		std::stringstream ss;
		ss << v;
		attr.label(v) = ss.str();
		if(m_isTerminal[v]) {
			attr.fillColor(v) = Color::Red;
		}
	}

	FMMMLayout layout;
	layout.useHighLevelOptions(true);
	layout.unitEdgeLength(40);
	layout.newInitialPlacement(true);
	layout.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	layout.call(attr);
	GraphIO::drawSVG(attr, filename);
}

double STPSolver::weightOf(edge e) {
	double result = MAX_WEIGHT;
	
	if(e != NULL) {
		OGDF_ASSERT(m_mapping[e] != NULL);
		OGDF_ASSERT(m_mapping[e]->graphOf() == &m_originalGraph);
		result = m_originalWeights[m_mapping[e]];
	}
	
	return result;
}

bool STPSolver::validateMapping()
{
	edge e;
	forall_edges(e, m_graph) {
		OGDF_ASSERT(m_mapping[e] != NULL);
		OGDF_ASSERT(m_mapping[e]->graphOf() == &m_originalGraph);
		OGDF_ASSERT(m_originalWeights[m_mapping[e]]);
	}
	return true;
}

double STPSolver::solve(Graph &tree)
{
	return bnbInternal(0);
}

double STPSolver::bnbInternal(double prevCost)
{
	double result = MAX_WEIGHT;
	
	if(prevCost < m_upperBound) {
		if(m_terminals.size() < 2) {
			// all terminals are connected
			m_upperBound = prevCost;
			result = prevCost;
		} 
		else {
			edge branchingEdge = NULL;
			double maxPenalty = -1;

			// calculate penalties for nodes
			double sumOfMinWeights = 0; // b
			double sumOfMinTermWeights = 0; // c
			double absoluteMinTermWeight = MAX_WEIGHT;
			for(ListConstIterator<node> it = m_terminals.begin(); sumOfMinWeights < MAX_WEIGHT && it.valid(); ++it) {
				double minTermWeight = MAX_WEIGHT,
				       minWeight = MAX_WEIGHT,
				       secondMinWeight = MAX_WEIGHT;
				edge minEdge = NULL;

				// investigate all edges of each terminal
				// calculate lower boundary and find branching edge
				List<edge> adjEdges;
				m_graph.adjEdges(*it, adjEdges);
				for(ListConstIterator<edge> itEdge = adjEdges.begin(); itEdge.valid(); ++itEdge) {
					edge e = *itEdge;
					if(weightOf(e) < minWeight) {
						secondMinWeight = minWeight;
						minWeight = weightOf(e);
						minEdge = e;
					}
					else {
						if(weightOf(e) < secondMinWeight) {
							secondMinWeight = weightOf(e);
						}
					}
					
					if(m_isTerminal[e->opposite(*it)] && weightOf(e) < minTermWeight) {
						minTermWeight = weightOf(e);
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
				double maxCost = m_upperBound - prevCost;
				if(sumOfMinTermWeights < MAX_WEIGHT) {
					sumOfMinTermWeights -= absoluteMinTermWeight;
				}
				if(maxCost <= sumOfMinWeights && maxCost <= sumOfMinTermWeights) {
					branchingEdge = NULL;
				}
			}

			OGDF_ASSERT(branchingEdge == NULL || m_mapping[branchingEdge] != NULL);
			
			// branching edge has been found or there is no feasible solution
			if(branchingEdge != NULL && weightOf(branchingEdge) < MAX_WEIGHT) {	
				// remove branching edge
				node nodeToRemove = branchingEdge->source();
				node targetNode = branchingEdge->target();
				edge origBranchingEdge = m_mapping[branchingEdge];
				m_graph.delEdge(branchingEdge);
				branchingEdge->~EdgeElement();

				// first branch: Inclusion of the edge
				// remove source node of edge and calculate new edge weights	
				OGDF_ASSERT(targetNode != NULL);
				OGDF_ASSERT(nodeToRemove != NULL);
				OGDF_ASSERT(m_graph.searchEdge(targetNode, nodeToRemove) == NULL);
				OGDF_ASSERT(m_graph.searchEdge(nodeToRemove, targetNode) == NULL);
				
				List<node> delEdges, movedEdges;
				List<edge> origDelEdges;
				ListConstIterator<edge> itNext;
				
				List<edge> edges;
				m_graph.adjEdges(nodeToRemove, edges);
				for(ListConstIterator<edge> it = edges.begin(); it.valid(); it = itNext) {
					edge e = *it;
					itNext = it.succ();

					OGDF_ASSERT(e != branchingEdge);
					OGDF_ASSERT(e->target() == nodeToRemove || e->source() == nodeToRemove);
					OGDF_ASSERT(e->opposite(nodeToRemove) != targetNode);

					node w = e->opposite(nodeToRemove);
					edge f = m_graph.searchEdge(w, targetNode);
					if(f == NULL) {
						f = m_graph.searchEdge(targetNode, w);
					}
					bool deletedEdgeE = false;
					if(f != NULL) {
						if(weightOf(f) < weightOf(e)) {
							delEdges.pushFront(e->target());
							delEdges.pushFront(e->source());
							origDelEdges.pushFront(m_mapping[e]);
							m_graph.delEdge(e);
							e->~EdgeElement();

							deletedEdgeE = true;
						}
						else {
							delEdges.pushFront(f->target());
							delEdges.pushFront(f->source());
							origDelEdges.pushFront(m_mapping[f]);
							m_graph.delEdge(f);
							e->~EdgeElement();
						}
					}
					if(!deletedEdgeE) {	
						if(e->target() == nodeToRemove) {
							OGDF_ASSERT(e->source() != targetNode);
							movedEdges.pushFront(e->source());
							m_graph.moveTarget(e, targetNode);
						}
						else {
							OGDF_ASSERT(e->source() == nodeToRemove)
							OGDF_ASSERT(e->target() != targetNode)
							movedEdges.pushFront(e->target());
							m_graph.moveSource(e, targetNode);
						}
					}
				}
				// nodeToRemove is isolated at this point
				// thus no need to actually remove it
				// (easier to keep track of CopyGraph mapping)
				
				// remove node from terminals too
				bool remNodeIsTerminal = m_isTerminal[nodeToRemove],
                                     targetNodeIsTerminal = m_isTerminal[targetNode];
				if(remNodeIsTerminal) {
					m_terminals.del(m_terminals.search(nodeToRemove));
					m_isTerminal[nodeToRemove] = false;
				}
				if(!targetNodeIsTerminal) {
					m_terminals.pushFront(targetNode);
					m_isTerminal[targetNode] = true;
				}

				// calculate result on modified graph
				result = bnbInternal(weightOf(branchingEdge) + prevCost);
				
				// restore previous graph	
				if(remNodeIsTerminal) {
					m_terminals.pushFront(nodeToRemove);
					m_isTerminal[nodeToRemove] = true;
				}
				if(!targetNodeIsTerminal) {
					m_terminals.del(m_terminals.search(targetNode));
					m_isTerminal[targetNode] = false;
				}
			      
				// restore moved edges 
				while(!movedEdges.empty()) {
					node v = movedEdges.popFrontRet();

					edge e = m_graph.searchEdge(v, targetNode);
					if(e == NULL) {
						e = m_graph.searchEdge(targetNode, v);
					}
					OGDF_ASSERT(e != NULL);
					OGDF_ASSERT(e->opposite(targetNode) != nodeToRemove);
		
					if(e->source() == v) {
						m_graph.moveTarget(e, nodeToRemove);
					}
					else {
						m_graph.moveSource(e, nodeToRemove);
					}
				}

				// restored deleted edges
				while(!delEdges.empty()) {
					node source = delEdges.popFrontRet();
					node target = delEdges.popFrontRet();
					edge e = m_graph.newEdge(source, target);
					OGDF_ASSERT(!origDelEdges.empty());
					m_mapping[e] = origDelEdges.popFrontRet();
				}
				OGDF_ASSERT(origDelEdges.empty());
				
				// sencond branch: Exclusion of the edge
				double exEdgeResult = bnbInternal(prevCost);

				// decide which branch returned best result
				if(exEdgeResult < result) {
					result = exEdgeResult;
				}

				// finally: restore the branching edge
				edge f = m_graph.newEdge(nodeToRemove, targetNode);
				m_mapping[f] = origBranchingEdge;
			}
		}
		//OGDF_ASSERT(validateMapping(graph, mapping, origWeights));
	}
	return result;
}