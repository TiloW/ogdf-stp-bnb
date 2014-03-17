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

void STPSolver::writeSVG(const std::string &filename) const
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

double STPSolver::weightOf(edge e) const
{
	double result = MAX_WEIGHT;
	
	if(e != NULL) {
		OGDF_ASSERT(m_mapping[e] != NULL);
		OGDF_ASSERT(m_mapping[e]->graphOf() == &m_originalGraph);
		result = m_originalWeights[m_mapping[e]];
	}
	
	return result;
}

bool STPSolver::validateMapping() const
{
	edge e;
	forall_edges(e, m_graph) {
		OGDF_ASSERT(m_mapping[e] != NULL);
		OGDF_ASSERT(m_mapping[e]->graphOf() == &m_originalGraph);
		OGDF_ASSERT(m_originalWeights[m_mapping[e]]);
	}
	return true;
}

edge STPSolver::deleteEdge(edge e)
{
	edge result = m_mapping[e];
	m_graph.delEdge(e);
	e->~EdgeElement();
	return result;
}

edge STPSolver::newEdge(node source, node target, edge e)
{
	edge result = m_graph.newEdge(source, target);
	m_mapping[result] = e;
	return result;
}

void STPSolver::moveSource(edge e, node v)
{
	m_graph.moveSource(e, v);
}

void STPSolver::moveTarget(edge e, node v)
{
	m_graph.moveTarget(e, v);
}

void STPSolver::setTerminal(node v, bool isTerminal)
{
	if(isTerminal && !m_isTerminal[v]) {
		m_terminals.pushFront(v);
		m_isTerminal[v] = true;
	}
	else {
		if(!isTerminal && m_isTerminal[v]) {
			m_terminals.del(m_terminals.search(v));
			m_isTerminal[v] = false;
		}
	}
}

double STPSolver::solve(Graph &tree)
{
	return bnbInternal(0);
}

edge STPSolver::determineBranchingEdge(double prevCost)
{
	edge result = NULL;
	double maxPenalty = -1;

	// calculate penalties for nodes
	double sumOfMinWeights = 0; // b
	double sumOfMinTermWeights = 0; // c
	double absoluteMinTermWeight = MAX_WEIGHT;
	for(ListConstIterator<node> it = m_terminals.begin(); 
	  sumOfMinWeights < MAX_WEIGHT && it.valid(); 
	  ++it) {
		double minTermWeight = MAX_WEIGHT,
		       minWeight = MAX_WEIGHT,
		       secondMinWeight = MAX_WEIGHT;
		edge minEdge = NULL;

		// investigate all edges of each terminal
		// calculate lower boundary and find branching edge
		List<edge> adjEdges;
		m_graph.adjEdges(*it, adjEdges);
		for(ListConstIterator<edge> itEdge = adjEdges.begin(); 
		  itEdge.valid(); 
		  ++itEdge) {
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
			result = minEdge;
			if(minWeight == MAX_WEIGHT) {
				sumOfMinWeights = MAX_WEIGHT;
			}
		} else {
			sumOfMinWeights += minWeight;
			// update branching edge if need be
			double penalty = secondMinWeight - minWeight;
			if(penalty > maxPenalty) {
				maxPenalty = penalty;
				result = minEdge;
			}
		}
	}

	// compare bounds for this graph
	if(result != NULL) {
		double maxCost = m_upperBound - prevCost;
		if(sumOfMinTermWeights < MAX_WEIGHT) {
			sumOfMinTermWeights -= absoluteMinTermWeight;
		}
		if(maxCost <= sumOfMinWeights && maxCost <= sumOfMinTermWeights) {
			result = NULL;
		}
	}

	return result;
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
			edge branchingEdge = determineBranchingEdge(prevCost);
			
			// branching edge has been found or there is no feasible solution
			if(branchingEdge != NULL && weightOf(branchingEdge) < MAX_WEIGHT) {	
				// chose node to remove
				node nodeToRemove = branchingEdge->source();
				node targetNode = branchingEdge->target();
			
				// This seems to speed up things.
				// Always chosing only terminal nodes to remove 
				// (or only steiner nodes if possible)
				// yields even worse running times than no swapping at all.
				if(randomDouble(0, 1) < 0.5) {
					nodeToRemove = branchingEdge->target();
					targetNode = branchingEdge->source();
				}
				
				// remove branching edge
				edge origBranchingEdge = deleteEdge(branchingEdge);

				// first branch: Inclusion of the edge
				// remove source node of edge and calculate new edge weights	
				OGDF_ASSERT(targetNode != NULL);
				OGDF_ASSERT(nodeToRemove != NULL);
				OGDF_ASSERT(m_graph.searchEdge(targetNode, nodeToRemove) == NULL);
				OGDF_ASSERT(m_graph.searchEdge(nodeToRemove, targetNode) == NULL);
				
				List<node> delEdges, movedEdges;
				List<edge> origDelEdges;
				
				List<edge> edges;
				m_graph.adjEdges(nodeToRemove, edges);
				for(ListConstIterator<edge> it = edges.begin(); it.valid();) {
					edge e = *it;
					it = it.succ();

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
							deleteEdge(e);

							deletedEdgeE = true;
						}
						else {
							delEdges.pushFront(f->target());
							delEdges.pushFront(f->source());
							origDelEdges.pushFront(m_mapping[f]);
							deleteEdge(f);
						}
					}
					if(!deletedEdgeE) {	
						if(e->target() == nodeToRemove) {
							OGDF_ASSERT(e->source() != targetNode);
							movedEdges.pushFront(e->source());
							moveTarget(e, targetNode);
						}
						else {
							OGDF_ASSERT(e->source() == nodeToRemove)
							OGDF_ASSERT(e->target() != targetNode)
							movedEdges.pushFront(e->target());
							moveSource(e, targetNode);
						}
					}
				}
				// nodeToRemove is isolated at this point
				// thus no need to actually remove it
				// (easier to keep track of CopyGraph mapping)
				
				// remove node from terminals too
				bool targetNodeIsTerminal = m_isTerminal[targetNode],
				     nodeToRemoveIsTerminal = m_isTerminal[nodeToRemove];

				OGDF_ASSERT(targetNodeIsTerminal || nodeToRemoveIsTerminal);
				setTerminal(nodeToRemove, false);
				setTerminal(targetNode, true);

				// calculate result on modified graph
				result = bnbInternal(weightOf(branchingEdge) + prevCost);
				
				// restore previous graph

				// restore terminals
				setTerminal(nodeToRemove, nodeToRemoveIsTerminal);
				setTerminal(targetNode, targetNodeIsTerminal);	
			      
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
						moveTarget(e, nodeToRemove);
					}
					else {
						moveSource(e, nodeToRemove);
					}
				}

				// restored deleted edges
				while(!delEdges.empty()) {
					OGDF_ASSERT(!origDelEdges.empty());
					
					node source = delEdges.popFrontRet();
					node target = delEdges.popFrontRet();

					newEdge(source, target, origDelEdges.popFrontRet());
				}
				OGDF_ASSERT(origDelEdges.empty());
				
				// sencond branch: Exclusion of the edge
				double exEdgeResult = bnbInternal(prevCost);

				// decide which branch returned best result
				if(exEdgeResult < result) {
					result = exEdgeResult;
				}

				// finally: restore the branching edge
				newEdge(nodeToRemove, targetNode, origBranchingEdge);
			}
		}
		OGDF_ASSERT(validateMapping());
	}
	return result;
}
