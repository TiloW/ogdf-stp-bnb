#include "STPSolver.h"

// tmp
#include <iostream>

STPSolver::STPSolver(
  const Graph &graph,
  const EdgeArray<double> &weights, 
  const List<node> &terminals) : 
    m_originalGraph(graph),
    m_originalWeights(weights),
    m_originalTerminals(terminals)
{
	m_upperBound = MAX_WEIGHT;
	m_isTerminal.init(m_originalGraph, NULL);
	int nodeCount = m_originalGraph.numberOfNodes();
	m_edges = Array2D<edge>(0, nodeCount, 0, nodeCount, NULL);

	edge e;
	forall_edges(e, m_originalGraph) {
		setEdge(e->source(), e->target(), e);
	}

	forall_listiterators(node, it, m_originalTerminals) {
		setTerminal(*it, true);
	} 
}
/*
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
		if(isTerminal(v)) {
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
*/

double STPSolver::weightOf(const node u, const node v) const
{
	double result = MAX_WEIGHT;

	edge e = lookupEdge(u, v);
	if(e != NULL) {
		OGDF_ASSERT(e->graphOf() == &m_originalGraph);
		result = m_originalWeights[e];
	}
	
	return result;
}

double STPSolver::weightOf(const NodeTuple &e) const
{
	return weightOf(e.x1(), e.x2());
}

edge STPSolver::lookupEdge(const node u, const node v) const
{
	edge result = NULL;
	if(u != NULL && v != NULL) {
		OGDF_ASSERT(u->graphOf() == &m_originalGraph);
		OGDF_ASSERT(v->graphOf() == &m_originalGraph);
		result = m_edges(u->index(), v->index());
	}
	return result;
}

void STPSolver::setEdge(const node u, const node v, const edge e)
{
	m_edges(u->index(), v->index()) = 
	  m_edges(v->index(), u->index()) = e;
}

void STPSolver::delEdge(const node u, const node v)
{
	setEdge(u, v, NULL);
}

void STPSolver::moveEdge(node u, node v, node w)
{
	OGDF_ASSERT(u != w);
	OGDF_ASSERT(lookupEdge(u, v) != NULL);
	OGDF_ASSERT(lookupEdge(v, w) == NULL);
	setEdge(w, v, lookupEdge(u, v));
	setEdge(u, v, NULL);
}

void STPSolver::setTerminal(const node v, bool makeTerminal)
{
	if(makeTerminal && !isTerminal(v)){
		m_isTerminal[v] = m_terminals.pushFront(v);
	}
	else {
		if(!makeTerminal && isTerminal(v)) {
			m_terminals.del(m_isTerminal[v]);
			m_isTerminal[v] = NULL;
		}
	}
}

bool STPSolver::isTerminal(const node v) const
{
	return m_isTerminal[v].valid();
}

double STPSolver::solve(Graph &tree)
{
	return bnbInternal(0);
}

STPSolver::NodeTuple STPSolver::determineBranchingEdge(double prevCost) const
{
	NodeTuple result = NO_EDGE;
	double maxPenalty = -1;

	// calculate penalties for nodes
	double sumOfMinWeights = 0; // b
	double sumOfMinTermWeights = 0; // c
	double absoluteMinTermWeight = MAX_WEIGHT;
	for(ListConstIterator<node> it = m_terminals.begin(); 
	  sumOfMinWeights < MAX_WEIGHT && it.valid(); 
	  ++it) {
		const node t = *it;
		double minTermWeight = MAX_WEIGHT,
		       minWeight = MAX_WEIGHT,
		       secondMinWeight = MAX_WEIGHT;
		NodeTuple minEdge = NO_EDGE;

		// investigate all edges of each terminal
		// calculate lower boundary and find branching edge
		node v;
		forall_nodes(v, m_originalGraph) {
			if(weightOf(t, v) < minWeight) {
				secondMinWeight = minWeight;
				minWeight = weightOf(t, v);
				minEdge = NodeTuple(t, v);
			}
			else {
				if(weightOf(t, v) < secondMinWeight) {
					secondMinWeight = weightOf(t, v);
				}
			}
			
			if(isTerminal(v) && weightOf(t, v) < minTermWeight) {
				minTermWeight = weightOf(t, v);
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
	if(result != NO_EDGE) {
		double maxCost = m_upperBound - prevCost;
		if(sumOfMinTermWeights < MAX_WEIGHT) {
			sumOfMinTermWeights -= absoluteMinTermWeight;
		}
		if(maxCost <= sumOfMinWeights && maxCost <= sumOfMinTermWeights) {
			result = NO_EDGE;
		}
	}

	return result;
}

double STPSolver::bnbInternal(double prevCost)
{
	double result = MAX_WEIGHT;

	Array2D<edge> tmpEdges(m_edges);
	
	if(prevCost < m_upperBound) {
		if(m_terminals.size() < 2) {
			// all terminals are connected
			m_upperBound = prevCost;
			result = prevCost;
		} 
		else {
			NodeTuple branchingEdge = determineBranchingEdge(prevCost);

			// branching edge has been found or there is no feasible solution
			if(weightOf(branchingEdge) < MAX_WEIGHT) {	
				// chose node to remove
				node nodeToRemove = branchingEdge.x1();
				node targetNode = branchingEdge.x2();
			/*
				// This seems to speed up things.
				if(nodeToRemove->degree() < targetNode->degree()) {
					nodeToRemove = branchingEdge.x2();
					targetNode = branchingEdge.x1();
				}
			*/	
				// remove branching edge
				edge origBranchingEdge = lookupEdge(nodeToRemove, targetNode);
				OGDF_ASSERT(origBranchingEdge != NULL);
				delEdge(nodeToRemove, targetNode);

				// first branch: Inclusion of the edge
				// remove source node of edge and calculate new edge weights	
				OGDF_ASSERT(targetNode != NULL);
				OGDF_ASSERT(nodeToRemove != NULL);
				OGDF_ASSERT(lookupEdge(targetNode, nodeToRemove) == NULL);
				
				List<node> delEdges, movedEdges;
				List<edge> origDelEdges;
				
				node v;
				forall_nodes(v, m_originalGraph) {
					// TODO: Only consider non-isolated nodes
					if(weightOf(v, nodeToRemove) < MAX_WEIGHT) {
						bool doMoveEdge = true;
						OGDF_ASSERT(v != nodeToRemove);
						if(weightOf(v, targetNode) < weightOf(v, nodeToRemove)) {
							delEdges.pushFront(v);
							delEdges.pushFront(nodeToRemove);
							origDelEdges.pushFront(lookupEdge(v, nodeToRemove));
							delEdge(v, nodeToRemove);

							doMoveEdge = false;
						}
						else {
							delEdges.pushFront(v);
							delEdges.pushFront(targetNode);
							origDelEdges.pushFront(lookupEdge(v, targetNode));
							delEdge(v, targetNode);
						}
		
						if(doMoveEdge) {
							movedEdges.pushFront(v);
							moveEdge(nodeToRemove, v, targetNode);	
						}
					}
				}
				// nodeToRemove is isolated at this point
				// thus no need to actually remove it
				// (easier to keep track of CopyGraph mapping)
				
				// remove node from terminals too
				bool targetNodeIsTerminal = isTerminal(targetNode),
				     nodeToRemoveIsTerminal = isTerminal(nodeToRemove);

				OGDF_ASSERT(targetNodeIsTerminal || nodeToRemoveIsTerminal);
				setTerminal(nodeToRemove, false);
				setTerminal(targetNode, true);

				// calculate result on modified graph
				result = bnbInternal(m_originalWeights[origBranchingEdge] + prevCost);
				
				// restore previous graph

				// restore terminals
				setTerminal(nodeToRemove, nodeToRemoveIsTerminal);
				setTerminal(targetNode, targetNodeIsTerminal);	
			      
				// restore moved edges 
				while(!movedEdges.empty()) {
					node v = movedEdges.popFrontRet();
					moveEdge(targetNode, v, nodeToRemove);
				}

				// restore deleted edges
				while(!delEdges.empty()) {
					OGDF_ASSERT(!origDelEdges.empty());
					
					node u = delEdges.popFrontRet();
					node v = delEdges.popFrontRet();

					setEdge(u, v, origDelEdges.popFrontRet());
				}
				OGDF_ASSERT(origDelEdges.empty());
				
				// sencond branch: Exclusion of the edge
				double exEdgeResult = bnbInternal(prevCost);

				// decide which branch returned best result
				if(exEdgeResult < result) {
					result = exEdgeResult;
				}

				// finally: restore the branching edge
				setEdge(nodeToRemove, targetNode, origBranchingEdge);
			}
		}
		//OGDF_ASSERT(validateMapping());
	}

	int nodeCount = m_originalGraph.numberOfNodes();
	for(int i = 0; i < nodeCount; i++)
		for(int ii = 0; ii < nodeCount; ii++)
			OGDF_ASSERT(tmpEdges(i, ii) == m_edges(i, ii));

	return result;
}
