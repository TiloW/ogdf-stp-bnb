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
	m_mapping.init(m_graph);
	m_isTerminal.init(m_graph, NULL);
	int nodeCount = m_originalGraph.numberOfNodes();
	m_edges = Array2D<edge>(0, nodeCount, 0, nodeCount, NULL);
	
	NodeArray<node> copiedNodes(m_originalGraph);

	node v;
	forall_nodes(v, m_originalGraph) {
		node copiedNode = m_graph.newNode();
		copiedNodes[v] = copiedNode;
	}

	edge e;
	forall_edges(e, m_originalGraph) {
		node source = copiedNodes[e->source()],
		     target = copiedNodes[e->target()];
		
		newEdge(source, target, e);
	}

	forall_listiterators(node, it, m_originalTerminals) {
		setTerminal(copiedNodes[*it], true);
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

double STPSolver::weightOf(const edge e) const
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

edge STPSolver::lookupEdge(const node u, const node v) const
{
	return m_edges(u->index(), v->index());
}

void STPSolver::setEdgeLookup(const node u, const node v, const edge e)
{
	m_edges(u->index(), v->index()) = 
	  m_edges(v->index(), u->index()) = e;
}

edge STPSolver::deleteEdge(edge e)
{
	edge result = m_mapping[e];
	setEdgeLookup(e->source(), e->target(), NULL);
	m_graph.delEdge(e);
	e->~EdgeElement();
	
	return result;
}

edge STPSolver::newEdge(node source, node target, edge e)
{
	edge result = m_graph.newEdge(source, target);
	m_mapping[result] = e;
	setEdgeLookup(source, target, result);

	return result;
}

void STPSolver::moveSource(edge e, node v)
{
	OGDF_ASSERT(e != NULL);
	setEdgeLookup(e->source(), e->target(), NULL);
	setEdgeLookup(v, e->target(), e);
	m_graph.moveSource(e, v);
}

void STPSolver::moveTarget(edge e, node v)
{
	OGDF_ASSERT(e != NULL);
	setEdgeLookup(e->source(), e->target(), NULL);
	setEdgeLookup(e->source(), v, e);
	m_graph.moveTarget(e, v);
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

edge STPSolver::determineBranchingEdge(double prevCost) const
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
		const node t = *it;
		double minTermWeight = MAX_WEIGHT,
		       minWeight = MAX_WEIGHT,
		       secondMinWeight = MAX_WEIGHT;
		edge minEdge = NULL;

		// investigate all edges of each terminal
		// calculate lower boundary and find branching edge
		for(adjEntry adj = t->firstAdj(); adj; adj = adj->succ())  {
			edge e = adj->theEdge();
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
			
			if(isTerminal(adj->twinNode()) && weightOf(e) < minTermWeight) {
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
				if(nodeToRemove->degree() < targetNode->degree()) {
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
				
				adjEntry adjNext;
				for(adjEntry adj = nodeToRemove->firstAdj(); adj; adj = adjNext) {
					adjNext = adj->succ();
					edge e = adj->theEdge();

					OGDF_ASSERT(e != branchingEdge);
					OGDF_ASSERT(e->target() == nodeToRemove || e->source() == nodeToRemove);
					OGDF_ASSERT(adj->twinNode() != targetNode);

					edge f = lookupEdge(targetNode, adj->twinNode());
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
				bool targetNodeIsTerminal = isTerminal(targetNode),
				     nodeToRemoveIsTerminal = isTerminal(nodeToRemove);

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

				// restore deleted edges
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
