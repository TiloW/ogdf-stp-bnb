#ifndef STP_SOLVER_BNB
#define STP_SOLVER_BNB

#include <limits>
#include <string>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/fileformats/GraphIO.h>

#define MAX_WEIGHT std::numeric_limits<double>::max()

using namespace ogdf;

class STPSolver
{
private:
	const Graph &m_originalGraph;
	const EdgeArray<double> &m_originalWeights;
	const List<node> &m_originalTerminals;
	
	Graph m_graph;
	List<node> m_terminals;
	NodeArray<bool> m_isTerminal;
	EdgeArray<edge> m_mapping;
	
    	double m_upperBound;

	/**
	 * Prints the given Steiner tree problem.
	 * Used solely for debugging.
	 * 
	 * \param filename
	 * 	the name of the SVG file to be created
	 */
	void writeSVG(const std::string &filename) const;

	/**
	 * Used to validate the current mapping of edges to orignal edges
	 * Used solely for debugging.
	 * The mapping is validated by using OGDF_ASSERT.
	 * 
	 * \return
	 * 	always returns true
	 */
	bool validateMapping() const;
		
	/**
	 * Returns the cost of the specified edge.
	 * Looks up the corresponding edge in the original graph
	 * and retrieves its weight.
	 * 
	 * \return
	 *   weight of e
	 */
	double weightOf(edge e) const;

	/**
	 * Calculates the optimal Steinter tree recursivly.
	 * Should not be called directly but by STPSolver::solve.
	 * Each edge is either included or excluded, which gives rise to up to two new branches in each step.
	 * 
	 * \param prevCost
	 *	the cost accumulated in previous recursion steps (previously included edges)
	 * 
	 * \return
	 *	the total weight of the optimal solution (including prevCost)
	 *	note: This might be higher then the actual solution if no solution
	 *	satisfying the upper bound can be found.
	 */
	double bnbInternal(double prevCost);

	/*
	 * Removes the specified edge from the graph.
	 * The corresponding original edge is returned.
	 *
	 * \return
	 *	the original edge, according to m_mapping
	 */
	edge deleteEdge(edge e);

	/*
	 * Creates a new edge.
	 *
	 * \param source
	 *	the source node of the new edge
	 * \param target
	 *	the target node of the new edge
	 * \param originalEdge
	 *	the corresponding edge in the original graph
	 */
	edge newEdge(node source, node target, edge originalEdge);

	/*
	 * Moves the source of the edge to the specified node.
	 *
	 * \param e
	 *	the edge to be moved
	 * \param v
	 * 	the new source of e
	 */
	void moveSource(edge e, node v);

	/*
	 * Moves the target of the edge to the specified node.
	 *
	 * \param e
	 *	the edge to be moved
	 * \param v
	 * 	the new target of e
	 */
	void moveTarget(edge e, node v);

	/*
	 * Updates the status of the given node to
	 * either terminal or steiner node.
	 * No side-effects occur even if status is already correct.
	 *
	 * \param v
	 *	the node to be updated
	 * \param isTerminal
	 * 	true to set it to terminal
	 *	false to set it to steiner
	 */
	void setTerminal(node v, bool isTerminal);

	/**
	 * Decides which edge to branch on.
	 * Might return NULL if current upper bound can not be reached
	 * with this graph.
	 * 
	 * \param prevCost
	 *	the cost of previously chosen edges
	 *	used for comparing to current upper bound
	 *
	 * \return
	 *	edge to branch on next
	 *	might be NULL if current upper bound is not reachable
	 */
	edge determineBranchingEdge(double prevCost);

public:
	/**
	 * Creates an STPSolver for the given STP instance.
	 *
	 * \param graph
	 *	the Graph containing both steiner and terminal nodes
	 * \param weights
	 *	the weight of each edge
	 * \param terminals
	 *	the set of terminal nodes, as opposed to steiner nodes
	 */
	STPSolver(const Graph &graph, const EdgeArray<double> &weights, const List<node> &terminals);
	
	/**
	 * Solves the current STP instance.
	 * Will return the total cost of the optimal solution.
	 * 
	 * \param tree
	 * 	will hold the included edges, not yet implemented
	 *	TODO: Consider EdgeArray<bool>
	 */
	double solve(Graph &tree);
};

#endif
