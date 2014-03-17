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
	void writeSVG(const std::string &filename);

	/**
	 * Used to validate the current mapping of edges to orignal edges
	 * Used solely for debugging.
	 * The mapping is validated by using OGDF_ASSERT.
	 * 
	 * \return
	 * 	always returns true
	 */
	bool validateMapping();
	
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
	
	/**
	 * Returns the cost of the specified edge.
	 * Looks up the corresponding edge in the original graph
	 * and retrieves its weight.
	 * 
	 * \return
	 *   weight of e
	 */
	double weightOf(edge e);

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