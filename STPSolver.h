#ifndef STP_SOLVER_BNB
#define STP_SOLVER_BNB

#include <limits>
#include <string>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/tuples.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/fileformats/GraphIO.h>

#define MAX_WEIGHT std::numeric_limits<double>::max()

using namespace ogdf;

class STPSolver
{
	typedef Tuple2<node,node> NodeTuple;
private:
	const NodeTuple NO_EDGE = NodeTuple(NULL, NULL);
	const Graph &m_originalGraph;
	const EdgeArray<double> &m_originalWeights;
	const List<node> &m_originalTerminals;
	
	List<node> m_terminals;
	NodeArray<ListIterator<node>> m_isTerminal;
	
    	double m_upperBound;

	Array2D<edge> m_edges;

	/**
	 * Prints the given Steiner tree problem.
	 * Used solely for debugging.
	 * 
	 * \param filename
	 * 	the name of the SVG file to be created
	 */
	//void writeSVG(const std::string &filename) const;

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
	 * \param u
	 * 	one of the nodes incident to the edge to weight
	 * \param v
	 * 	the twin node of u
	 * 
	 * \return
	 *   weight of e
	 */
	double weightOf(const node u, const node v) const;
	
	/**
	 * \param e
	 * 	the edge as a tuple of two nodes
	 */
	double weightOf(const NodeTuple &e) const;
	
	/**
	 * \param e
	 *	the edge in the original graph
	 */
	double weightOf(const edge e) const;

	/**
	 * Calculates the optimal Steinter tree recursivly.
	 * Should not be called directly but by STPSolver::solve.
	 * Each edge is either included or excluded, 
	 * which gives rise to up to two new branches in each step.
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
	 * Removes the edge.
	 * 
	 * \param u
	 * 	one of the nodes incident to the edge to weight
	 * \param v
	 * 	the twin node of u
	 */
	void delEdge(const node u, const node v);

	/**
	 * Sets the edge incident to both node u and v to e.
	 *
	 * \param u
	 * 	one of the nodes incident to the edge to weight
	 * \param v
	 * 	the twin node of u
	 * \param e
	 *	the edge in the original graph
	 *	with associated weight
	 *	set to NULL to remove the edge
	 */
	void setEdge(const node u, const node v, const edge e);

	/**
	 * Move the edge adjacent to node u and v from u to w.
	 * There must be no edge adjacent to both node v and w. 
	 *
	 * \param u
	 * 	one of the nodes incident to the edge to weight
	 * \param v
	 * 	the twin node of u
	 * \param w
	 *	the node replacing u
	 */
	void moveEdge(const node u, const node v, const node w);

	/**
	 * Updates the status of the given node to
	 * either terminal or steiner node.
	 * No side-effects occur even if status is already correct.
	 *
	 * \param v
	 *	the node to be updated
	 * \param makeTerminal
	 * 	true to set it to terminal
	 *	false to set it to steiner
	 */
	void setTerminal(const node v, bool makeTerminal);

	/**
	 * Returns whether this node is a terminal or
	 * a steiner node.
	 *
	 * \param v
	 *	the node to check
	 *
	 * \return
	 *	true if v is a terminal
	 * 	false otherwise
	 *
	 */
	bool isTerminal(const node v) const;

	/**
	 * Retrieves the edge incident to both node u and v.
	 *
	 * \param u
	 *	one of the nodes of the undirected edge
	 * \param v
	 *	the opposite node of u
	 *
	 * \return
	 * 	the edge between u and v
	 *	NULL if it does not exist
	 */ 
	edge lookupEdge(const node u, const node v) const;

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
	NodeTuple determineBranchingEdge(double prevCost) const;

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
	 */
	double solve(Graph &tree);
};

#endif
