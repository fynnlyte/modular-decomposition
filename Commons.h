#ifndef COMMONS_H
#define COMMONS_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>

struct vertex_prop;

struct corrV_t {
  typedef boost::vertex_property_tag kind;
};

struct facPos_t {
  typedef boost::vertex_property_tag kind;
};

struct vmap_t {
  typedef boost::vertex_property_tag kind;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, vertex_prop, boost::vecS> Tree; // Tree, we choose VertexList = listS because some deleting is going on
typedef typename boost::graph_traits<Tree>::vertex_descriptor TVertex;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
			    boost::property<boost::vertex_name_t, std::string, boost::property<corrV_t, TVertex, boost::property<facPos_t, int> > >,
			    boost::property<boost::edge_index_t, int> > Graph; // graph
			   
typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
			   
typedef boost::subgraph< Graph > SubGraph;
typedef typename boost::graph_traits<Tree>::edge_descriptor TEdge;
typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
typedef typename boost::graph_traits<SubGraph>::vertex_descriptor SVertex;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
typedef boost::graph_traits<Graph>::edge_iterator eedge_iter;
typedef boost::graph_traits<Graph>::adjacency_iterator aiter;
typedef boost::graph_traits<SubGraph>::vertex_iterator svertex_iter;
typedef boost::graph_traits<Tree>::vertex_iterator tvertex_iter;
typedef boost::graph_traits<Graph>::out_edge_iterator g_edge_iter;
typedef boost::graph_traits<Tree>::edge_iterator tedge_iter;
typedef boost::graph_traits<Tree>::out_edge_iterator edge_iter; // OUT! JEEZ!
typedef boost::graph_traits<Tree>::adjacency_iterator taiter;
  
#endif