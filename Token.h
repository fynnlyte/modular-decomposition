#ifndef TOKEN_H
#define TOKEN_H

#include <iostream>
#include <string>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include "Commons.h"

enum class Type
{
  LBracket, // "(" 
  Node, // "Real" vertex
  RBracket // ")"
};

class Token
{
public:
  //Factory methods
  static Token createLBracket();
  static Token createRBracket();
  static Token createNode(Vertex v);

  const std::string toString(const Graph& g) const;
  Type getType() const;
  Vertex getVertex() const;
private:
  Token(Type t, Vertex v);

  Type type;
  Vertex vertex;
};

#endif
