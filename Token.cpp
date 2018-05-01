#include "Token.h"

Token Token::createLBracket()
{
  Vertex v;
  return Token(Type::LBracket, v);
}
Token Token::createRBracket()
{
  Vertex v;
  return Token(Type::RBracket, v);
}
Token Token::createNode(Vertex v)
{
  return Token(Type::Node, v);
}

Type Token::getType() const
{
  return type;
}

const std::string Token::toString(const Graph& g) const
{
  const auto& name = boost::get(boost::vertex_name, g);
  switch (getType())
  {
    case Type::LBracket : return "(";
    case Type::RBracket : return ")";
    case Type::Node : return boost::get(name, getVertex());
  }
}

Vertex Token::getVertex() const
{
  return vertex;
}

Token::Token(Type t, Vertex v)
{
  vertex = v;
  type = t;
}
