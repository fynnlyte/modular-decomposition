#ifndef DYCKWORD_H
#define DYCKWORD_H

#include "Token.h"

class DyckWord
{
public:
  DyckWord(const Graph& g, const std::map<Vertex, int>& lbrackets, const std::map<Vertex, int>& rbrackets, const std::vector<Vertex>& fac);
  int size() const;
  Token getToken(int i) const;
  const std::string toString(const Graph& g) const;
private:
  int numBrackets(const std::map<Vertex, int>& brackets, Vertex v);
  std::vector<Token> dyck_word;
};

#endif
