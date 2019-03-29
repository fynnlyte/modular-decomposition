#include "DyckWord.h"

int DyckWord::numBrackets(const std::map<Vertex, int>& brackets, Vertex v) {
  auto res = brackets.find(v);
  
  if(res == brackets.end()) {
    return 0;
  }
  
  return res->second;
}

DyckWord::DyckWord(const Graph& g, const std::map<Vertex, int>& lbrackets, const std::map<Vertex, int>& rbrackets, const std::vector<Vertex>& fac)
{
  const auto& name = boost::get(boost::vertex_name,g);
  for (const auto& v : fac)
  {
    std::generate_n(std::back_inserter(dyck_word), numBrackets(lbrackets, v), [] { return Token::createLBracket(); });

    dyck_word.push_back(Token::createNode(v)); // vertex

    std::generate_n(std::back_inserter(dyck_word), numBrackets(rbrackets, v), [] { return Token::createRBracket(); });
  }
}

int DyckWord::size() const
{
  return dyck_word.size();
}

Token DyckWord::getToken(int i) const
{
  return dyck_word[i];
}

const std::string DyckWord::toString(const Graph& g) const
{
  std::string ret;
  for (const Token& t : dyck_word)
    ret += t.toString(g);
  return ret;
}