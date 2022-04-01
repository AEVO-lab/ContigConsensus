#pragma once

#include <array>
#include <bitset>
#include <iostream>

using namespace std;


class Nuc //: //public array<unsigned char,4>
{
  unsigned char A : 2;
  unsigned char T : 2;
  unsigned char C : 2;
  unsigned char G : 2;
  friend Nuc operator+(const Nuc&, const Nuc &);
    
  Nuc() : A(0), T(0), C(0), G(0) {}
  //array<unsigned char, 4>{{0,0,0,0}} {}
    
public:
  Nuc(char c) : Nuc()
  {
    *this=c; 
  }

  void operator=(char c)
  {
    switch (c) {
    case 'N':
      A = A==4 ? 4 : A+1;
      T = T==4 ? 4 : T+1;
      C = C==4 ? 4 : C+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'A':
      A = A==4 ? 4 : A+1;
      break;
    case 'T':
      T = T==4 ? 4 : T+1;
      break;
    case 'C':
      C = C==4 ? 4 : C+1;
      break;
    case 'G':
      G = G==4 ? 4 : G+1;
      break;
    }
  }

  bool operator==(const Nuc &n)
  {
    return A && n.A || T && n.T || C && n.C || G && n.G;
  }

  void operator+=(const Nuc&n)
  {
    A = A + n.A >4 ? 4 : A +n.A;
    T = T + n.T >4 ? 4 : T +n.T;
    C = C + n.C >4 ? 4 : C +n.C;
    G = G + n.G >4 ? 4 : G +n.G;	
  }

  operator char() const
  {

    char m='A';
    if(T>m)
      m='T';
    if(C>m)
      m='C';
    if(G>m)
      m='G';
    return m;
  }

   
};               

Nuc operator+(const Nuc&n1, const Nuc &n2)
{
  Nuc n;
  n.A = n1.A + n2.A >4 ? 4 : n1.A +n2.A;
  n.T = n1.T + n2.T >4 ? 4 : n1.T +n2.T;
  n.C = n1.C + n2.C >4 ? 4 : n1.C +n2.C;
  n.G = n1.G + n2.G >4 ? 4 : n1.G +n2.G;	

  return n;
}
