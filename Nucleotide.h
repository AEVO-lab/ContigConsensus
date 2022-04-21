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
    case 'R':
      A = A==4 ? 4 : A+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'Y':
      T = T==4 ? 4 : T+1;
      C = C==4 ? 4 : C+1;
      break;
    case 'K':
      T = T==4 ? 4 : T+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'M':
      A = A==4 ? 4 : A+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'S':
      C = C==4 ? 4 : C+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'W':
      A = A==4 ? 4 : A+1;
      T = T==4 ? 4 : T+1;
      break;
    case 'B':
      T = T==4 ? 4 : T+1;
      C = C==4 ? 4 : C+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'D':
      A = A==4 ? 4 : A+1;
      T = T==4 ? 4 : T+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'H':
      A = A==4 ? 4 : A+1;
      T = T==4 ? 4 : T+1;
      C = C==4 ? 4 : C+1;
      break;
    case 'V':
      A = A==4 ? 4 : A+1;
      C = C==4 ? 4 : C+1;
      G = G==4 ? 4 : G+1;
      break;
    case 'N':
      A = A==4 ? 4 : A+1;
      T = T==4 ? 4 : T+1;
      C = C==4 ? 4 : C+1;
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
    unsigned m = max(max(A,T),max(C,G));
    bool bA = A==m;
    bool bT = T==m;
    bool bC = C==m;
    bool bG = G==m;
    if(bA && bT && bC && bG)
      return 'N';
    if(bA && !bT && !bC && !bG)
      return 'A';

    if(!bA && bT && bC && bG)
      return 'T';

    if(!bA && !bT && bC && !bG)
      return 'C';

    if(!bA && !bT && !bC && bG)
      return 'G';

    if(bA && !bT && !bC && bG)
      return 'R';

    if(!bA && bT && bC && !bG)
      return 'Y';

    if(!bA && bT && bC && !bG)
      return 'K';

    if(bA && !bT && !bC && bG)
      return 'M';

    if(!bA && !bT && bC && bG)
      return 'S';

    if(bA && bT && !bC && !bG)
      return 'W';

    if(!bA && bT && bC && bG)
      return 'B';

    if(bA && bT && !bC && bG)
      return 'D';

    if(bA && bT && bC && !bG)
      return 'H';

    return 'V';
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
