#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <iostream>
#include <random>

using namespace std;


class Nuc //: //public array<unsigned char,4>
{
  // prioritized
  unsigned char A : 4;
  unsigned char T : 4;
  unsigned char C : 4;
  unsigned char G : 4;

  // absorb
  unsigned char a : 4;
  unsigned char t : 4;
  unsigned char c : 4;
  unsigned char g : 4;

  friend Nuc operator+(const Nuc&, const Nuc &);
    
  Nuc() : A(0), T(0), C(0), G(0), a(0), t(0) ,c(0), g(0) {}
  //array<unsigned char, 4>{{0,0,0,0}} {}


  static char convert_char(unsigned char A, unsigned char T, unsigned char C, unsigned char G){
    unsigned m = max(max(A,T),max(C,G));
    bool bA = A==m;
    bool bT = T==m;
    bool bC = C==m;
    bool bG = G==m;
    
    if(bA && bT && bC && bG){
      return 'N';
    }
      
    if(bA && !bT && !bC && !bG)
      return 'A';

    if(!bA && bT && !bC && !bG)
      return 'T';

    if(!bA && !bT && bC && !bG)
      return 'C';

    if(!bA && !bT && !bC && bG)
      return 'G';

    //return 'N';

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
    
public:

  static char Complement(char v)
  {
    switch (v) {
    case 'A':
      return 'T';
      break;
    case 'T':
      return 'A';
      break;
    case 'C':
      return 'G';
      break;
    case 'G':
      return 'C';
      break;
    case 'R':
      return 'Y';
      break;
    case 'Y':
      return 'R';
      break;
    case 'K':
      return 'M';
      break;
    case 'M':
      return 'K';
      break;
    case 'B':
      return 'V';
      break;
    case 'V':
      return 'B';
      break;
    case 'D':
      return 'H';
      break;
    }
    return v;
  }
  
  Nuc(char v, bool prioritized=true) : Nuc()
  {
    // todo correct it
    switch (v) {
    case 'A':
      if(prioritized) A = 1;
      a=1;
      break;
    case 'T':
      if(prioritized) T = 1;
      t=1;
      break;
    case 'C':
      if(prioritized) C = 1;
      c=1;
      break;
    case 'G':
      if(prioritized) G =1;
      g=1;
      break;
    
    // case 'R':
    //   A = A==3 ? 3 : A+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'Y':
    //   T = T==3 ? 3 : T+1;
    //   C = C==3 ? 3 : C+1;
    //   break;
    // case 'K':
    //   T = T==3 ? 3 : T+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'M':
    //   A = A==3 ? 3 : A+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'S':
    //   C = C==3 ? 3 : C+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'W':
    //   A = A==3 ? 3 : A+1;
    //   T = T==3 ? 3 : T+1;
    //   break;
    // case 'B':
    //   T = T==3 ? 3 : T+1;
    //   C = C==3 ? 3 : C+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'D':
    //   A = A==3 ? 3 : A+1;
    //   T = T==3 ? 3 : T+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'H':
    //   A = A==3 ? 3 : A+1;
    //   T = T==3 ? 3 : T+1;
    //   C = C==3 ? 3 : C+1;
    //   break;
    // case 'V':
    //   A = A==3 ? 3 : A+1;
    //   C = C==3 ? 3 : C+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    // case 'N':
    //   A = A==3 ? 3 : A+1;
    //   T = T==3 ? 3 : T+1;
    //   C = C==3 ? 3 : C+1;
    //   G = G==3 ? 3 : G+1;
    //   break;
    }

  }

  Nuc reverse() const
  {
    Nuc n;
    n.A=T;
    n.T=A;
    n.C=G;
    n.G=C;
    n.a=t;
    n.t=a;
    n.c=g;
    n.g=c;
    return n;
  }

  
  bool is_equal(const Nuc &n, bool isReversed=false) const
  {
    if(isReversed)
      return a && n.t || t && n.a || c && n.g || c && n.g;

    return a && n.a || t && n.t || c && n.c || g && n.g;
  }
  
  // bool operator==(const Nuc &n)
  // {
  //   return A && n.A || T && n.T || C && n.C || G && n.G;
  // }

  void operator+=(const Nuc&n)
  {
    A = A + n.A >3 ? 3 : A +n.A;
    T = T + n.T >3 ? 3 : T +n.T;
    C = C + n.C >3 ? 3 : C +n.C;
    G = G + n.G >3 ? 3 : G +n.G;
    
    a = a + n.a >3 ? 3 : a +n.a;
    t = t + n.t >3 ? 3 : t +n.t;
    c = c + n.c >3 ? 3 : c +n.c;
    g = g + n.g >3 ? 3 : g +n.g;	

  }
  

  char disp2(unsigned c) const
  {
    switch (c) {
    case 0 :
      return '0';
    case 1 :
      return '1';
    case 2:
      return '2';
    case 3:
      return '3';
    }
    return '+';
  }
  void display() const
  {
    cout << "(" << disp2(A) << "," << disp2(T) << "," << disp2(C) << "," << disp2(G) << ")";
    //    cout << "(" << disp2(a) << "," << disp2(t) << "," << disp2(c) << "," << disp2(g) << ")";
  }


  
  operator char() const
  {
    
    unsigned m = max(max(A,T),max(C,G));
    bool bA = A==m;
    bool bT = T==m;
    bool bC = C==m;
    bool bG = G==m;
       
    if(bA && !bT && !bC && !bG)
      return 'A';

    if(!bA && bT && !bC && !bG)
      return 'T';

    if(!bA && !bT && bC && !bG)
      return 'C';

    if(!bA && !bT && !bC && bG)
      return 'G';

    if(a*bA+t*bT+c*bC+g*bG==0)
      return convert_char(A,T,C,G);
    return convert_char(a*bA, t*bT, c*bC, g*bG);

    
    // unsigned m = max(max(a,t),max(c,g));
    // bool ba = a==m;
    // bool bt = t==m;
    // bool bc = c==m;
    // bool bg = g==m;
       
    // if(ba && !bt && !bc && !bg)
    //   return 'A';

    // if(!ba && bt && !bc && !bg)
    //   return 'T';

    // if(!ba && !bt && bc && !bg)
    //   return 'C';

    // if(!ba && !bt && !bc && bg)
    //   return 'G';

    // if(A*ba+T*bt+C*bc+G*bg==0)
    //   return convert_char(a,t,c,g);
    // return convert_char(A*ba, T*bt, C*bc, G*bg);
  }

  char output(bool non_random) const
  {

    return (char)(*this);

    
    unsigned m = max(max(A,T),max(C,G));
    bool bA = A==m;
    bool bT = T==m;
    bool bC = C==m;
    bool bG = G==m;
    std::random_device dev;
    std::mt19937 rng(dev());
    string r;
    
    if(bA && bT && bC && bG){
      if(non_random)
	return 'N';
      r ="ATCG";
    }
      
    if(bA && !bT && !bC && !bG)
      return 'A';

    if(!bA && bT && !bC && !bG)
      return 'T';

    if(!bA && !bT && bC && !bG)
      return 'C';

    if(!bA && !bT && !bC && bG)
      return 'G';

    //    return 'N';

    if(bA && !bT && !bC && bG){
      if(non_random)
	return 'R';
      r="AG";
    }

    if(!bA && bT && bC && !bG){
      if(non_random)
	return 'Y';
      r="TC";

    }

    if(!bA && bT && !bC && bG){
      if(non_random)
	return 'K';
      r="TG";
    }
      

    if(bA && !bT && bC && !bG){
      if(non_random)
	return 'M';
      r="AC";
    }

    if(!bA && !bT && bC && bG){
      if(non_random)
	return 'S';
      r="CG";
    }

    if(bA && bT && !bC && !bG){
      if(non_random)
	return 'W';
      r="AT";
    }

    if(!bA && bT && bC && bG){
      if(non_random)
	return 'B';
      r="TCG";
    }

    if(bA && bT && !bC && bG){
      if(non_random)
	return 'D';
      r="ATG";
    }
    if(bA && bT && bC && !bG){
      if(non_random)
	return 'H';
      r="ATC";
    }
    if(non_random)
      return 'V';
    else r="TCG";
    // std::sample(r.begin(), r.end(), std::back_inserter(r),
    //             1, std::mt19937{std::random_device{}()});
    
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(r.begin(), r.end(), g);
    return r[0];
  }

  unsigned score()
  {
    return A+T+C+G;//+a+t+c+g;
  }

   
};               

Nuc operator+(const Nuc&n1, const Nuc &n2)
{
  Nuc n;
  n.A = n1.A + n2.A >3 ? 3 : n1.A +n2.A;
  n.T = n1.T + n2.T >3 ? 3 : n1.T +n2.T;
  n.C = n1.C + n2.C >3 ? 3 : n1.C +n2.C;
  n.G = n1.G + n2.G >3 ? 3 : n1.G +n2.G;	

  n.a = n1.a + n2.a >3 ? 3 : n1.a +n2.a;
  n.t = n1.t + n2.t >3 ? 3 : n1.t +n2.t;
  n.c = n1.c + n2.c >3 ? 3 : n1.c +n2.c;
  n.g = n1.g + n2.g >3 ? 3 : n1.g +n2.g;	

  return n;
}
