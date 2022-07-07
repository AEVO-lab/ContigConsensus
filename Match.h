#pragma once
#include <algorithm>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "Contig.h"
#include "CostMap.hpp"


using namespace std;


//Note: this is a struct, so everything is public
struct Match
{
  static bool shorter(const Match *m1,const Match *m2)
  {
    return m1->length(0)+m1->length(1) < m2->length(0)+m2->length(1);
  }

  typedef set<const Match*,decltype(&shorter)> SortedLengthSet;
  typedef map<const Contig*,map<const Contig*,Match*>> MM_map;


  /**
   * Filters matches to return only those that are type 1 or 2 
   **/
  static SortedLengthSet match_filter(unsigned t, vector<Match> & matches)
  {
    SortedLengthSet retval(&Match::shorter);

    for (int i = 0; i < matches.size(); ++i) 
      if (matches[i].is_of_type(t))
	retval.insert(&matches[i]);
    return retval;
  }

  struct Mcontig {
    Contig* contig;
    size_t start, end;
  };

  array<Mcontig, 2> contigs;

  //size_t length;
  unsigned score;

  Match() : score(0) {} // empty match
  Match(Contig *c1, Contig *c2, size_t start1, size_t end1, size_t start2, size_t end2, unsigned score)
    : contigs{{{c1, start1, end1},{c2,start2, end2}}}, score(score) {}

  Match(Contig *c1, Contig *c2, size_t initial_start_c1, bool is_c1_reversed, size_t initial_start_c2, bool is_c2_reversed) : score(0)
  {

    // cout <<"Match with: " << *c1 << " and " << *c2 << endl;
    // cout <<"Parameters: " << initial_start_c1 << "|" << is_c1_reversed << "|" << initial_start_c2 << "|" << is_c2_reversed << endl;
    size_t tail_c1=initial_start_c1, head_c1 = initial_start_c1,
      tail_c2=initial_start_c2, head_c2=initial_start_c2;

    do{
      // cout << "tail_c1:" << tail_c1 << ":"<< (char)c1->getNuc(tail_c1)
      // 	   << " tail_c2:" << tail_c2 << ":"<< (char)c2->getNuc(tail_c2)<< endl;
      score+=c1->getNuc(tail_c1)==c2->getNuc(tail_c2);
      tail_c1 = is_c1_reversed ? tail_c1-1 : tail_c1+1;
      tail_c2 = is_c2_reversed ? tail_c2-1 : tail_c2+1;
    }while(tail_c1<c1->size() && tail_c2<c2->size());

    tail_c1 = is_c1_reversed ? tail_c1+1 : tail_c1-1;
    tail_c2 = is_c2_reversed ? tail_c2+1 : tail_c2-1;


    head_c1 = is_c1_reversed ? head_c1+1 : head_c1-1;
    head_c2 = is_c2_reversed ? head_c2+1 : head_c2-1;
    while(head_c1<c1->size() && head_c2<c2->size()){
      // cout << "head_c1:" << head_c1 << ":"<< (char)c1->getNuc(head_c1)
      // 	   << " head_c2:" << head_c2 << ":"<< (char)c2->getNuc(head_c2)<< endl;

      score+=c1->getNuc(head_c1).is_equal(c2->getNuc(head_c2),is_c1_reversed!=is_c2_reversed);
      head_c1 = is_c1_reversed ? head_c1+1 : head_c1-1;
      head_c2 = is_c2_reversed ? head_c2+1 : head_c2-1;
    }
    head_c1 = is_c1_reversed ? head_c1-1 : head_c1+1;
    head_c2 = is_c2_reversed ? head_c2-1 : head_c2+1;


    unsigned a = c1->get_set_id()<c2->get_set_id() ? 0 : 1;
    unsigned b = c1->get_set_id()<c2->get_set_id() ? 1 : 0;

    contigs[a].contig=c1;
    contigs[a].start=head_c1;
    contigs[a].end=tail_c1;
      
    contigs[b].contig=c2;
    contigs[b].start=head_c2;
    contigs[b].end=tail_c2;

    //      display_contig_names();
  }


  bool is_of_type(unsigned t) const
  {
    return projected_start(t)==0
      || (projected_start((t+1)%2)==0 && projected_end((t+1)%2)==contig((t+1)%2)->size());
  }

  Contig * contig(unsigned t)  const
  {
    return contigs[t].contig;
  }

  //  TODO vérifier que c'est bien le bon match ?
  Contig * contig(const Contig * c) const
  {
    return c==contigs[0].contig ? contigs[1].contig : contigs[0].contig;
  }

  size_t start(unsigned t) const
  {
    return contigs[t].start;
  }

  size_t end(unsigned t) const
  {
    return contigs[t].end;
  }

  size_t length(unsigned t) const
  {
    return (is_reverse(t) ? start(t) - end(t) : end(t) - start(t))+1;
  }

  bool is_reverse(unsigned t) const
  {
    return contigs[t].start > contigs[t].end;
  }

  bool is_full(unsigned t)  const
  {
    return length(t) == (contigs[t].contig->size());
  }

  size_t projected_start(unsigned t) const
  {
    // TODO is that just min(start,end)?
    return min(start(t),end(t));
  }

  size_t projected_end(unsigned t) const
  {
    return max(start(t),end(t));
  }

  // start of the contig if we reverse the direction of the contig
  size_t relative_start(unsigned t) const
  {
    return contigs[t].contig->size()-1-contigs[t].start;
  }

  bool intersect(const Match * m) const
  {
    if(this->contig((unsigned)0)==m->contig((unsigned)0) && this->contig(1)==m->contig(1))
      return true;

    for(size_t i =0; i <=1; i++){
      if (this->contig(i) == m->contig(i) &&
	  ( this->is_reverse(i)!= m->is_reverse(i) ||
	    (this->start(i) >= min(m->start(i),m->end(i))
	     && this->start(i) <= max(m->start(i),m->end(i))) ||
	    (this->end(i) >= min(m->start(i),m->end(i))
	     && this->end(i) <= max(m->start(i),m->end(i))))) {
	return true;
      }
    }
    return false;
  }

  bool intersect_debug(const Match * m) const
  {
    cout << "1): " << (this->contig((unsigned)0)==m->contig((unsigned)0) && this->contig(1)==m->contig(1)) << endl;

    

    
    for(size_t i =0; i <=1; i++){
      cout << "i: " << i  << endl;
      cout << "this->contig(i) == m->contig(i): " << (this->contig(i) == m->contig(i)) << endl;
      cout << "this->is_reverse(i)!= m->is_reverse(i): " << (this->is_reverse(i)!= m->is_reverse(i)) << endl;
      cout << "(this->start(i) >= m->start(i) && this->start(i) <= m->end(i)): " << (this->start(i) >= m->start(i) && this->start(i) <= m->end(i)) << endl;
      cout << "(this->end(i) >= m->start(i) && this->end(i) <= m->end(i)): " << (this->start(i) >= m->start(i) && this->start(i) <= m->end(i)) << endl;
    }
    return false;
  }


  bool intersect(vector<const Match*> v) const
  {
    return any_of(v.begin(), v.end(),
		   [&](const Match *m) {return this->intersect(m);});
  }

  bool contains(const Contig *c) const // for debug
  {
    return contig((unsigned)0)==c || contig(1)==c;
  }

  // for debug
  void display_contig_names() const
  {
    cout << this->contigs[0].contig->getName()
	 << "(" << contigs[0].start <<"," << contigs[0].end<< ", s:"<< contigs[0].start <<")"
	 << "|" << this->contigs[1].contig->getName()
	 << "(" << contigs[1].start <<"," << contigs[1].end<< ", s:"<< contigs[0].start <<")"
	 << "-> " << score
	 << endl;
  }

  void display_alignment() const
  {

    cout << contigs[0].contig->getName() <<"(" << contigs[0].contig->size() << ")"
	 << " vs " << contigs[1].contig->getName() <<"(" << contigs[1].contig->size() << ")" << endl;
    unsigned it1_a=start((unsigned)0),it1_b=start((unsigned)0),it1_c=start((unsigned)0) ;
    unsigned it2_a=start(1),it2_b=start(1),it2_c=start(1);
    unsigned count=0;
    bool cont_a=true, cont_b=true;

    while(cont_a){
      cont_a=false;
      cout << it1_a << "|" << it2_a << "    " << endl;
      
      cont_b=true;
      count=0;
      while (cont_b) {
	count++;
	if(count==80){
	  cont_a=true;
	  cont_b=false;
	}
	cout << contigs[0].contig->getNuc(it1_a,is_reverse((unsigned)0));
	if(it1_a==end((unsigned)0))
	  cont_b=false;
	it1_a = is_reverse((unsigned)0) ? it1_a-1 : it1_a+1;
      }
      std::cout << std::endl;


      cont_b=true;
      count=0;


      while (cont_b) {
	count++;
	if(contigs[0].contig->getNuc(it1_b,is_reverse((unsigned)0))==contigs[1].contig->getNuc(it2_b,is_reverse(1)))
	  cout << "|";
	else cout << " ";
	if(it1_b==end((unsigned)0)||count==80)
	  cont_b=false;
	it1_b = is_reverse((unsigned)0) ? it1_b-1 : it1_b+1;
	it2_b = is_reverse(1) ? it2_b-1 : it2_b+1;
      }
      std::cout << std::endl;
    
      cont_b=true;
      count=0;
      while (cont_b) {
	count++;
	cout << contigs[1].contig->getNuc(it2_a,is_reverse(1));
	if(it2_a==end(1)||count==80)
	  cont_b=false;
	it2_a = is_reverse(1) ? it2_a-1 : it2_a+1;
      }
      std::cout << std::endl;
      std::cout << std::endl;
    

      
    }
  }

};
typedef map<unsigned, map<unsigned,tuple<vector<Match>,unsigned,vector<const Match*>>>> MatchMatrix;


// ostream &operator<<(ostream& os, const Match& m)
// {
//   auto f = max((size_t)22,max(m.start(1)+m.contig(0)->size(),m.start(0)+m.contig(1)->size()));
//   std::cout << endl << string(f,'-')<<endl;
//   cout << /*"Len: " <<m.length<<*/ " - Score: "<< m.score <<"\n";
//   std::cout << string(f,'=')<<endl;

//   for(int c : {0,1}){
//     if(m.is_reverse(c)) cout << "R| ";
//     else cout << " | ";
//     std::cout << string(m.start((c+1)%2),' ');
//     for (int i = 0; i < m.contig(c)->size(); i++) {
//       if (i >= m.start(c) && i < m.start(c) + m.length) {
// 	if(m.contig(c)->at(i)==m.contig((c+1)%2)->at(i+m.start((c+1)%2)-m.start(c)))
// 	  cout << "\033[1;32m";
//         else cout << "\033[1;31m";
//       }

//       cout << (char)m.contig(c)->at(i, m.is_reverse(c));
//       cout << "\033[0m";
//     }
//     cout << endl;
//   }
//   std::cout << string(f,'=')<<endl;
//   return os;
// }






