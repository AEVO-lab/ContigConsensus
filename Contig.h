#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Nucleotide.h"

using namespace std;

class Contig
{
  friend ostream& operator<<(ostream&, const Contig&);
private:

  vector<Nuc> sequence;
  const string name;
  const unsigned set_id;
  size_t sequence_size;


  struct Component{
    string name;
    unsigned set_id;
    bool is_reversed;
    size_t shift, contig_size;

    Component(size_t contig_size,string name, unsigned set_id, size_t shift,bool is_reversed=false) :  name(name), set_id(set_id), is_reversed(is_reversed), shift(shift), contig_size(contig_size) {}

    //for debug
    void display()
    {
      cout << "{" << name << "," << set_id << "," << is_reversed << "," << shift << "}" << endl;
    }
  };
  vector<Component> component_contigs;

  
public:
  Contig(const Contig & c) : name(c.name), set_id(c.set_id), sequence(c.sequence), sequence_size(c.sequence_size), component_contigs(c.component_contigs)
  {}

  Contig(unsigned set_id, string &&name,vector<Nuc> && sequence) : name(name), set_id(set_id), sequence(move(sequence)) {
    sequence_size=this->sequence.size();
    component_contigs.emplace_back(sequence_size,name,set_id,0);
  }

  Contig(Contig && c, unsigned set_id) : name(c.name), set_id(set_id), sequence(c.sequence), component_contigs(c.component_contigs) {
    //name(move(c.name)), set_id(set_id), sequence(move(c.sequence)), component_contigs(move(c.component_contigs)) {
    sequence_size=this->sequence.size();
  }

    // merge two contigs
  Contig(Contig && c1, Contig && c2, bool is_c1_reversed, bool is_c2_reversed, size_t overlapping_size, unsigned set_id) :  name("("+ c1.name + "|" + c2.name +")"), set_id(set_id), sequence_size(c1.size())
  {

    if(!is_c1_reversed){
      sequence = c1.sequence;//move(c1.sequence);
      component_contigs = move(c1.component_contigs);
    }
    else {
      for(unsigned i=c1.size();i>0;--i)
	sequence.push_back(c1.sequence[i-1]);
      for(auto &C : c1.component_contigs)
	component_contigs.emplace_back(C.contig_size,C.name,C.set_id,c1.sequence_size-1-C.shift,!C.is_reversed);
    }
    

    unsigned shift=c1.size()-overlapping_size;
    if(!is_c2_reversed)
       for(auto &C : c2.component_contigs)
	 component_contigs.emplace_back(C.contig_size,C.name,C.set_id,shift+C.shift,C.is_reversed);
    else
      for(auto &C : c2.component_contigs)
	component_contigs.emplace_back(C.contig_size,C.name,C.set_id,shift+c2.sequence_size-1-C.shift,!C.is_reversed);

    if(is_c2_reversed){
      for(size_t i=0;i<overlapping_size;++i){
	if(c2.size()<=i) return;
	sequence[sequence.size()-overlapping_size+i]+=c2.sequence[c2.sequence.size()-1-i];
      }
      for(size_t i=overlapping_size;i<c2.sequence.size();++i)
	sequence.push_back(c2.sequence[c2.sequence.size()-1-i]);
    }
    else { 
      for(size_t i=0;i<overlapping_size;++i){
	if(c2.size()<=i) return;
	assert(i<c2.sequence.size());
	assert(sequence.size()-overlapping_size+i>=0 && sequence.size()-overlapping_size+i<=sequence.size());
	sequence[sequence.size()-overlapping_size+i]+=c2.sequence[i];
      }
      
      for(size_t i=overlapping_size;i<c2.sequence.size();++i)
	sequence.push_back(move(c2.sequence[i]));
      
    }
    sequence_size=sequence.size();
  }



  unsigned get_set_id() const
  {
    return set_id;
  }
  
  bool operator<(const Contig &c) const
  {
    return name<c.getName();
  }
  

  const string & getName() const
  { 
    return name;
  } 


  char at(size_t i, bool isReverse = false)
  {
    if (!isReverse)
      return sequence[i];
    else
      return sequence[sequence.size() - 1 - i];
  }

  char getNuc(size_t i)
  {
    assert(i<=sequence.size());
    return sequence[i];
  }


  size_t size() const
  {
    return sequence_size;
  }

  void display_sequence() const // for debug
  {
    for(auto &N : this->sequence)
      cout << (char)N;
    cout << endl;
  }
  void display_component() const // for debug
  {
    for(auto c : this->component_contigs)
      c.display();
  }

  vector<Component>& getComponent() {
    return component_contigs;
  }
};



bool operator<(const unique_ptr<Contig> &c, const string &s)
{
  return c->getName() < s;
}

bool operator<(const string &s,const unique_ptr<Contig> &c)
{
  return s< c->getName();
}

bool operator<(const unique_ptr<Contig> &c1, const unique_ptr<Contig> &c2)
{
  return c1->getName()< c2->getName();
}

typedef map<unsigned, set<unique_ptr<Contig>,less<>>> AssemblySet;




inline ostream &operator<<(ostream& os, const Contig& contig)
{
  os << ">"<<contig.name << endl;
  unsigned length=0;
  for(Nuc c : contig.sequence){
    os<<(char)c;
    length++;
    if(length==79){
      length=0;
      os << "\n";
    }
  }
  os << endl;
  return os;
}
