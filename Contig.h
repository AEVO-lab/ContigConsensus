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
  friend void write_contig(ostream&, const Contig&,unsigned,bool,bool);
private:

  vector<Nuc> sequence;
  const string name;
  unsigned long contig_id; //id within the set 
  const unsigned set_id;
  size_t sequence_size;
  bool prioritize;

  struct Component{
    string name;
    unsigned set_id;
    bool is_reversed;
    size_t shift, contig_size;

    Component(size_t contig_size,string name, unsigned set_id, size_t shift,  bool is_reversed=false) :  name(name), set_id(set_id), is_reversed(is_reversed), shift(shift), contig_size(contig_size) {}

    //for debug
    void display()
    {
      cout << "{" << name << "," << set_id << "," << is_reversed << "," << shift << "}" << endl;
    }
  };
  vector<Component> component_contigs;

  
public:
  Contig(const Contig & c,unsigned long contig_id, bool prioritize) : name(c.name), contig_id(contig_id),set_id(c.set_id), sequence(c.sequence), sequence_size(c.sequence_size), component_contigs(c.component_contigs), prioritize(prioritize)
  {}//unused

  Contig(unsigned set_id, unsigned long contig_id, string &&name, vector<Nuc> && sequence, bool prioritize) : name(name), contig_id(contig_id), set_id(set_id), sequence(move(sequence)), prioritize(prioritize)
  {
    sequence_size=this->sequence.size();
    component_contigs.emplace_back(sequence_size,name,set_id,0);
  }

  Contig(Contig && c, unsigned set_id, unsigned long contig_id) : name(move(c.name)), contig_id(contig_id),set_id(set_id), sequence(move(c.sequence)), component_contigs(move(c.component_contigs)), prioritize(c.prioritize) {
    sequence_size=this->sequence.size();
  } // unused

    // merge two contigs
  Contig(Contig && c1, Contig && c2, bool is_c1_reversed, bool is_c2_reversed, size_t overlapping_size, unsigned set_id, unsigned long contig_id) :  name("("+ c1.name + "|" + c2.name +")"), contig_id(contig_id), set_id(set_id), sequence_size(c1.size()), prioritize(c1.prioritize || c2.prioritize)
  {
    if(!is_c1_reversed){
      sequence = c1.sequence;//move(c1.sequence);
      component_contigs = move(c1.component_contigs);
    }
    else {
      for(unsigned i=c1.size();i>0;--i)
	sequence.push_back(c1.sequence[i-1].reverse());
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

    sequence_size=sequence.size();
    if(is_c2_reversed){
      for(size_t i=0;i<overlapping_size;++i){
	if(c2.size()<=i) break;
	sequence[sequence.size()-overlapping_size+i]+=c2.sequence[c2.sequence.size()-1-i].reverse();
      }
      for(size_t i=overlapping_size;i<c2.sequence.size();++i)
	sequence.push_back(c2.sequence[c2.sequence.size()-1-i].reverse());
    }
    else {
      for(size_t i=0;i<overlapping_size;++i){
	if(c2.sequence.size()<=i) break;
	assert(i<c2.sequence.size());
	assert(sequence.size()-overlapping_size+i>=0 && sequence.size()-overlapping_size+i<=sequence.size());
	sequence[sequence.size()-overlapping_size+i]+=c2.sequence[i];
      }
      
      for(size_t i=overlapping_size;i<c2.sequence.size();++i)
	sequence.push_back(move(c2.sequence[i]));
      
    }
    c1.sequence.clear();
    c1.sequence.shrink_to_fit();
    c2.sequence.clear();
    c2.sequence.shrink_to_fit();
    sequence_size=sequence.size();
  }


  void absorb(Contig &c, unsigned shift, bool reverse)
  {
    if(c.prioritize)
      prioritize=true;
    
    assert(shift+c.size()<=sequence_size);
    //c.display_sequence();

    if(!reverse){
      for(unsigned i=0; i<c.size(); ++i){	  
	sequence[i+shift]+=c.sequence[i];
      }
    }
    else {

      for(unsigned i=0; i<c.size(); ++i){
	sequence[i+shift]+=c.sequence[c.sequence_size-1-i].reverse();
      }
     
    }
    c.sequence.clear();
    c.sequence.shrink_to_fit();
  }


  unsigned get_set_id() const
  {
    return set_id;
  }

  unsigned long get_contig_id() const
  {
    return contig_id;
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
    
    return sequence[sequence.size() - 1 - i];
  }

  Nuc getNuc(size_t i, bool isReverse = false) const
  {
    if(!(i>=0 && i<sequence.size())){
      cout << "error: " << contig_id << endl;
      cout << "i: " << i << " size: " << sequence_size << "|" << sequence.size() << endl; 
    }
    assert(i>=0 && i<sequence.size());
    if (!isReverse)
      return sequence[i];

    //    return sequence[sequence.size() - 1 - i];
    return (sequence[i]).reverse();

  }


  size_t size() const
  {
    return sequence_size;
  }

  void display_sequence(bool reverse=false) const // for debug
  {
    unsigned count=0;
    if(reverse){
      for(size_t i=sequence_size;i>0;i--){
	cout << (char)getNuc(i-1,reverse);
	count++;
	if(count==100){
	  count=0;
	  cout << endl;
	}
      }
      cout << endl;
      return;
    }

      
    for(auto &N : this->sequence){
      cout << (char)N;
      	count++;
	if(count==100){
	  count=0;
	  cout << endl;
	}

    }
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



// bool operator<(const unique_ptr<Contig> &c, const string &s)
// {
//   return c->getName() < s;
// }

// bool operator<(const string &s,const unique_ptr<Contig> &c)
// {
//   return s< c->getName();
// }
bool operator<(const unique_ptr<Contig> &c, const unsigned long &s)
{
  return c->get_contig_id() < s;
}

bool operator<(const unsigned long &s,const unique_ptr<Contig> &c)
{
  return s< c->get_contig_id();
}
bool operator<(const unique_ptr<Contig> &c1, const unique_ptr<Contig> &c2)
{
  return c1->get_contig_id()< c2->get_contig_id();
}

typedef map<unsigned, set<unique_ptr<Contig>,less<>>> AssemblySet;




inline ostream &operator<<(ostream& os, const Contig& contig)
{
  os << ">" << contig.name << endl;
  unsigned length=0;
  for(Nuc c : contig.sequence){
    //if((char)c=='N') continue;
    os<<(char)c;
    length++;
    if(length==99){
      length=0;
      os << "\n";
    }
  }
  os << endl;
  return os;
}


void write_contig(ostream& os, const Contig& contig,unsigned id, bool random_when_tie, bool prioritize)
{
  if(prioritize && !contig.prioritize) return;
  //os << ">"<< id << "_" << contig.size() << "_" << contig.sequence.size() << endl;// << contig.name << endl;
  os << ">" << contig.name << endl;
  unsigned tmp=0;
  unsigned length=0;
  for(Nuc c : contig.sequence){
    //if((char)c=='N') continue;
    os<<(char)c;
    length++;
    if(length==99){
      length=0;
      os << "\n";
    }
  }
  os << endl;
}

