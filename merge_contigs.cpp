#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <time.h>
#include <tuple>
#include <unistd.h>
#include <utility>

#include "approximation.cpp"

#define ABSORBED 0

using namespace std;

struct Data
{
  set<unique_ptr<Match>,less<>> all_matches;
  set<unique_ptr<Contig>,less<>> contigs;
  
  AssemblySet assembly_sets;
  MatchMatrix matches;
  map<string,unsigned> ids;
  float percentage;
  unsigned long contig_id_count;

  Data(float percentage=100.0) : percentage(percentage), contig_id_count(1) {} 
};   


void update_match(Data &data, Match &m, pair<size_t,bool> shift, Contig * old_contig, Contig * new_contig, vector<unique_ptr<Match>> &insert)
{
  size_t shift_c1=0, shift_c2=0;
  bool is_c1_reversed=false, is_c2_reversed=false;
  
  Contig * c1 = m.contig((unsigned)0);
  unsigned long id_c1 = c1->get_contig_id();

  Contig * c2 = m.contig(1);
  unsigned long id_c2=c2->get_contig_id();

  if(c1==old_contig){
    c1=new_contig;
    is_c1_reversed=shift.second;
    shift_c1=shift.first;
  } else {
    c2 = new_contig;
    is_c2_reversed=shift.second;
    shift_c2=shift.first;
  }
  
  unsigned start_c1= is_c1_reversed ? shift_c1 + m.relative_start((unsigned)0) : shift_c1 + m.start((unsigned)0);
  unsigned start_c2= is_c2_reversed ? shift_c2 + m.relative_start(1) : shift_c2 + m.start(1);

  unique_ptr<Match> tmp = make_unique<Match>(c1,c2,start_c1,is_c1_reversed!=m.is_reverse((unsigned)0),
				       start_c2,is_c2_reversed!=m.is_reverse(1));
  if(tmp->score*100/(float)tmp->length(0) >=data.percentage)
    insert.push_back(move(tmp));

}



void update_match(Data &data, Match &m, map<unsigned long,tuple<unsigned long,size_t,bool>>& m_contig, unsigned other_set_id, unsigned new_set_id)
{
  auto & new_set=data.assembly_sets[new_set_id];
  size_t shift_c1=0, shift_c2=0;
  bool is_c1_reversed=false , is_c2_reversed=false;
  
  Contig * c1 = m.contig((unsigned)0);
  unsigned long id_c1 = c1->get_contig_id();

  while(m_contig.find(id_c1)!=m_contig.end()){
    auto &r = m_contig[id_c1];
    id_c1=get<0>(r);

    if(id_c1==ABSORBED) return; // absorbed in full match
    shift_c1+=get<1>(r);
    is_c1_reversed=(is_c1_reversed!=get<2>(r));
  }

  if(id_c1!=c1->get_contig_id()){
    assert(new_set.find(id_c1)!=new_set.end());
    c1=new_set.find(id_c1)->get();
  }
  
  unsigned start_c1= is_c1_reversed ? shift_c1 + m.relative_start((unsigned)0) : shift_c1 + m.start((unsigned)0);
  
  Contig * c2 = m.contig(1);
  unsigned long id_c2=c2->get_contig_id();

  while(m_contig.find(id_c2)!=m_contig.end()){
    auto &r = m_contig[id_c2];
    id_c2= get<0>(r);

    if(id_c2==ABSORBED) return; // absorbed in a full match
    shift_c2+=get<1>(r);
    is_c2_reversed=(is_c2_reversed!=get<2>(r));
  }

  if(id_c2!=c2->get_contig_id()){
    assert(new_set.find(id_c2)!=new_set.end());
    c2=new_set.find(id_c2)->get();
  }

  unsigned start_c2= is_c2_reversed ? shift_c2 + m.relative_start(1) : shift_c2 + m.start(1);

  get<0>(data.matches[other_set_id][new_set_id]).emplace_back(Match(c1,c2,start_c1,is_c1_reversed!=m.is_reverse((unsigned)0),
							       start_c2,is_c2_reversed!=m.is_reverse(1)));

  auto & back =get<0>(data.matches[other_set_id][new_set_id]).back();

  
  if(back.score*100.0/(float)back.length((unsigned)0) < data.percentage)
    get<0>(data.matches[other_set_id][new_set_id]).pop_back();
}


Contig* merge_full_match(Data &data, Match *m,map<unsigned long,tuple<unsigned long,size_t,bool>> &m_contig)
{
   bool reverse=m->is_reverse(0)!=m->is_reverse(1);
   size_t shift;
   Contig *c1, *c2;
   if(m->is_full(1)){
     c1 = m->contig((unsigned)0);
      c2 = m->contig(1);
      shift=m->projected_start(0);
    } else {
      c1 = m->contig(1);
      c2 = m->contig((unsigned)0);
      shift=m->projected_start(1);
    } 
   m_contig[c2->get_contig_id()]=make_tuple(ABSORBED,0,false);
    
    return c2;
}

vector<const Match*> merge_full_matches(Data &data,vector<const Match *> &selected_matches, vector<string> &trash,
					map<unsigned long,tuple<unsigned long,size_t,bool>> &m_contig, set<unsigned long> &absorbent,
					unsigned new_set_id, unsigned set_id1, unsigned set_id2)
{
  auto &new_set =data.assembly_sets[new_set_id];

  vector<const Match*> other_matches;
  for(auto &m : selected_matches){
    if(!m->is_full((unsigned) 0) && !m->is_full(1)){
      other_matches.push_back(m);
      continue;
    }
    Contig * c1 = m->contig((unsigned)0);
    unsigned long id_c1 = c1->get_contig_id();
    
    Contig * c2 = m->contig((unsigned)1);
    unsigned long id_c2 = c2->get_contig_id();
 

    unsigned overlap;
    size_t shift;

    bool reverse=m->is_reverse(0)!=m->is_reverse(1);
    if(!m->is_full(1)){
      overlap = m->contig(1)->size()-m->projected_start(1);//m->is_reverse(1) ?  m->contig(1)->size()-m->relative_start(1) : m->contig(1)->size()-m->start(1);
      swap(c1,c2);
      swap(id_c1,id_c2);
      shift=m->projected_start(1);
    } else {
      shift = m->projected_start(0);
      overlap = m->contig((unsigned)0)->size()-m->projected_start(0);
    } 


    c1->absorb(*c2, shift, reverse);
    m_contig[id_c2]=make_tuple(ABSORBED,0,false);
     
  }

  return other_matches;

  
}

Contig* merge_non_full_match(Data &data, Match *m, pair<size_t,bool> &shift_c1, pair<size_t,bool> &shift_c2)
{
  Contig * c1 = m->contig((unsigned)0);
  bool is_c1_reversed=m->is_reverse((unsigned)0);
  unsigned long id_c1=c1->get_contig_id();

  Contig * c2 = m->contig(1);
  bool is_c2_reversed=m->is_reverse(1);
  unsigned long id_c2=c2->get_contig_id();


  unsigned relative_start_c1=is_c1_reversed ? c1->size()-1-m->projected_end((unsigned)0) :
    is_c1_reversed!=m->is_reverse((unsigned)0)? m->contig((unsigned)0)->size()-m->start((unsigned)0) : m->start((unsigned)0);
  unsigned relative_start_c2=is_c2_reversed ? c2->size()-1-m->projected_end(1) :
    is_c2_reversed!=m->is_reverse(1)? m->contig(1)->size()-m->start(1) : m->start(1);//
  m->start(1);


  if(relative_start_c1>relative_start_c2) {
    Contig * s = data.contigs.emplace(move(make_unique<Contig>(move(*c1),move(*c2),
							       is_c1_reversed,
							       is_c2_reversed,
							       m->length(1),0, data.contig_id_count++))).first->get();
    shift_c1 = make_pair(0, is_c1_reversed);
    shift_c2 = make_pair(s->size()-c2->size(), is_c2_reversed);
    return s;
		      
     

  } 

  Contig * s = data.contigs.emplace(move(make_unique<Contig>(move(*c2),move(*c1),
							  is_c2_reversed, is_c1_reversed,
							  m->length(1),
							  0, data.contig_id_count++))).first->get();

  shift_c1=make_pair(s->size()-c1->size(), is_c1_reversed);
  shift_c2=make_pair(0, is_c2_reversed);
  return s;



}

void merge_other_matches(Data &data, vector<const Match *> &selected_matches, vector<string> &trash,
			 map<unsigned long,tuple<unsigned long,size_t, bool>>& m_contig,//map<Contig*,tuple<Contig*,size_t, bool>>& m_contig,
			 unsigned new_set_id,unsigned set_id1,
			 unsigned set_id2)
{

  unsigned it=1;
  auto &new_set =data.assembly_sets[new_set_id];

  for(auto &m : selected_matches){    
    Contig * c1 = m->contig((unsigned)0);
    it++;
    
    unsigned shift_c1=0;
    bool is_c1_reversed=m->is_reverse((unsigned)0);
    unsigned long id_c1=c1->get_contig_id();
    
    while(m_contig.find(id_c1)!=m_contig.end()){
      auto &tmp = m_contig[id_c1];
      id_c1=get<0>(tmp);
      shift_c1+=get<1>(tmp);
      is_c1_reversed=(is_c1_reversed!=get<2>(tmp));
    }
    if(id_c1==0) continue;
    if(c1->get_contig_id()!=id_c1){
      assert(new_set.find(id_c1)!=new_set.end());
      c1=new_set.find(id_c1)->get();
    }

    Contig * c2 = m->contig(1);
    unsigned shift_c2=0;
    bool is_c2_reversed=m->is_reverse(1);
    unsigned long id_c2=c2->get_contig_id();
    while(m_contig.find(id_c2)!=m_contig.end()){
      auto &tmp= m_contig[id_c2];
      id_c2= get<0>(tmp);
      shift_c2+=get<1>(tmp);
      is_c2_reversed=(is_c2_reversed!=get<2>(tmp));
    }
    if(id_c2==0) continue;
    if(c2->get_contig_id()!=id_c2){
      assert(new_set.find(id_c2)!=new_set.end());
      c2=new_set.find(id_c2)->get();
    }


    unsigned relative_start_c1=is_c1_reversed ? c1->size()-1-m->projected_end((unsigned)0) :
      is_c1_reversed!=m->is_reverse((unsigned)0)? m->contig((unsigned)0)->size()-m->start((unsigned)0) : m->start((unsigned)0);
    unsigned relative_start_c2=is_c2_reversed ? c2->size()-1-m->projected_end(1) :
      is_c2_reversed!=m->is_reverse(1)? m->contig(1)->size()-m->start(1) : m->start(1);//
    m->start(1);


    if(relative_start_c1>relative_start_c2) {
      Contig * s = new_set.emplace(move(make_unique<Contig>(move(*c1),move(*c2),
							    is_c1_reversed,
							    is_c2_reversed,
							    m->length(1),new_set_id, data.contig_id_count++))).first->get();
     
      m_contig[c1->get_contig_id()]=make_tuple(s->get_contig_id(), 0, is_c1_reversed);
      m_contig[c2->get_contig_id()]=make_tuple(s->get_contig_id(), s->size()-c2->size(), is_c2_reversed);
     
      auto iter = new_set.find(c1->get_contig_id());
      if(iter!=new_set.end())
	new_set.erase(iter);
      iter = new_set.find(c2->get_contig_id());
      if(iter!=new_set.end())
	new_set.erase(iter);

    } else {
      Contig * s = new_set.emplace(move(make_unique<Contig>(move(*c2),move(*c1),
							    is_c2_reversed, is_c1_reversed,
							    m->length(1),
							    new_set_id, data.contig_id_count++))).first->get();


      m_contig[id_c2]=make_tuple(s->get_contig_id(), 0, is_c2_reversed);
      m_contig[id_c1]=make_tuple(s->get_contig_id(), s->size()-c1->size(), is_c1_reversed);

      auto iter = new_set.find(c1->get_contig_id());
      if(iter!=new_set.end())
	new_set.erase(iter);
      iter = new_set.find(c2->get_contig_id());
      if(iter!=new_set.end())
	new_set.erase(iter);
      
    }
  }
}




void merge_match(Data &data, unsigned new_set_id,unsigned set_id1,
                 unsigned set_id2)
{
  auto &selected_matches = get<2>(data.matches[set_id1][set_id2]);
  vector<string> trash;
  set<unsigned long> absorbent;
  map<unsigned long, tuple<unsigned long, size_t, bool>> m_contig;
  auto v =  merge_full_matches(data,selected_matches, trash, m_contig,absorbent, new_set_id, set_id1, set_id2);
  merge_other_matches(data, v, trash, m_contig, new_set_id, set_id1, set_id2);

  auto &new_set =data.assembly_sets[new_set_id];

  // add contigs not selected in a match
  unsigned set_id[2] = {set_id1,set_id2};
  for(unsigned id : set_id){
    for(auto it=data.assembly_sets[id].begin();it!=data.assembly_sets[id].end();++it){
      if(m_contig.find(it->get()->get_contig_id())==m_contig.end()){

	Contig * s = new_set.emplace(move(make_unique<Contig>(move(*(it->get())),new_set_id,data.contig_id_count++))).first->get();
	  m_contig[it->get()->get_contig_id()]=make_tuple(s->get_contig_id(), 0, false);
      }
    }
  }

  
  
  cout << "Update matches" << endl;
  // update other matches
  for(auto &a : data.matches){
    if(a.first==set_id1){
      for(auto &b : a.second){
	if (b.first == set_id2)
          continue;
	for(auto &m : get<0>(b.second)){
	  update_match(data,m, m_contig, b.first, new_set_id);
	}
      }
    }
    else if(a.first<set_id1){
      for(auto &m : get<0>(a.second[set_id1])){
	update_match(data,m, m_contig, a.first, new_set_id);
      }
    }
    if(a.first==set_id2){
      for(auto &b : a.second){
        if (b.first == set_id1)
          continue;
	for(auto &m : get<0>(b.second)){
	  update_match(data,m, m_contig, b.first, new_set_id);}
      }
    }
    else if(a.first<set_id2){
      if(a.first==set_id1)
	continue;
      for(auto &m : get<0>(a.second[set_id1])){
	update_match(data, m, m_contig, a.first, new_set_id);
      }
    }
  }



  // for(auto it = new_set.begin(); it!=new_set.end();){
  //   if(m_contig.find(it->get())!=m_contig.end()){
  //     // remove intermediate contigs
  //     it=new_set.erase(it);
  //   } else ++it;
    
  // }

  data.assembly_sets.erase(set_id1);
  data.assembly_sets.erase(set_id2);
  

  data.matches.erase(set_id1);
  data.matches.erase(set_id2);

  for(auto &a : data.matches){
    a.second.erase(set_id1);
    a.second.erase(set_id2);
  }
}

void merge_algorithm(Data &data, bool greedy){

  for(size_t i=0; i<data.ids.size();++i)
    for(size_t j=i+1;j<data.ids.size();++j){
      auto &m = data.matches[i][j];
      if(greedy)
	get<1>(m) =  greedy_fill(get<2>(m),get<0>(m));
      else
	get<1>(m) =  algo_approx(get<2>(m),get<0>(m));
    }

  unsigned new_id=data.ids.size();
  unsigned cpt_id=data.ids.size();
  for(unsigned cpt=0;cpt<cpt_id-1;++cpt){
  
    unsigned max_score=0;
    unsigned i=data.matches.begin()->first,
      j=data.matches.begin()->second.begin()->first;
    
    // select two sets with a maximum score
    for(auto &a : data.matches)
      for(auto &b : a.second) {
	if(max_score<= get<1>(b.second)){
	  max_score=get<1>(b.second);
	  i=a.first;
	  j=b.first;
	}
      }
    string si, sj;
    for(auto &p : data.ids){
      if(p.second==i)
	si=p.first;
      if(p.second==j)
	sj=p.first;
    }
    cout << "Merge " << si << " with " << sj << endl;
    string new_set_name = si;
    new_set_name.append("|");
    new_set_name.append(sj);
    data.ids[new_set_name]=new_id;

    merge_match(data, new_id, i ,j);
    
    
    for(auto & a : data.assembly_sets){
      if(a.first==new_id) continue;
      auto &m = data.matches[a.first][new_id];
      if(greedy)
	get<1>(m) =  greedy_fill(get<2>(m), get<0>(m));
      else
	get<1>(m) = algo_approx(get<2>(m), get<0>(m));
    }
    new_id++;
  }
}

void merge_algorithm_general_greedy(Data &data)
{
  map<unsigned long, tuple<unsigned long, size_t, bool>> m_contig;

  while(!data.all_matches.empty()){
    const unique_ptr<Match> &p = *(data.all_matches.begin());
    Match * m = data.all_matches.begin()->get();
    if(m->is_full()){
      Contig * c = merge_full_match(data, m, m_contig);
      // remove matches containing the absorbed contig
      for(auto it=next(data.all_matches.begin()), last=data.all_matches.end(); it != last;){	
	if(it->get()->contains(c)) {
	  it=data.all_matches.erase(it);
	}
	else ++it;
      }
      data.contigs.erase(data.contigs.find(c->get_contig_id()));
    }
    else {
      pair<size_t,bool> shift_c1,shift_c2;
      Contig * c = merge_non_full_match(data, m,shift_c1,shift_c2);
      vector<unique_ptr<Match>> insert;
      for(auto it=next(data.all_matches.begin()), last = data.all_matches.end(); it != last;){
	Match * m2 = it->get();
	unsigned y=0;
	if(it->get()->contains(m->contig((unsigned)0))  
	   && it->get()->contains(m->contig(1))){
	  it=data.all_matches.erase(it);
	}
	else if(it->get()->contains(m->contig((unsigned)0))){
	  update_match(data,*it->get(),shift_c1,m->contig((unsigned)0),c,insert);
	  it=data.all_matches.erase(it);
	  
	}
	else if(it->get()->contains(m->contig(1))){
	  update_match(data,*it->get(),shift_c2,m->contig(1),c,insert);
	  it=data.all_matches.erase(it);
	}
	else ++it;
      }
      for(auto &u : insert){
	data.all_matches.insert(move(u));
      }
      data.contigs.erase(data.contigs.find(m->contig((unsigned)0)->get_contig_id()));
      data.contigs.erase(data.contigs.find(m->contig(1)->get_contig_id()));
    }
    data.all_matches.erase(data.all_matches.find(p));
  }
}
