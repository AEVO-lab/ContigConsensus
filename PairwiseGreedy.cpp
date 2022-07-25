#include "PairwiseGreedy.hpp"

set<unique_ptr<Contig>,less<>> & PairwiseGreedy::get_consensus_result()
{
  return  assembly_sets.begin()->second;
}
  
Contig * PairwiseGreedy::construct_contig(unsigned set_id,unsigned long contig_id_count, string && name,vector<Nuc> && seq, bool prioritize)
{
  return assembly_sets[set_id].emplace(make_unique<Contig>(set_id,contig_id_count,move(name),move(seq),prioritize)).first->get();
}

void PairwiseGreedy::construct_match(Contig *c1, Contig *c2,size_t start1, size_t end1, size_t start2, size_t end2, unsigned score)
{
  unsigned set_id1 = min(c1->get_set_id(),c2->get_set_id());
  unsigned set_id2 = max(c1->get_set_id(),c2->get_set_id());
  get<0>(matches[set_id1][set_id2]).emplace_back(c1,c2,start1,end1,start2,end2,score);
}


unsigned PairwiseGreedy::greedy_fill(vector<const Match*>& selected_matches, vector<Match>& all_matches){

  unsigned added_score=0;
  sort(all_matches.begin(),all_matches.end(),[](const Match &m1, const Match& m2){
    return m1.score>m2.score;
  });

  for(Match & m : all_matches){
    if (!m.intersect(selected_matches)) {
      selected_matches.push_back(&m);
      added_score += m.score;
    }
  }
  return added_score;
}


vector<const Match*> PairwiseGreedy::merge_full_matches(vector<const Match *> &selected_matches,
							map<unsigned long,tuple<unsigned long,size_t,bool>> &m_contig, set<unsigned long> &absorbent,
							unsigned new_set_id, unsigned set_id1, unsigned set_id2)
{
  auto &new_set =assembly_sets[new_set_id];

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
      overlap = m->contig(1)->size()-m->projected_start(1);
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

void PairwiseGreedy::merge_other_matches(vector<const Match *> &selected_matches,
					 map<unsigned long,tuple<unsigned long,size_t, bool>>& m_contig,
					 unsigned new_set_id,unsigned set_id1,
					 unsigned set_id2)
{

  unsigned it=1;
  auto &new_set =assembly_sets[new_set_id];

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
							    m->length(1),new_set_id, contig_id_count++))).first->get();
     
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
							    new_set_id, contig_id_count++))).first->get();


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




void PairwiseGreedy::merge_match(unsigned new_set_id,unsigned set_id1,
				 unsigned set_id2)
{
  auto &selected_matches = get<2>(matches[set_id1][set_id2]);
  set<unsigned long> absorbent;
  map<unsigned long, tuple<unsigned long, size_t, bool>> m_contig;
  auto v =  merge_full_matches(selected_matches, m_contig,absorbent, new_set_id, set_id1, set_id2);
  merge_other_matches(v, m_contig, new_set_id, set_id1, set_id2);

  auto &new_set =assembly_sets[new_set_id];

  // add contigs not selected in a match
  unsigned set_id[2] = {set_id1,set_id2};
  for(unsigned id : set_id){
    for(auto it=assembly_sets[id].begin();it!=assembly_sets[id].end();++it){
      if(m_contig.find(it->get()->get_contig_id())==m_contig.end()){

	Contig * s = new_set.emplace(move(make_unique<Contig>(move(*(it->get())),new_set_id,contig_id_count++))).first->get();
	m_contig[it->get()->get_contig_id()]=make_tuple(s->get_contig_id(), 0, false);
      }
    }
  }

  cout << "Update matches" << endl;
  // update other matches
  for(auto &a : matches){
    if(a.first==set_id1){
      for(auto &b : a.second){
	if (b.first == set_id2)
          continue;
	for(auto &m : get<0>(b.second)){
	  update_match(m, m_contig, b.first, new_set_id);
	}
      }
    }
    else if(a.first<set_id1){
      for(auto &m : get<0>(a.second[set_id1])){
	update_match(m, m_contig, a.first, new_set_id);
      }
    }
    if(a.first==set_id2){
      for(auto &b : a.second){
        if (b.first == set_id1)
          continue;
	for(auto &m : get<0>(b.second)){
	  update_match(m, m_contig, b.first, new_set_id);}
      }
    }
    else if(a.first<set_id2){
      if(a.first==set_id1)
	continue;
      for(auto &m : get<0>(a.second[set_id1])){
	update_match(m, m_contig, a.first, new_set_id);
      }
    }
  }

  assembly_sets.erase(set_id1);
  assembly_sets.erase(set_id2);
  

  matches.erase(set_id1);
  matches.erase(set_id2);

  for(auto &a : matches){
    a.second.erase(set_id1);
    a.second.erase(set_id2);
  }
}


void PairwiseGreedy::update_match(Match &m, map<unsigned long,tuple<unsigned long,size_t,bool>>& m_contig, unsigned other_set_id, unsigned new_set_id)
{
  auto & new_set=assembly_sets[new_set_id];
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

  get<0>(matches[other_set_id][new_set_id]).emplace_back(Match(c1,c2,start_c1,is_c1_reversed!=m.is_reverse((unsigned)0),
							       start_c2,is_c2_reversed!=m.is_reverse(1)));

  auto & back =get<0>(matches[other_set_id][new_set_id]).back();

  
  if(back.score*100.0/(float)back.length((unsigned)0) < percentage)
    get<0>(matches[other_set_id][new_set_id]).pop_back();
}


PairwiseGreedy::PairwiseGreedy(float percentage) : ConsensusConstruction(percentage) {}

void PairwiseGreedy::merge_algorithm(){

  for(size_t i=0; i<ids.size();++i)
    for(size_t j=i+1;j<ids.size();++j){
      auto &m = matches[i][j];
      get<1>(m) =  greedy_fill(get<2>(m),get<0>(m));
    }

  unsigned new_id=ids.size();
  unsigned cpt_id=ids.size();
  for(unsigned cpt=0;cpt<cpt_id-1;++cpt){
  
    unsigned max_score=0;
    unsigned i=matches.begin()->first,
      j=matches.begin()->second.begin()->first;
    
    // select two sets with a maximum score
    for(auto &a : matches)
      for(auto &b : a.second) {
	if(max_score<= get<1>(b.second)){
	  max_score=get<1>(b.second);
	  i=a.first;
	  j=b.first;
	}
      }
    string si, sj;
    for(auto &p : ids){
      if(p.second==i)
	si=p.first;
      if(p.second==j)
	sj=p.first;
    }
    cout << "Merge " << si << " with " << sj << endl;
    string new_set_name = si;
    new_set_name.append("|");
    new_set_name.append(sj);
    ids[new_set_name]=new_id;

    merge_match(new_id, i ,j);
    
    
    for(auto & a : assembly_sets){
      if(a.first==new_id) continue;
      auto &m = matches[a.first][new_id];
      get<1>(m) =  greedy_fill(get<2>(m), get<0>(m));

    }
    new_id++;
  }
}

  

