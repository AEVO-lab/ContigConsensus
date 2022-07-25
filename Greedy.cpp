#include "Greedy.hpp"


set<unique_ptr<Contig>,less<>> & Greedy::get_consensus_result()
{
  return contigs;
}

Contig * Greedy::construct_contig(unsigned set_id,unsigned long contig_id_count, string && name,vector<Nuc> && seq, bool prioritize){
  return contigs.emplace(make_unique<Contig>(set_id,contig_id_count,move(name),move(seq),prioritize)).first->get();
}

void Greedy::construct_match(Contig *c1, Contig *c2,size_t start1, size_t end1, size_t start2, size_t end2, unsigned score)
{
  all_matches.emplace(move(make_unique<Match>(c1,c2,start1,end1,start2,end2,score)));
}



Contig* Greedy::merge_full_match(Match *m,map<unsigned long,tuple<unsigned long,size_t,bool>> &m_contig)
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
    c1->absorb(*c2, shift, reverse);
    return c2;
  }


  Contig* Greedy::merge_non_full_match(Match *m, pair<size_t,bool> &shift_c1, pair<size_t,bool> &shift_c2)
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
      Contig * s = contigs.emplace(move(make_unique<Contig>(move(*c1),move(*c2),
							    is_c1_reversed,
							    is_c2_reversed,
							    m->length(1),0, contig_id_count++))).first->get();
      shift_c1 = make_pair(0, is_c1_reversed);
      shift_c2 = make_pair(s->size()-c2->size(), is_c2_reversed);
      return s;
		      
    } 

    Contig * s = contigs.emplace(move(make_unique<Contig>(move(*c2),move(*c1),
							  is_c2_reversed, is_c1_reversed,
							  m->length(1),
							  0, contig_id_count++))).first->get();

    shift_c1=make_pair(s->size()-c1->size(), is_c1_reversed);
    shift_c2=make_pair(0, is_c2_reversed);
    return s;

  }

  void Greedy::update_match(Match &m, pair<size_t,bool> shift, Contig * old_contig, Contig * new_contig, vector<unique_ptr<Match>> &insert)
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
    if(tmp->score*100/(float)tmp->length(0) >=percentage)
      insert.push_back(move(tmp));

  }
    
  Greedy::Greedy(float percentage) : ConsensusConstruction(percentage) {}

  void Greedy::merge_algorithm()
{
  map<unsigned long, tuple<unsigned long, size_t, bool>> m_contig;


  while(!all_matches.empty()){
    const unique_ptr<Match> &p = *(all_matches.begin());
    Match * m = all_matches.begin()->get();


    if(m->is_full()){
      Contig * c = merge_full_match(m, m_contig);
      // remove matches containing the absorbed contig
      for(auto it=next(all_matches.begin()), last=all_matches.end(); it != last;){	
	if(it->get()->contains(c)) {
	  it=all_matches.erase(it);
	}
	else ++it;
      }
      contigs.erase(contigs.find(c->get_contig_id()));
    }
    else {
      pair<size_t,bool> shift_c1,shift_c2;
      Contig * c = merge_non_full_match(m,shift_c1,shift_c2);
      vector<unique_ptr<Match>> insert;
      for(auto it=next(all_matches.begin()), last = all_matches.end(); it != last;){
	Match * m2 = it->get();
	unsigned y=0;
	if(it->get()->contains(m->contig((unsigned)0))  
	   && it->get()->contains(m->contig(1))){
	  it=all_matches.erase(it);
	}
	else if(it->get()->contains(m->contig((unsigned)0))){
	  update_match(*it->get(),shift_c1,m->contig((unsigned)0),c,insert);
	  it=all_matches.erase(it);
	  
	}
	else if(it->get()->contains(m->contig(1))){
	  update_match(*it->get(),shift_c2,m->contig(1),c,insert);
	  it=all_matches.erase(it);
	}
	else ++it;
      }
      for(auto &u : insert){
	all_matches.insert(move(u));
      }
      contigs.erase(contigs.find(m->contig((unsigned)0)->get_contig_id()));
      contigs.erase(contigs.find(m->contig(1)->get_contig_id()));
    }
    all_matches.erase(all_matches.find(p));
  }
}
