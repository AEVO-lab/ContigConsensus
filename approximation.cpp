#pragma once

#include <iostream>
#include <map>

#include "Contig.h"
#include "Hungarian_algo.hpp"
#include "Match.h"
#include "Parser.hpp"
#include "MatchesPerContig.hpp"
#include "Nucleotide.h"

using namespace std;

unsigned greedy_fill(vector<const Match*>& selected_matches, vector<Match>& all_matches){

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


queue<const Match*> construct_BEO(Match::SortedLengthSet& matches)
{

  queue<const Match*> beoMatches;
  unsigned nb_non_full_contigs=0;
  map<const Contig*,MatchesPerContig> mpc;
  set<const Match*> non_full_matches;

  for (const Match *m : matches){
    mpc.emplace(m->contig((unsigned)0), MatchesPerContig(0));
    mpc.emplace(m->contig(1), MatchesPerContig(1));

    nb_non_full_contigs+=mpc.at(m->contig((unsigned)0)).add(m);
    nb_non_full_contigs+=mpc.at(m->contig(1)).add(m);
  }
  while (!matches.empty()){
    if(nb_non_full_contigs==0){
      for (const Match *m : matches) 
        beoMatches.push(m);
      
      return beoMatches;
    }

		
    auto p = *find_if(mpc.begin(),mpc.end(),[](const pair<const Contig *, MatchesPerContig> &p){
      return !p.second.only_full();
    });

            
    const Contig * cI = p.first;
    const Match* IJ = p.second.early();
    const Contig *cJ; // equivalent to w

    while(true){

      cJ = IJ->contig(cI);
	
      if (mpc.at(cJ).is_early_or_late(IJ)) 
	break;
      

      if(IJ->projected_end(cJ->get_set_id())==cJ->size()-1){
	cI=cJ; 
	IJ= mpc.at(cI).late();
      }
      else {
	cI=cJ;
	IJ= mpc.at(cI).early();
      }
    }
    beoMatches.push(IJ);

    nb_non_full_contigs-= mpc.at(cI).remove(IJ);
    nb_non_full_contigs-= mpc.at(cJ).remove(IJ);
    matches.erase(IJ);
  }
  return beoMatches;
}

vector<const Match*> local_ratio_rec(queue<const Match*> &beo, vector<pair<const Match*,int>> &sigma, unsigned& score)
{
  if(beo.empty())
    return vector<const Match*>();
  
  const  Match * F = beo.front();
  beo.pop();

  unsigned score_f = F->score;
  for(auto &p : sigma){
    if(F->intersect(p.first)) // what if p.second<0?
      score_f-= p.second;
  }

  if(score_f>0)
    sigma.push_back(make_pair(F, score_f));
  
  vector<const Match*> M = local_ratio_rec(beo,sigma,score);
  
  if (score_f>0 && !F->intersect(M)) {
    M.push_back(F);
    score+=F->score;
  }

  return M;
}


vector<const Match*> local_ratio(queue<const Match*> &beo, unsigned& score)
{
  vector<pair<const Match*,int>> sigma;
  return local_ratio_rec(beo,sigma,score);
}


vector<const Match*> greedy_maximum_matching(const Match::MM_map& matches, unsigned &score)
{
  score=0;
  vector<const Match*> selected_matches;
  for(auto p1 : matches)
    for(auto p2 : p1.second)
      selected_matches.push_back(p2.second);

  sort(selected_matches.begin(),selected_matches.end(),[](const Match *m1, const Match* m2){
    return m1->score>m2->score;
  });

  vector<const Match*> rtrn;
  for(const Match * m : selected_matches){
    if (!m->intersect(rtrn)) {
      rtrn.push_back(m);
      score += m->score;
    }
  }
  return rtrn;
}

unsigned algo_approx(vector<const Match*>& selected_matches, vector<Match>& all_matches)
{
  Match::MM_map max_matches;
  unsigned maximum_value=0;


  vector<const Match*> selected_matches_best;
  unsigned score_max=0;
  unsigned score_tmp;
  vector<const Match*> selected_matches_tmp;
  
  for(unsigned i=0;i<2;++i){
  
    auto m = Match::match_filter(i, all_matches);
    auto beo = construct_BEO(m);
  

    score_tmp=0;
    selected_matches_tmp = local_ratio(beo,score_tmp);

    score_tmp+= greedy_fill(selected_matches_tmp, all_matches);
    if(score_tmp>score_max){
      score_max=score_tmp;
      selected_matches_best=move(selected_matches_tmp);
    }
  }
  
  
  for(Match &m : all_matches){
    try{
      if(max_matches.at(m.contig((unsigned)0)).at(m.contig(1))->score < m.score) {
	max_matches[m.contig((unsigned)0)][m.contig(1)] = &m;
	maximum_value = max(maximum_value, m.score);
      }
    }
    catch(out_of_range) {
      max_matches[m.contig((unsigned)0)][m.contig(1)] = &m;
      maximum_value = max(maximum_value, m.score);
    }
  }
  
  score_tmp=0;
  
  selected_matches_tmp = greedy_maximum_matching(max_matches, score_tmp);  //hungarian_algorithm(S, T, max_matches, maximum_value,score_m3);
  score_tmp+= greedy_fill(selected_matches_tmp, all_matches);
  

  if(score_tmp>score_max){
    score_max=score_tmp;
    selected_matches_best=move(selected_matches_tmp);
  }

  selected_matches=move(selected_matches_best);
  return score_max;
}
