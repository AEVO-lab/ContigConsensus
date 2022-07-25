#include "ConsensusConstruction.hpp"
#include "Contig.h"


class PairwiseGreedy : public ConsensusConstruction
{
  AssemblySet assembly_sets;
  MatchMatrix matches;


  virtual set<unique_ptr<Contig>,less<>> & get_consensus_result();
  
  virtual Contig * construct_contig(unsigned set_id,unsigned long contig_id_count, string && name,vector<Nuc> && seq, bool prioritize);

  virtual void construct_match(Contig *c1, Contig *c2,size_t start1, size_t end1, size_t start2, size_t end2, unsigned score);


  unsigned greedy_fill(vector<const Match*>& selected_matches, vector<Match>& all_matches);


  vector<const Match*> merge_full_matches(vector<const Match *> &selected_matches,
					map<unsigned long,tuple<unsigned long,size_t,bool>> &m_contig, set<unsigned long> &absorbent,
					  unsigned new_set_id, unsigned set_id1, unsigned set_id2);


  void merge_other_matches(vector<const Match *> &selected_matches,
			 map<unsigned long,tuple<unsigned long,size_t, bool>>& m_contig,
			 unsigned new_set_id,unsigned set_id1,
			   unsigned set_id2);


  void merge_match(unsigned new_set_id,unsigned set_id1,
		   unsigned set_id2);


  void update_match(Match &m, map<unsigned long,tuple<unsigned long,size_t,bool>>& m_contig, unsigned other_set_id, unsigned new_set_id);

public:

  PairwiseGreedy(float percentage);

  void merge_algorithm();

  
} ;   
