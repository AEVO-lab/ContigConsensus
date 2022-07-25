#include "ConsensusConstruction.cpp"

class Greedy : public ConsensusConstruction
{
  set<unique_ptr<Match>,less<>> all_matches;
  set<unique_ptr<Contig>,less<>> contigs;

  virtual set<unique_ptr<Contig>,less<>> & get_consensus_result();

  virtual Contig * construct_contig(unsigned set_id,unsigned long contig_id_count, string && name,vector<Nuc> && seq, bool prioritize);

  virtual void construct_match(Contig *c1, Contig *c2,size_t start1, size_t end1, size_t start2, size_t end2, unsigned score);


  Contig* merge_full_match(Match *m,map<unsigned long,tuple<unsigned long,size_t,bool>> &m_contig);

  Contig* merge_non_full_match(Match *m, pair<size_t,bool> &shift_c1, pair<size_t,bool> &shift_c2);

  void update_match(Match &m, pair<size_t,bool> shift, Contig * old_contig, Contig * new_contig, vector<unique_ptr<Match>> &insert);
  
public:

  Greedy(float percentage);

  virtual ~Greedy() {}
  void merge_algorithm();


};                  
