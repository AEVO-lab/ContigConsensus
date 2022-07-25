#pragma once
#include <algorithm>
#include <dirent.h>
#include <memory>


#include "Contig.h"
#include "Match.h"


using namespace std;

#define DEBUG false
#define ABSORBED 0

#define SEPARATOR_CHAR '!'
#define SEPARATOR_STR "!"


class ConsensusConstruction
{
  map<string,Contig*> m_contig;

protected:
  unsigned long contig_id_count;
  map<string,unsigned> ids;
  float percentage;
  
protected:
  ConsensusConstruction(float percentage) : percentage(percentage), contig_id_count(1) {}
  
  virtual set<unique_ptr<Contig>,less<>> & get_consensus_result()=0;
  
  virtual Contig * construct_contig(unsigned set_id,unsigned long contig_id_count, string && name,vector<Nuc> && seq, bool prioritize)=0;

  virtual void construct_match(Contig *c1, Contig *c2,size_t start1, size_t end1, size_t start2, size_t end2, unsigned score)=0;


private:
  
  void parseFastaFile(const char* fileName, unsigned set_id, bool prioritize);

  void parseMatches(const char * fileName, unsigned id_set1, unsigned id_set2);


public:

  

  virtual ~ConsensusConstruction() {}

  virtual void merge_algorithm()=0;
  
  void treatFastaDirectory(const char *dirName, bool only_read_names, string prioritize);

  void treatMatchDirectory(const char *dirName, bool display=true);

  void output_contigs(const char* fileName, bool random_when_tie, bool prioritize);

  void output_contig_ordering(const char * fileName, bool prioritize);

  void createMatchDirectory(const char *m_dirName, const char *f_dirName);


};   
