#include "Contig.h"
#include "Match.h"
#include "approximation.cpp"
#include "Parser.hpp"
#include "options.hpp"
#include <memory>


using namespace std;

AssemblySet copy_pair(AssemblySet & assembly_set, unsigned id1, unsigned id2)
{
  AssemblySet tmp;
  for(auto &c : assembly_set[id1])
    tmp[id1].insert(make_unique<Contig>(*c.get()));
  for(auto &c : assembly_set[id2])
    tmp[id2].insert(make_unique<Contig>(*c.get()));

  return tmp;
}


int main(int argc, char *argv[])
{

  Options options(argc,argv);
  Data data(options.percentage_of_match);

  
  map<string,Contig*> m_contig;
  treatFastaDirectory(options.contig_directory.c_str(), data.assembly_sets, data.ids,options.every_pair,data.contig_id_count,m_contig);
  if(options.generate_matches)
    createMatchDirectory(options.matches_directory.c_str(), options.contig_directory.c_str(), data.ids);
  

  std::cout << "Consensus construction" << std::endl;

  if(options.every_pair){
    for(auto &p : data.ids){
      for(auto &p2 : data.ids){
	unsigned long contig_id_count=0;
	if(p2.second<=p.second) continue;
	AssemblySet a_s;
	map<string,Contig*> m_contig;
	string tmp = options.contig_directory;
	tmp.append("/");
	tmp.append(p.first);
	tmp.append(".fasta");

	
	parseFastaFile(tmp.c_str(), 0, a_s,contig_id_count,m_contig);

	tmp = options.contig_directory;
	tmp.append("/");
	tmp.append(p2.first);
	tmp.append(".fasta");
	parseFastaFile(tmp.c_str(), 1, a_s,contig_id_count,m_contig);

	map<string,unsigned> ids_tmp;
	ids_tmp[p.first]=0;
	ids_tmp[p2.first]=1;

	MatchMatrix m;
	treatMatchDirectory(options.matches_directory.c_str(), a_s, ids_tmp, m,options.percentage_of_match,m_contig,false);

	
	merge_algorithm(data);
	tmp = options.output_folder;
	tmp.append("/");
	tmp.append(p.first);
	tmp.append("!");
	tmp.append(p2.first);

	string output = tmp;
	output.append(".fasta");
	output_contig(output.c_str(), a_s, options.random_when_tie);

	output = tmp;
	output.append(".cons");
	output_contig_ordering(output.c_str(), a_s,ids_tmp);
	
      }
    }
  }else {
    treatMatchDirectory(options.matches_directory.c_str(), data.assembly_sets, data.ids, data.matches, options.percentage_of_match,m_contig);
    m_contig.clear();
    merge_algorithm(data);

    output_contig(options.output_fast.c_str(),data.assembly_sets,options.random_when_tie);
    output_contig_ordering(options.output_cons.c_str(),data.assembly_sets, data.ids);
  }
  return 0;
}

