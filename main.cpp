#include <algorithm>
#include <memory>
#include <regex>
#include <string>

#include "Contig.h"
#include "merge_contigs.cpp"
#include "Parser.hpp"
#include "options.hpp"


using namespace std;

string compact(const string &s1, const string &s2)
{
  for(size_t i1 = 0, i2=0; i1!=s1.size() &&  i2!=s2.size();){
    string str_tmp1, str_tmp2;
    for(;str_tmp1.back()!='_' && i1<s1.size();i1++)
      str_tmp1.push_back(s1[i1]);
    for(;str_tmp2.back()!='_' && i2<s2.size();i2++)
      str_tmp2.push_back(s2[i2]);
    if(str_tmp1!=str_tmp2){
      str_tmp1.pop_back();
      str_tmp2.pop_back();
      str_tmp1.append("!");
      str_tmp1.append(str_tmp2);
      return str_tmp1;
    }
  }
  return "null";
}

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

if(options.mode==Options::greedy_general)
  treatFastaDirectory(options.contig_directory.c_str(), data.ids, data.contigs,
		      options.pairwise, data.contig_id_count, m_contig,options.prioritize);
else treatFastaDirectory(options.contig_directory.c_str(), data.assembly_sets, data.ids,options.pairwise,data.contig_id_count,m_contig);
  if(options.generate_matches)
    createMatchDirectory(options.matches_directory.c_str(), options.contig_directory.c_str(), data.ids);
  
  std::cout << "Consensus construction" << std::endl;

  if(options.pairwise){
    for(auto &p : data.ids){
      for(auto &p2 : data.ids){
	unsigned long contig_id_count=0;
	if(p2.second<=p.second) continue;

	Data data_tmp(options.percentage_of_match);

	map<string,Contig*> m_contig;
	string tmp = options.contig_directory;
	tmp.append("/");
	tmp.append(p.first);
	tmp.append(".fasta");

	if(options.mode==Options::greedy_general){
	  parseFastaFile(tmp.c_str(), 0, data_tmp.contigs, data_tmp.contig_id_count,m_contig);
	}
	else parseFastaFile(tmp.c_str(), 0, data_tmp.assembly_sets ,contig_id_count,m_contig);

	tmp = options.contig_directory;
	tmp.append("/");
	tmp.append(p2.first);
	tmp.append(".fasta");

	if(options.mode==Options::greedy_general)
	  parseFastaFile(tmp.c_str(), 1, data_tmp.contigs, data_tmp.contig_id_count,m_contig);
	else parseFastaFile(tmp.c_str(), 1, data_tmp.assembly_sets,contig_id_count,m_contig);

	data_tmp.ids[p.first]=0;
	data_tmp.ids[p2.first]=1;

	MatchMatrix m;
	if(options.mode==Options::greedy_general){
	  treatMatchDirectory(options.matches_directory.c_str(), data_tmp.ids, data_tmp.all_matches, options.percentage_of_match,m_contig,options.prioritize, false);
	}
	else treatMatchDirectory(options.matches_directory.c_str(), data_tmp.assembly_sets, data_tmp.ids, m,options.percentage_of_match,m_contig,false);
	m_contig.clear();
	
	if(options.mode==Options::greedy_general)
	  merge_algorithm_general_greedy(data_tmp);
	else merge_algorithm(data_tmp,options.mode==Options::greedy_pairwise);
	tmp = options.output_folder;
	tmp.append("/");


	string comp = compact(p.first, p2.first);

	tmp.append(comp);
	if(options.mode==Options::greedy_general)
	  tmp.append("_gg_");
	else tmp.append("_pg_");
	tmp.append(to_string((int)options.percentage_of_match));
	
	// tmp.append(p.first);
	// tmp.append("!");
	// tmp.append(p2.first);

	string output = tmp;
	output.append(".fasta");

	if(options.mode==Options::greedy_general)
	  output_contig(output.c_str(),data_tmp.contigs,options.random_when_tie);
	else output_contig(output.c_str(), data_tmp.assembly_sets, options.random_when_tie);

	output = tmp;
	output.append(".cons");
	//	output_contig_ordering(output.c_str(), data_tmp.assembly_sets,data_tmp.ids);
	
      }
    }
  }else {
    if(options.mode==Options::greedy_general)
      treatMatchDirectory(options.matches_directory.c_str(), data.ids, data.all_matches, options.percentage_of_match,m_contig, options.prioritize);
    else treatMatchDirectory(options.matches_directory.c_str(), data.assembly_sets, data.ids, data.matches, options.percentage_of_match,m_contig);
    m_contig.clear();
    if(options.mode==Options::greedy_general)
      merge_algorithm_general_greedy(data);
    else merge_algorithm(data,options.mode==Options::greedy_pairwise);

    if(options.mode==Options::greedy_general)
      output_contig(options.output_fast.c_str(),data.contigs,options.random_when_tie,options.prioritize!="-1");
    else output_contig(options.output_fast.c_str(),data.assembly_sets,options.random_when_tie);
    //output_contig_ordering(options.output_cons.c_str(),data.assembly_sets, data.ids);
  }
  return 0;
}

