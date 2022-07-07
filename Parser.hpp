#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <unistd.h>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Contig.h"
#include "Match.h"

#define DEBUG false

#define SEPARATOR_CHAR '!'
#define SEPARATOR_STR "!"

using namespace std;


void parseFastaFile(const char* fileName, unsigned set_id, AssemblySet &assembly_sets, unsigned long &contig_id_count,map<string,Contig*> &m_contig)
{
  size_t max_size=0;
  string max_name;
  ifstream file(fileName);
  if(!file)
    {
      cout << "Error: could not open " << fileName << endl;
      exit(EXIT_FAILURE);
    }
  
  auto &v = assembly_sets[set_id];
  file.get(); // do not read the first '>'
  while(!file.eof())
    {
      vector<Nuc> seq;
       
      // ignore commentary line
      //		file.ignore(numeric_limits<streamsize>::max(),'\n');
      string line;
      getline(file,line);

      string name = line.substr(0,line.find_first_of(" \n"));

      int c;
      do {
	c = file.get();
	if(c == '\n' || c == EOF|| c == '>' || c == '*') continue;
	if((char)c != 'N' && (char)c != 'A' && (char)c != 'T' && (char)c != (char)'C' && c != (char)'G') continue;
	seq.push_back(c);
			
      } while (c!='>' && c != EOF);
      if(seq.size()>max_size){
	max_size=seq.size();
	max_name=name;
      }

      if(seq.size()>0)
	m_contig[name]=v.emplace(make_unique<Contig>(set_id,contig_id_count++,move(name),move(seq))).first->get();


		
    }

#if DEBUG==true
  for(const auto &cont : v)
    cout << *cont.get() << "\n";
#endif
  cout << fileName << ": " << max_name << endl;
  cout << "size: " << max_size << endl;
  

}


void parseMatches(const char * fileName, AssemblySet & assembly_sets, MatchMatrix & matches,
		  unsigned id_set1, unsigned id_set2,float percentage,map<string,Contig*> &m_contig)
{
  ifstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }

  // insert keys
  auto &v = get<0>(matches[id_set1][id_set2]);
  
  size_t score, startS, endS, lengthS, startT, endT, lengthT;
  string nameS, nameT;
  while (file >> score) {

    file >> nameS;
    file >> startS;
    file >> endS;
    file >> lengthS;
    file >> nameT;
    file >> startT;
    file >> endT;
    file >> lengthT;


    auto tmp = m_contig.find(nameS);
    if(tmp==m_contig.end()){
      	cout << "error " << nameS << " doesn't exist in "<< id_set1 << endl;
	exit(EXIT_FAILURE);
    }
    Contig *s = tmp->second;
    assert(s->get_set_id()==id_set1 || s->get_set_id()==id_set2);

    tmp = m_contig.find(nameT);
    if(tmp==m_contig.end()){
      	cout << "error " << nameT << " doesn't exist in "<< id_set2 << endl;
	exit(EXIT_FAILURE);
    }
    Contig *t = tmp->second;
    assert(t->get_set_id()==id_set2 || t->get_set_id()==id_set1);
    // TODO use pair<string,unsigned> for the map key

    while(startS != 1 && startS != lengthS && startT !=1 && startT !=lengthT){
      startS = startS<endS ? startS-1 : startS+1;
      startT = startT<endT ? startT-1 : startT+1;      
      if(s->getNuc(startS-1).is_equal(t->getNuc(startT-1),startT>endT))
	score++;
    }
    while(endS != 1 && endS != lengthS && endT !=1 && endT !=lengthT){
      endS = startS<endS ? endS+1 : endS-1;
      endT = startT<endT ? endT+1 : endT-1;
      if(s->getNuc(endS-1).is_equal(t->getNuc(endT-1),startT>endT))
	score++;
    }

    unsigned lengthS = startS<endS ? endS-startS : startS-endS;
    unsigned lengthT = startT<endT ? endT-startT : startT-endT;
    assert(lengthS==lengthT);
  
    if((float)score*100.0/(lengthS+1)>=percentage){
      if(s->get_set_id()<t->get_set_id())
	v.emplace_back(s,t,startS-1,endS-1,startT-1,endT-1,score);
      else
	v.emplace_back(t,s,startT-1,endT-1,startS-1,endS-1,score);

      //v.back().display_alignment();

      //  if(t->getName()=="NODE_31159_length_688_cov_50.574127" && s->getName()=="333"){
      // 	v.back().display_alignment();
      // 	exit(EXIT_SUCCESS);
      //  }
      // if(t->getName()=="NODE_31159_length_688_cov_50.574127" && s->getName()=="384"){
      // 	v.back().display_alignment();

      // }

    }
    //if(it==15) exit(EXIT_SUCCESS);
    // if(lengthS <=40){
    //   std::cout << lengthS << endl;
    //   Match & m = v.back();
    //   m.display_contig_names();
    //   m.display_alignment();
    //   exit(EXIT_SUCCESS);
    // }

    
    
  }
  //  std::cout << v.size() << "/" << it << std::endl;
}


void treatFastaDirectory(const char *dirName, AssemblySet &assembly_sets, map<string,unsigned> &ids,
			 bool only_read_names, unsigned long& contig_id_count, map<string,Contig*> &m_contig)
{
  unsigned id=0;
  DIR* dir;
  struct dirent *ent;
  dir=opendir(dirName);
  if (dir == NULL) {
    cout << "Error: could not open directory " << dirName << "\n";
    exit(EXIT_FAILURE);
  }

  while((ent=readdir(dir))!=NULL){
    if(ent->d_type!=DT_REG)
      continue;
		
    string f_name = ent->d_name;
    if(f_name.find_last_of(".")!= string::npos && f_name.substr(f_name.find_last_of(".")) != ".fasta")
      continue;
		
    cout << "Read fasta file: " << f_name <<"\n";
    string file = dirName;
    file.append("/");
    file.append(ent->d_name);
    ids[f_name.substr(0,f_name.find_last_of("."))]=id;
    if(!only_read_names) parseFastaFile(file.c_str(), id, assembly_sets, contig_id_count,m_contig);
    id++;
  }
}


void treatMatchDirectory(const char *dirName, AssemblySet &assembly_sets, map<string,unsigned> &ids,
			 MatchMatrix &matches, float percentage,map<string,Contig*>& m_contig, bool display=true)
{
  DIR* dir;
  struct dirent *ent;
  dir=opendir(dirName);
  if (dir == NULL) {
    cout << "Error: could not open directory " << dirName << "\n";
    exit(EXIT_FAILURE);
  }
  while((ent=readdir(dir))!=NULL){
    if(ent->d_type!=DT_REG)
      continue;
    
    string f_name = ent->d_name;
    if(f_name.substr(f_name.find_last_of(".")) != ".txt")
      continue;
		
   
    string s1 = f_name.substr(0,f_name.find_first_of(SEPARATOR_CHAR));
    string s2 = f_name.substr(f_name.find_first_of(SEPARATOR_CHAR)+1);
    s2 = s2.substr(0,s2.find_first_of('.'));

    auto id1 = ids.find(s1);
    if(id1==ids.end()){
      if(display) cout << "Error, there is no " << s1 <<".fasta file.\n Ignoring " << f_name << endl;
      continue;
    }
    auto id2 = ids.find(s2);
    if(id2==ids.end()){
      if(display) cout << "Error, there is no " << s2 <<".fasta file.\n Ignoring " << f_name << endl;
      continue;
    }
    
    
    string file = dirName;
    file.append("/");
    file.append(ent->d_name);

    cout << "Read match file: " << f_name <<"\n";
    if(id1->second<id2->second)
      parseMatches(file.c_str(), assembly_sets, matches,id1->second,id2->second,percentage,m_contig);
    else parseMatches(file.c_str(), assembly_sets, matches,id2->second,id1->second,percentage,m_contig);
  }
  
}

void createMatchDirectory(const char *m_dirName,const char *f_dirName,map<string,unsigned> &ids){
  //std::filesystem::create_directory("test");

  string cmd_mkdir = "mkdir ";
  cmd_mkdir.append(m_dirName);
  system(cmd_mkdir.c_str());

  cout << "Generating match files" << endl;
  for(auto it_a=ids.begin();it_a!=ids.end();++it_a){
    for(auto it_b=next(it_a);it_b!=ids.end();++it_b){
      string cmd_blast="blastn -task megablast -query "; //-strand plus
      cmd_blast.append(f_dirName);
      cmd_blast.append("/");
      cmd_blast.append(it_a->first);
      cmd_blast.append(".fasta -subject ");
      cmd_blast.append(f_dirName);
      cmd_blast.append("/");
      cmd_blast.append(it_b->first);
      cmd_blast.append(".fasta -ungapped -out \"");
      cmd_blast.append(m_dirName);
      cmd_blast.append("/");
      cmd_blast.append(it_a->first);
      cmd_blast.append(SEPARATOR_STR);
      cmd_blast.append(it_b->first);
      cmd_blast.append(".txt\" -outfmt \"6 score qseqid qstart qend qlen sseqid sstart send slen\"");
      system(cmd_blast.c_str());
    }
  }
}

void output(const char * fileName, vector<const Match*>& matches)
{

  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  for(size_t i = 0; i<matches.size(); ++i){
    file << "#Match" << i+1<<"-score-" << matches[i]->score<<endl;
    for(unsigned t = 0; t<=1; ++t){
      file <<">" << matches[i]->contig(t)->getName() << ": " 
	   << matches[i]->start(t) << " " << matches[i]->end(t) 
	   <<endl;

    }
    file << endl;
  }
}

void output_contig(const char* fileName, AssemblySet & assembly_sets,bool random_when_tie)
{
  unsigned max_size=0;
  string max_name;

  unsigned id=1;
  
  auto & contigs = assembly_sets.begin()->second;

  vector<Contig*> v_contigs;
  for(auto &c : contigs)
    v_contigs.push_back(c.get());

  sort(v_contigs.begin(),v_contigs.end(),[](Contig *c1,Contig *c2){return c1->size()<c2->size();});



  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  // write_contig(file, *v_contigs.back(), 1, !random_when_tie);
  // return;
  
  cout << "Writing fasta result in " << fileName << endl;
  for(auto &c : v_contigs){
    if(c->size()>max_size){
      max_size=c->size();
      max_name=c->getName();
    }
    if(c->getName().find("k119_79")!=string::npos) cout << "size k119: " << c->size() <<endl;
    //file << *c;
    write_contig(file, *c, id, !random_when_tie);
    id++;
  }
  cout << max_size << ": " << endl; //<< max_name << endl;
}

void output_contig_ordering(const char * fileName, AssemblySet & assembly_sets, map<string,unsigned> &ids)
{

  map<unsigned,string> r_map;
  for(auto &a : ids)
    r_map[a.second]=a.first;
  
  auto & contigs = assembly_sets.begin()->second;
  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Writing ordering result in " << fileName << endl;
  for(auto &c : contigs){
    file << c->getName() << ":"<<endl;
    for(auto & C : c->getComponent()){
      file << "name: " << C.name << " ";
      file << "file: " << r_map[C.set_id] << " ";
      if(C.is_reversed){
	file << "position: " << C.shift+1-C.contig_size << " ";
	file << "reversed";
      }
      else {
	file << "position: " << C.shift << " ";
	file << "non-reversed";
      }
      file <<endl;
    }
  }
}
