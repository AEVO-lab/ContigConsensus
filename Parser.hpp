#pragma once

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


void parseFastaFile(const char* fileName, unsigned set_id, AssemblySet &assembly_sets)
{
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
      if(seq.size()>0)
	v.emplace(make_unique<Contig>(set_id,move(name),move(seq)));
		
    }

#if DEBUG==true
  for(const auto &cont : v)
    cout << *cont.get() << "\n";
#endif

}


void parseMatches(const char * fileName, AssemblySet & assembly_sets, MatchMatrix & matches, unsigned id_set1, unsigned id_set2)
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

    auto a = assembly_sets[id_set1].find(nameS);
    if(a==assembly_sets[id_set1].end()){
      a = assembly_sets[id_set2].find(nameS);
      if(a==assembly_sets[id_set2].end()){
	cout << "error " << nameS << " doesn't exist in "<< id_set1 << endl;
	for(auto &p : assembly_sets[id_set1]){
	  cout << p->getName() << endl;
	  exit(EXIT_FAILURE);
	}
	
	exit(EXIT_FAILURE);
      }
    }
    Contig *s = a->get();
    
    a = assembly_sets[id_set2].find(nameT);
    if(a==assembly_sets[id_set2].end()){
      a = assembly_sets[id_set1].find(nameT);
      if(a==assembly_sets[id_set1].end()){
	cout << "error " << nameT << " doesn't exist" << endl;
	exit(EXIT_FAILURE);
      }
    }
    Contig *t = a->get();


    while(startS != 1 && startS != lengthS && startT !=1 && startT !=lengthT){
      startS = startS<endS ? startS-1 : startS+1;
      startT = startT<endT ? startT-1 : startT+1;
      if(s->getNuc(startS-1)==t->getNuc(startT-1))
	score++;
    }
    while(endS != 1 && endS != lengthS && endT !=1 && endT !=lengthT){
      endS = startS<endS ? endS+1 : endS-1;
      endT = startT<endT ? endT+1 : endT-1;
      if(s->getNuc(endS-1)==t->getNuc(endT-1))
	score++;
    }

    unsigned lengthS = startS<endS ? endS-startS : startS-endS;
    unsigned lengthT = startT<endT ? endT-startT : startT-endT;
    assert(lengthS==lengthT);
    
    if(s->get_set_id()<t->get_set_id())
      v.emplace_back(s,t,startS-1,endS-1,startT-1,endT-1,score);
    else
      v.emplace_back(t,s,startT-1,endT-1,startS-1,endS-1,score);
    
    
  }

}



void treatFastaDirectory(const char *dirName, AssemblySet &assembly_sets, map<string,unsigned> &ids,bool only_read_names)
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
    if(f_name.substr(f_name.find_last_of(".")) != ".fasta")
      continue;
		
    cout << "Read fasta file: " << f_name <<"\n";
    string file = dirName;
    file.append("/");
    file.append(ent->d_name);
    ids[f_name.substr(0,f_name.find_last_of("."))]=id;
    if(!only_read_names) parseFastaFile(file.c_str(), id, assembly_sets);
    id++;
  }
}


void treatMatchDirectory(const char *dirName, AssemblySet &assembly_sets, map<string,unsigned> &ids,MatchMatrix &matches, bool display=true)
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
      parseMatches(file.c_str(), assembly_sets, matches,id1->second,id2->second);
    else parseMatches(file.c_str(), assembly_sets, matches,id2->second,id1->second);
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
      string cmd_blast="blastn -task megablast -query ";
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

void output_contig(const char* fileName, AssemblySet & assembly_sets)
{
  auto & contigs = assembly_sets.begin()->second;
  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Writing fasta result in " << fileName << endl;
  for(auto &c : contigs){
    file << *c;
  }
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
