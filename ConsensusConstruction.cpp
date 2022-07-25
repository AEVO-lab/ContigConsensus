#include "ConsensusConstruction.hpp"
#include <cstdlib>
#include <iostream>

void ConsensusConstruction::parseFastaFile(const char* fileName, unsigned set_id, bool prioritize)
{
  ifstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  
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
	m_contig[name]=construct_contig(set_id, contig_id_count++, move(name), move(seq),prioritize);
    }

}


void ConsensusConstruction::parseMatches(const char * fileName,
					 unsigned id_set1, unsigned id_set2)
{
  ifstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }

  // insert keys
  //auto &v = get<0>(matches[id_set1][id_set2]);
  
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
	construct_match(s,t,startS-1,endS-1,startT-1,endT-1,score);
      else
	construct_match(t,s,startT-1,endT-1,startS-1,endS-1,score);

    }
  }
}

void ConsensusConstruction::treatFastaDirectory(const char *dirName, bool only_read_names, string prioritize)
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
    if(!only_read_names) parseFastaFile(file.c_str(), id,f_name==prioritize);
    id++;
  }
}
  
void ConsensusConstruction::treatMatchDirectory(const char *dirName, bool display)
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
      parseMatches(file.c_str(), id1->second, id2->second);
    else parseMatches(file.c_str(),id2->second,id1->second);
  }
  
}

void ConsensusConstruction::output_contigs(const char* fileName,bool random_when_tie, bool prioritize)
{
  unsigned id=1;
  
  auto & contigs = get_consensus_result();

  vector<Contig*> v_contigs;

  // TODO optimize sorting
  for(auto &c : contigs)
    v_contigs.push_back(c.get());

  sort(v_contigs.begin(),v_contigs.end(),[](Contig *c1,Contig *c2){return c1->size()<c2->size();});



  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  
  cout << "Writing fasta result in " << fileName << endl;
  for(auto &c : v_contigs){
    // if(c->getName()=="NODE_20075_length_38_cov_96.078949"){
    // 	c->write_contig_debug(cout, id, !random_when_tie, 10);
    // 	exit(EXIT_SUCCESS);
    // }

    //      c->write_contig(cout, id, !random_when_tie,prioritize);
    c->write_contig_with_filter(file, id, !random_when_tie, 10);
    id++;
  }
}

void ConsensusConstruction::output_contig_ordering(const char * fileName, bool prioritize)
{
  map<unsigned,string> r_map;
  for(auto &a : ids)
    r_map[a.second]=a.first;

  auto & contigs = get_consensus_result();
  
  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Writing ordering result in " << fileName << endl;
  for(auto &c : contigs){
    //    if(prioritize &&) TODO fix prioritize
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

void ConsensusConstruction::createMatchDirectory(const char *m_dirName, const char *f_dirName){


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
      // cmd_blast.append(".txt\"");
      cmd_blast.append(".txt\" -outfmt \"6 score qseqid qstart qend qlen sseqid sstart send slen\"");
      system(cmd_blast.c_str());
    }
  }
}
