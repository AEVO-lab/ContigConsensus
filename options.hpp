#include <iostream>
#include <unistd.h>
#include <getopt.h>

using namespace std;

struct  Options {

  enum Mode {
    greedy_pairwise, approx, greedy_general
  } ;
  
  string matches_directory, contig_directory, output_fast, output_cons, output_folder;
  bool generate_matches, allow_overlap, random_when_tie, pairwise;
  int opt;
  float percentage_of_match;
  Mode mode;
  string prioritize;

  Options(int argc, char * argv[]) : generate_matches(false), matches_directory("-1"), contig_directory("-1"), output_fast("output.fasta"), output_cons("output.cons"), output_folder("results"), allow_overlap(false), random_when_tie(false), pairwise(false), percentage_of_match(70.0), mode(greedy_general), prioritize("-1")
  {
    option long_arg[] = {
      {"matches", required_argument, nullptr, 'm'},
      {"contigs", required_argument,nullptr, 'c'},
      {"fastaout", required_argument, nullptr, 'f'},
      {"orderout", required_argument,nullptr, 'o'},
      {"folderout", required_argument, nullptr, 'd'},
      {"method", required_argument, nullptr, 'a'},
      {"percentage",required_argument, nullptr, 'e'},
      {"prioritize", required_argument, nullptr, 'b'},
      {"pairwise", no_argument, nullptr,'p'},
      {"genmatch", no_argument, nullptr,'g'},
      {"random", no_argument, nullptr, 'r'},
      {"help", no_argument, nullptr, 'h'}
    };

    string tmp;
    while((opt = getopt_long(argc, argv, "m:c:f:o:d:a:e:b:pgrh",long_arg,nullptr)) !=-1) {
      switch (opt)
	{
	case 'm':
	  matches_directory=optarg;
	  break;
	case 'c':
	  contig_directory=optarg;
	  break;
	case 'f':
	  output_fast=optarg;
	  break;
	case 'o': 
	  output_cons=optarg;
	  break;
	case 'd': 
	  output_folder=optarg;
	  break;
	case 'a':
	  tmp=optarg;
	  if(tmp=="greedy_pairwise")
	    mode=greedy_pairwise;
	  else if(tmp=="approx")
	    mode=approx;
	  else if(tmp=="greedy_general")
	    mode=greedy_general;
	  else {
	    cout << "Error, the mode is unspecified" << endl;
	    cout << tmp << endl;
	    exit(EXIT_SUCCESS);
	  }
	  break;
	case 'e':
	  percentage_of_match=stof(optarg);
	  break;
	case 'p':
	  pairwise=true;
	  break;
	case 'g':
	  generate_matches=true;
	  break;
	case 'b':
	  prioritize=optarg;
	  break;

	case 'h':
	default:
	  cout << "Available options:\n";
	  cout << "--matches     -m:   directory containing the match files\n";
          cout << "                    Input files must be generated by blastn\n";
          cout << "                    using the command: \n\n";
	  cout << "                    \t blastn -task megablast -query file1.fasta\n";
	  cout << "                    \t -subject file2.fasta -ungapped -out inputFileName\n";
	  cout << "                    \t -outfmt \"6 score qseqid qstart qend qlen sseqid sstart send slen\"\n\n";

	  cout << "--contigs     -c:   directory containing the fasta files\n";
	  cout << "--fastaout    -f:   result fasta file\n";
	  cout << "--orderout    -o:   file containing the relative positions of the input contigs\n";
	  cout << "--folderout   -d:   result directory when pairwise option is selected\n";
	  cout << "--method      -a:   method used to construct the consensus\n";
	  cout << "                    Available options:\n";
	  cout << "                    greedy_pairwise, approx, greedy_general (default)\n";
	  cout << "--percentage  -e:   percentage minimum of correspondance in a match\n";
	  cout << "--pairwise    -p:   generate a result for each pair of input fasta file\n";
	  cout << "--genmatch    -g:   generate matches files with blastn\n";
	  cout << "--prioritize    :   todo\n";
	  cout << endl;
	  exit(EXIT_SUCCESS);
	  break;
	}
    }
    if(contig_directory=="-1") {
      cout << "Error, fasta directory not provided\n";
      exit(EXIT_FAILURE);
    }
    if(matches_directory=="-1") {
      cout << "Error, match directory not provided\n";
      exit(EXIT_FAILURE);
    }
  }
};
