// B657 assignment 3 skeleton code, D. Crandall
//
// Compile with: "make"
//
// This skeleton code implements nearest-neighbor classification
// using just matching on raw pixel values, but on subsampled "tiny images" of
// e.g. 20x20 pixels.
//
// It defines an abstract Classifier class, so that all you have to do
// :) to write a new algorithm is to derive a new class from
// Classifier containing your algorithm-specific code
// (i.e. load_model(), train(), and classify() methods) -- see
// NearestNeighbor.h for a prototype.  So in theory, you really
// shouldn't have to modify the code below or the code in Classifier.h
// at all, besides adding an #include and updating the "if" statement
// that checks "algo" below.
//
// See assignment handout for command line and project specifications.
//
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <sys/types.h>
#include <dirent.h>
#include <map>
#include <numeric>
#include <sstream>
#include <fstream>
#include <iterator>
#include <cstdlib>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

// Dataset data structure, set up so that e.g. dataset["bagel"][3] is
// filename of 4th bagel image in the dataset
typedef map<string, vector<string> > Dataset;

#include <Classifier.h>



typedef map<string, vector< vector<SiftDescriptor> > > Data;
typedef map<string,vector< vector<int> > > Data_Histograms;
typedef map<string,vector< vector<double> > > Data_Haar;
typedef map<string,vector< string> > Data_Deep;

typedef map<string,int> class_encode;


#include "HaarFeatures.h"
#include "KMeans.h"
#include "BagOfVisualWords.h"
#include "DeepFeatures.h"
#include "SVM.h"
#define k 50
// Figure out a list of files in a given directory.
//
vector<string> files_in_directory(const string &directory, bool prepend_directory = false)
{
  vector<string> file_list;
  DIR *dir = opendir(directory.c_str());
  if(!dir)
    throw std::string("Can't find directory " + directory);
  
  struct dirent *dirent;
  while ((dirent = readdir(dir))) 
    if(dirent->d_name[0] != '.')
      file_list.push_back((prepend_directory?(directory+"/"):"")+dirent->d_name);

  closedir(dir);
  return file_list;
}

int main(int argc, char **argv)
{
  try {
    if(argc < 3)
      throw string("Insufficent number of arguments");

    string mode = argv[1];
    string algo = argv[2];
	
	class_encode class_encoding ;
	/*
	class_encoding["bagel"]= 1; 
	class_encoding["bread"]= 2;
	class_encoding["brownie"]= 3;
	class_encoding["chickennugget"]= 4;
	class_encoding["churro"]=5;
	class_encoding["croissant"]=6;
	class_encoding["frenchfries"]=7;
	class_encoding["hamburger"]=8;
	class_encoding["hotdog"]=9;
	class_encoding["jambalaya"]=10;
	class_encoding["kungpaochicken"]=11;
	class_encoding["lasagna"]=12;
	class_encoding["muffin"]=13;
	class_encoding["paella"]=14;
	class_encoding["pizza"]=15;
	class_encoding["popcorn"]=16;
	class_encoding["pudding"]=17;
	class_encoding["salad"]=18;
	class_encoding["salmon"]=19;
	class_encoding["scone"]=20;
	class_encoding["spaghetti"]=21;
	class_encoding["sushi"]=22;
	class_encoding["taco"]=23;
	class_encoding["tiramisu"]=24;
	class_encoding["waffle"]=25;
	*/

    // Scan through the "train" or "test" directory (depending on the
    //  mode) and builds a data structure of the image filenames for each class.
    Dataset filenames; 
    vector<string> class_list = files_in_directory(mode);
    for(vector<string>::const_iterator c = class_list.begin(); c != class_list.end(); ++c)
      filenames[*c] = files_in_directory(mode + "/" + *c, true);
	
	//creates map of classes with index from class_list
	for(int i=0;i<class_list.size();i++){
		class_encoding[class_list[i]] = i+1;
	}
    // set up the classifier based on the requested algo
    Classifier *classifier=0;
	classifier = new SVM(class_list);
    if(algo == "baseline"){
		cout << "nn" << endl;
	}else if(algo == "haar"){
			HaarFeatures *haar = new HaarFeatures(1);
			Data_Haar data_haar = haar->get_haar_features(filenames);
			cout << "Writing Features to file Started..." << endl;
			cout << endl;
			ofstream out;
			if(mode == "train") out.open("data_haar.dat");
			else out.open("test_data.dat");
			for(Data_Haar::const_iterator it=data_haar.begin();it!=data_haar.end();++it){
				for(int i=0; i<it->second.size();i++){
					out << class_encoding[it->first] << " ";
					for(int j = 0; j < it->second[i].size(); j++)
					{
						out<<j+1<<":"<<it->second[i][j] << " ";
					}
					//std::copy(it->second[i].begin(), it->second[i].end(), ostream_iterator<double>(out, " "));
					out << endl;
				}
				
			}
			out.close();
			cout << "Writing Features to file Completed!" << endl;
			cout << endl;
		
		
	}else if(algo == "bow"){
		cout << "Creating Bag Of Words features..." << endl;
		if(mode == "train"){
			BagOfVisualWords bovw(k);
			Data_Histograms data_hist = bovw.get_bov(filenames);
			
			ofstream out("data_histogram.dat");
			for(Data_Histograms::const_iterator it=data_hist.begin();it!=data_hist.end();++it){
				for(int i=0; i<it->second.size();i++){
					out << class_encoding[it->first] << " ";
					for(int j = 0; j < it->second[i].size(); j++)
					{
						out<<j+1<<":"<<it->second[i][j]<< " ";
					}
					//std::copy(it->second[i].begin(), it->second[i].end(), ostream_iterator<int>(out, " "));
					out << endl;
				}
				
			}
			out.close();
			cout << endl;
			cout << "Bag Of Words Features Construction Completed!" << endl;
			ofstream out1("vocabulary.txt");
			vector< vector<float> > vocabulary = bovw.get_vocabulary();
			
			for(int i=0; i<vocabulary.size();i++){
				
				std::copy(vocabulary[i].begin(), vocabulary[i].end(), ostream_iterator<float>(out1, " "));
				out1 << endl;
			}
			out1.close();
		}
		
	}else if(algo == "deep"){
		
			DeepFeatures deep_feat;
			Data_Deep deep_features = deep_feat.get_deep_features(filenames);
			cout << "Features Extraction Completed..." << endl;
			ofstream out;
			if(mode == "train") out.open("deep_features.dat");
			else out.open("test_data.dat");
			for(Data_Deep::const_iterator it=deep_features.begin();it!=deep_features.end();++it){
				 // input stringstream to read from the line
				//istream_iterator<int> iterator;
				// create a vector containing inters in the line (left to right)
				//using iterator = std::istream_iterator<int> ;
				//std::vector<int> seq { iterator(stm), iterator() } ;
				for(int i=0;i<it->second.size();i++){
					int d = it->second.size();
					
					cout << "Size:" << d << endl;
					
					out << class_encoding[it->first] << " ";
						
					std::istringstream buf(it->second[i]);
					std::istream_iterator<std::string> beg(buf), end;
					
					std::vector<std::string> tokens(beg, end);
					
					for(int m=0;m<tokens.size();m++){
						float t = strtof(tokens[m].c_str(),0);
						out << m+1 << ":" << t << " ";
					}
					
					out << endl;
					/*
					 istringstream stm(it->second[i]) ;
					 vector<float> seq;
					 float number;
					 while ( stm >> number )
					 seq.push_back( number );
					 
					 for(int j=0; j< seq.size();i++){
					 
					 //std::copy(it->second[i].begin(), it->second[i].end(), ostream_iterator<string>(out, " "));
					 out << j+1 << ":" << seq[j] << endl;
					 }
					 */
					
				}		
			}
			out.close();
			cout << "Writing to file completed" << endl;
		
		
	}else if(algo == "eigen"){
		cout << "eigen" <<endl;
	}
    else
      throw std::string("unknown classifier " + algo);

    // now train or test!
    if(mode == "train"){
      classifier->train(filenames,algo);
	  cout << "train" << endl;
		}
    else if(mode == "test")
      classifier->test(filenames,algo);
    else
      throw std::string("unknown mode!");
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}