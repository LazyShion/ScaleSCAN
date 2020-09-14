
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<cstdlib>
#include<time.h>
#include<unordered_set>

using namespace std;

struct GRAPH{
  unsigned int nodemax;
  unsigned int edgemax;
  unsigned int *node;  
  unsigned int *edge;  
  unsigned int *edgef; 
  unsigned int *similarity;
  unsigned int *sd;   
  unsigned int *ed;   
  unsigned int *label;
  unsigned int *included_cluster; 
  unsigned int *cluset; 
  GRAPH(char *infilename);
  void power_low_output();
};


