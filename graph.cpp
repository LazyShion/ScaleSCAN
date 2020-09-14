#include"graph.h"

#define MYU 2
#define UNCLASSIFIED 0 

GRAPH::GRAPH(char *infilename)
{

  ifstream infile(infilename,ios_base::in | ios_base::binary);

  infile.read((char*)&nodemax, 4);
  infile.read((char*)&edgemax, 4);

  label = new unsigned int[nodemax];
  node = new unsigned int[nodemax+1];
  edge = new unsigned int[edgemax];
  similarity = new unsigned int[edgemax];
  sd = new unsigned int[nodemax];
  ed = new unsigned int[nodemax];
  edgef = new unsigned int[edgemax];
  cluset = new unsigned int[nodemax*(MYU-1)];
  included_cluster = new unsigned int[nodemax];
  
  infile.read((char*)&node[0],sizeof(unsigned int)*(nodemax+1));

  unsigned   int q = 0;
  for(unsigned int i = 0;i<nodemax;i++){
    label[i] = UNCLASSIFIED;
    sd[i] = 0;
    ed[i] = node[i+1]-node[i];
    included_cluster[i] = 0;
    for(unsigned int n=0;n<MYU-1;n++)cluset[(MYU-1)*i + n] = -1;
    for(unsigned int n=node[i];n<node[i+1];n++){
      edgef[q]=i;
      q++;
    }
  }
   infile.read((char*)&edge[0],sizeof(unsigned int)*edgemax);

   for(unsigned int i = 0;i <edgemax;i++){
     similarity[i] = UNCLASSIFIED;
   }
   
}

void GRAPH::power_low_output()
{
  ofstream of("power_low.txt",ios::app);
  unsigned int maxdegree = 0;
  unsigned int *degree;
  degree = new unsigned int[nodemax];
  for(unsigned int i=0;i<nodemax;i++){
    degree[i]=0;
  }

  for(unsigned int i=0;i<nodemax;i++){
    degree[node[i+1]-node[i]]  += 1;
    if(maxdegree < node[i+1]-node[i])maxdegree = node[i+1]-node[i];
  }

  for(unsigned int i=0;i<=maxdegree;i++){
    of << degree[i] << endl;
  }
  
}
