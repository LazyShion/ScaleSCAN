#include<iostream>
#include<fstream>
#include<map>
#include<vector>
#include<cstdlib>
#include<algorithm>
#include<set>

#define CH printf("check!!\n");

using namespace std;

int main(int argc,char *argv[])
{
  ifstream infile1(argv[1]);
  ofstream outfile(argv[2],ios_base::out | ios_base::binary);
  string str1,str2;
  unsigned int node1,node2;
  
  map<unsigned int , unsigned int> mp2;
  unsigned int node_count;

  node_count = 0;

  unsigned int index_size = 0;

  cout << "create_node_number_map" << endl;
  while(getline(infile1,str1,'\t')){
    getline(infile1,str2);
     node1 = atoi(str1.c_str());
     node2 = atoi(str2.c_str());
     if(node1 != node2){
        if(mp2.find(node1) == mp2.end()){
	  mp2[node1] = node_count;
	  node_count++;

	}

	index_size++;
	if(mp2.find(node2) == mp2.end()){
	  mp2[node2] = node_count;
	  node_count++;
	}
     }
      
  }

  
  ifstream infile2(argv[1]);
  unsigned int *edge_list = new unsigned int[index_size*2];

  unsigned int edge_count = 0;
  cout << "create_edge_list" << endl;
  while(getline(infile2,str1,'\t')){
     getline(infile2,str2);
     node1 = atoi(str1.c_str());
     node2 = atoi(str2.c_str());
     if(node1 != node2){
       edge_list[edge_count] = mp2[node1];
       edge_count++;
       edge_list[edge_count] = mp2[node2];
       edge_count++;
     }
  }
  mp2.clear();
  unsigned int c = 0;
  unsigned int *node = new unsigned int[node_count + 1];
  node[0] = 0;
  ofstream log("logfile.txt",ios_base::out | ios_base::binary);
  cout << "create_node_array" << endl;
  while(c < node_count){
    map<unsigned int,set<unsigned int>> mp;
    unsigned int c0 = c;
    unsigned int c2 = c+100000;

    if(c2 >= node_count){
      c2 = node_count;
    }
     
    
    for(unsigned int i = 0;i < index_size*2;i+=2){
      if(c <= edge_list[i] && edge_list[i] < c2){
	mp[edge_list[i]].insert(edge_list[i+1]);
      }
      
      if(c <= edge_list[i+1] && edge_list[i+1] < c2){
	mp[edge_list[i+1]].insert(edge_list[i]);
      }
    }
    
    while(c<c2){
      unsigned int index = mp[c].size();
      node[c+1] = node[c] + index;
      c++;
    }
    
    
    unsigned int edgesize = node[c] - node[c0];
    cout << edgesize << endl;
    unsigned int *edge_log = new unsigned int[edgesize];
    unsigned int n = 0;
    
    // cout << "aaa" << endl;
    while(c0 < c2){
      for(auto itr = mp[c0].begin();itr != mp[c0].end();++itr){
	edge_log[n] = *itr;
	n++;
      }
      c0++;
    }
    
    log.write((const char*)&edge_log[0],sizeof(unsigned int)*edgesize);
    delete[] edge_log;
    mp.clear();
    
  }

  cout << "write_file" << endl;
  unsigned int node_total = node_count + 1;
  unsigned int edge_total = node[node_count];
  outfile.write((const char*)&node_count,sizeof(int));
  outfile.write((const char*)&edge_total,sizeof(int));
  outfile.write((const char*)&node[0],sizeof(int)*(node_total));

  delete[] node;

  ifstream infile3("logfile.txt",ios_base::in | ios_base::binary);
  unsigned int *edge = new unsigned int[edge_total];
  infile3.read((char*)&edge[0],sizeof(int)*(edge_total));

  outfile.write((const char*)(edge),sizeof(int)*(edge_total));
  log.close();
  remove("logfile.txt");
   
}
