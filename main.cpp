/*
ScaleSCAN for AVX256
Last update: March 4th, 2019
Authors: T. Takahashi and H. Shiokawa
*/

#include "main.h"

#define CH printf("ch!\n");
#define CHT printf("KTKTKT\n");

unsigned int counter;
int main(int argc,char *argv[])
{
  counter = 0;
  unsigned int cluster_num = 0;
  unsigned int hub_num = 0;
  unsigned int out_num = 0;
  
  EPS = (float)atoi(argv[3]);
  MU = (unsigned int)atoi(argv[4]);
  POW_EPS = EPS*EPS;

  double t0 = omp_get_wtime();
  NUMT = atoi(argv[2]);
  GRAPH g(argv[1]);
  UnionFind uf(g.nodemax);
  
  double t1 = omp_get_wtime();

  clustering(&g, &uf);
  cout << "n =" << g.nodemax << endl;
  cout << "m =" << g.edgemax << endl;
 

}

inline void CheckCore(GRAPH *g, unsigned int node1,unsigned int node2,unsigned int i)
{
  if((*g).ed[node1] >= MU && (*g).sd[node1] < MU){
    unsigned int lower;
    if(node1 > node2){
      lower = BinarySearchEdge(g,node2,node1);
    }else{
      lower = i;
    }
    if(__sync_bool_compare_and_swap(&((*g).similarity[lower]),UNCLASSIFIED,CALC)){
      double ss = SS_calc(g,node1,node2);
      unsigned int mid = BinarySearchEdge(g,node2,node1);
	if(ss >= EPS){

	  (*g).similarity[i] = SS;
	  (*g).similarity[mid] = SS;
#pragma omp atomic
	  (*g).sd[node1]++;
#pragma omp atomic
	  (*g).sd[node2]++;
	  
	}else{
	  
	  (*g).similarity[i] = NONSS;
	  (*g).similarity[mid] = NONSS;
	  
#pragma omp atomic
	  (*g).ed[node1]--;
	  
#pragma omp atomic
	  (*g).ed[node2]--;
	
	}
      
      
	
    }
  }
  
}

inline bool pruning_rule(GRAPH *g,unsigned int node1, unsigned int node2){
  double x = (*g).node[node1+1] - (*g).node[node1] + 1 ;
  double y = (*g).node[node2+1] - (*g).node[node2] + 1;
  

  if(y<x*(POW_EPS) || x < y*(POW_EPS)){
    return true;
  }else{
    return false;
  }
}

inline double SS_calc(GRAPH *g,unsigned int node1,unsigned int node2){
  unsigned int common_node = 2;
  unsigned int node1size = ((*g).node[node1+1] - (*g).node[node1]) + 1;
  unsigned int node2size = ((*g).node[node2+1] - (*g).node[node2]) + 1;

  if(pruning_rule(g,node1,node2))return 0;
  if(2.0/sqrt(node1size * node2size)>=EPS)return EPS;

  unsigned int j,t,j2,t2;
  unsigned int to,jo;
  unsigned int cnv[PARA] = {0,0,0,0,0,0,0,0};
  __m256i ssecountplus = _mm256_set1_epi32(1);
  __m256i ssecnv = _mm256_load_si256((__m256i*)(cnv));
 
  if((*g).node[node1 + 1] - (*g).node[node1] < (*g).node[node2 + 1] - (*g).node[node2]){
    j=(*g).node[node1];
    t=(*g).node[node2];
    j2=(*g).node[node1+1];
    t2 =(*g).node[node2+1];
  }else{
    t=(*g).node[node1];
    j=(*g).node[node2];
    t2=(*g).node[node1 + 1];
    j2 =(*g).node[node2 + 1];
  }
  unsigned int size1 = (t2-t)/(j2-j);

  if(size1 > 2){
   	__m256i jnode = _mm256_set1_epi32((*g).edge[j]);
	__m256i tnode = _mm256_loadu_si256((__m256i*)((*g).edge + t));
	unsigned int vsize = ((t2 - t) / PARA)*PARA;
	unsigned int to = t;

	while(j<j2 && t<to+vsize){
	  __m256i mask =  _mm256_cmpeq_epi32(jnode,tnode);
	  mask = _mm256_and_si256(ssecountplus,mask);
	  ssecnv = _mm256_add_epi32(ssecnv,mask);
	
	  if((*g).edge[j] > (*g).edge[t+7]){
	    t+=PARA;
	    tnode = _mm256_loadu_si256((__m256i*)((*g).edge + t));
	  } else {
	    j++;
	    jnode = _mm256_set1_epi32((*g).edge[j]);    
	  }
	}
	_mm256_store_si256((__m256i*)cnv,ssecnv);
	
	for(unsigned int cnvplus=0;cnvplus < PARA;cnvplus++){
	  common_node += cnv[cnvplus];
	}
	
	t = (to+vsize);
	
  }else{
    unsigned int jsize = ((j2-j)/2)*2+j;
    unsigned int tsize = ((t2-t)/4)*4+t;
    __m256i jnode,tnode;
    
    jnode = _mm256_loadu_si256((__m256i*)((*g).edge + j));
    tnode = _mm256_loadu_si256((__m256i*)((*g).edge + t));
    
    while(j<jsize && t<tsize){
      
      __m256i jnodeA = _mm256_permutevar8x32_epi32(jnode,sj);
      __m256i tnodeA = _mm256_permutevar8x32_epi32(tnode,st);
      __m256i mask =  _mm256_cmpeq_epi32(jnodeA,tnodeA);
      mask = _mm256_and_si256(ssecountplus,mask);
      ssecnv = _mm256_add_epi32(ssecnv,mask);
      
      if((*g).edge[j+1] == (*g).edge[t+3]){
	j+=2;
	t+=4;
	jnode = _mm256_loadu_si256((__m256i*)((*g).edge + j));
	tnode = _mm256_loadu_si256((__m256i*)((*g).edge + t));
      }else if((*g).edge[j+1] > (*g).edge[t+3]){
	t+=4;
	
	tnode = _mm256_loadu_si256((__m256i*)((*g).edge + t));
	tnodeA = _mm256_permutevar8x32_epi32(tnode,st);
      }else{
	j+=2;
	jnode = _mm256_loadu_si256((__m256i*)((*g).edge + j));
      }
    }

    _mm256_store_si256((__m256i*)cnv,ssecnv);
    
    for(unsigned int cnvplus=0;cnvplus < PARA;cnvplus++){
      common_node += cnv[cnvplus];
    }
    
    
    if(j>=jsize){
      j = jsize;
    }else{
      t = tsize;
    }
  }

  while(j<j2 && t<t2){
    if((*g).edge[j] == (*g).edge[t]){
      common_node++;
      t++;
      j++;
    }else if((*g).edge[j] > (*g).edge[t]){
      t++;
    }else{
      j++; 
    }
  }
    
  return common_node/sqrt(node1size * node2size);
}

inline unsigned int BinarySearchEdge(GRAPH *g,unsigned int node1, unsigned int node2)
{
  unsigned int offset = (*g).node[node1];
  unsigned int low = 0;
  unsigned int high = (*g).node[node1 + 1]-(*g).node[node1];
  unsigned int mid;
  
  while(high >= low){
    mid = (low + high)/2;
    if((*g).edge[offset+mid] == node2){
      return offset+mid;
    }else if((*g).edge[offset+mid] > node2){
      high = mid-1;
    }else{
      low = mid + 1;
    }
  }
  return -1;
}

inline void ClusterCore(GRAPH *g,UnionFind *uf,unsigned int i)
{
  unsigned int x = (*g).edgef[i];
  unsigned int y = (*g).edge[i];
  //CH;
  if((*g).sd[x] >= MU && (*g).sd[y] >= MU){
    unsigned int rx;
    unsigned int ry;
    rx = (*uf).root(x);
    ry = (*uf).root(y);
    if(rx < ry){
	   unsigned int index = rx;
	   rx = ry;
	   ry = index;
    }
   
    if(rx != ry){
     if((*g).similarity[i] == SS){

       while(!__sync_bool_compare_and_swap(&((*uf).data[rx]),rx,ry)){
	 rx = (*uf).root(x);
	 ry = (*uf).root(y);
	
	 if(rx == ry)break;	
	 if(rx < ry){
	   unsigned int index = rx;
	   rx = ry;
	   ry = index;
	 }
       }
     }
     
     if((*g).similarity[i] == UNCLASSIFIED){
        unsigned int mid = BinarySearchEdge(g,y,x);

	 unsigned int lower;
	 if(x > y){
	   lower = mid;
	 }else{
	   lower = i;
	 }

	 if(__sync_bool_compare_and_swap(&((*g).similarity[lower]),UNCLASSIFIED,CALC)){
	   double ss = SS_calc(g,x,y);
	   if(ss >= EPS){
	     (*g).similarity[i] = SS;
	     (*g).similarity[mid] = SS;
	     
	     while(!__sync_bool_compare_and_swap(&((*uf).data[rx]),rx,ry)){
	       rx = (*uf).root(x);
	       ry = (*uf).root(y);
	       
	       if(rx == ry)break;
	       
	       if(rx < ry){
		 unsigned int index = rx;
		 rx = ry;
		 ry = index;
	       }
	     }
	     
	   }else{
	     (*g).similarity[i] = NONSS;
	     (*g).similarity[mid] = NONSS;
	   }
	 }

     }
    }
  }
  
}

inline void ClusterNoncore(GRAPH *g,UnionFind *uf)
{
  
#pragma omp parallel for num_threads(NUMT) schedule(dynamic ,50)
  for(unsigned int i = 0;i < (*g).nodemax;i++){
    for(unsigned int j = (*g).node[i];j < (*g).node[i+1];j++){
      unsigned int y = (*g).edge[j];
      if((*g).sd[i] >= MU || (*g).sd[y] < MU || (*g).included_cluster[i] >= (MU-1))continue;
      unsigned int ry = (*uf).root(y);
      unsigned int flag = 0;
      for(unsigned int n = 0;n < (*g).included_cluster[i];n++){
	if((*g).cluset[i*(MU-1) + n] == ry)flag = 1;
      }
      if(flag == 0){
	if((*g).similarity[j] == SS){
	  (*g).cluset[i*(MU-1)+(*g).included_cluster[i]] = ry;
	  (*g).included_cluster[i]++;
	}if((*g).similarity[j] == UNCLASSIFIED){
	  double ss = SS_calc(g,i,y);
	  if(ss >= EPS){
	    (*g).cluset[i*(MU-1)+(*g).included_cluster[i]] = ry;
	    (*g).included_cluster[i]++;
	  }
	}
	
      }
    }
  }
}

inline bool hub_check_uf(GRAPH *g,UnionFind *uf,unsigned int a)
{
  unsigned int counter = 0;
  set<int> c;
  
  for(unsigned int i=(*g).node[a];i<(*g).node[a+1];i++){
    unsigned int j = (*g).edge[i];
    if((*g).included_cluster[j] == 0 && (*g).sd[j] < MU)continue;
    if((*g).sd[j] >= MU){
    c.insert((*uf).root(j));
    if(c.size() >= 2)return true;
    }else if((*g).included_cluster[j] != 0){
      for(unsigned int n = 0;n<(*g).included_cluster[j];n++){
	c.insert((*g).cluset[j*(MU-1)+n]);
	if(c.size() >= 2)return true;
      }
    }
  }

  return false;
}

inline void clustering(GRAPH *g, UnionFind *uf){
  double t2 = omp_get_wtime();

  #pragma omp parallel for num_threads(NUMT) schedule(guided, 10000)
  for(unsigned int i = 0;i < (*g).edgemax;i++){
    CheckCore(g, (*g).edgef[i], (*g).edge[i],i);
  }

  double t3 = omp_get_wtime();

  #pragma omp parallel for num_threads(NUMT) schedule(guided, 10000)
  for(unsigned int i = 0;i < (*g).edgemax;i++){
    ClusterCore(g, uf,i);
  }
  ClusterNoncore(g, uf);
  
  double t4 = omp_get_wtime();
  cout << (double)(t4-t2)<< endl;
  
}
