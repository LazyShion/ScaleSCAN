#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <set>
#include <atomic> 
#include <cstdio>

#include"graph.h"
#include"unionfind.h"

#include <xmmintrin.h>	//SSE
#include <emmintrin.h>	//SSE  2
#include <pmmintrin.h>	//SSE  3
#include <tmmintrin.h>	//SSSE 3
#include <smmintrin.h>	//SSE  4.1
#include <nmmintrin.h>	//SSE  4.2
#include <immintrin.h>  //AVX~

#include <omp.h> //OpenMP

#define UNCLASSIFIED 0 
#define EXPLORED 1


#define NONMEMBER r 
#define CMEMBER 1 
#define HUB 2     
#define OUTLIER 3 
#define PARA 8 

#define SS 1
#define NONSS 2
#define CALC 3

float EPS = 0.5;
unsigned int MU = 2;
float POW_EPS = EPS*EPS;

const __m256i sj = _mm256_set_epi32(1,1,1,1,0,0,0,0);
const __m256i st = _mm256_set_epi32(3,2,1,0,3,2,1,0);
 

inline void clustering(GRAPH *g, UnionFind *uf);
inline void CheckCore(GRAPH *g, unsigned int node1, unsigned int node2, unsigned int i);
inline double SS_calc(GRAPH *g, unsigned int node1, unsigned int node2);
inline void ClusterCore(GRAPH *g, UnionFind *uf, unsigned int i);
inline void ClusterNoncore(GRAPH *g, UnionFind *uf);
inline bool hub_check_uf(GRAPH *g, UnionFind *uf, unsigned int a);
inline unsigned int BinarySearchEdge(GRAPH *g, unsigned int node1, unsigned int node2);
inline bool pruning_rule(GRAPH *g, unsigned int node1, unsigned int node2);

unsigned int NUMT;

using namespace std;
