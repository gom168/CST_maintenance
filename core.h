#pragma once
#include "Snap.h"
#include "Defines.h"
#include "UnionFind.h"
#include "Utility.h"
#include "Timer.h"
#include "edges.h"
#include "List.h"
#include "gadget/heap.h"
#include "gadget/treap.h"
#include <unordered_set>


using namespace std;



class core
{
#define ___assume(a,b,c) if(!c){auto x = a; a = b; b=x;}
private:
	std::string dir; // file path
	ui n; // number of vertices
	ui m;  // number of edges
	PUNGraph G;  // undirected graph
	ui *cores;  //coreness value of vertices
	ui *degree;  // deg(v)
	ui *peel_sequence;  // degeneracy order
	ui max_core;   // k_max
	std::vector<ui>Vc;  //  
	std::vector<ui>V_star; // vertices set whose coreness changes
	std::vector<std::vector<ui>>core_tree;
	vector<pair<pair<ui, ui>, ui> > spt;   // edges of CST
	edges edg;  // Kruskal order
	queue<pair<PII,ui>>sptt;  // save spt in T
	vector<int> mcd;  // mcd
	ui spt_num;   //number of the |CST(G)|
	ui cnt_tree;   
	PUNGraph TG; // CST graph

	UnionFind *uf_tree;  // union_find

	struct ListNode {
		int rem;
		int ext;
		int prev;
		int next;
	};

	int GetRank(const int v) {
		if (0 == rank_[v]) {
			rank_[v] = tree_.Rank(v);
			garbage_.push_back(v);
		}
		return rank_[v];
	}

	std::vector<int> head_;
	std::vector<int> tail_;
	std::vector<ListNode> node_;
	std::vector<int> mcd_;
	std::vector<int> deg_;
	std::vector<int> rank_;
	std::vector<int> root_;
	std::vector<bool> evicted_;
	std::vector<bool> visited_;
	gadget::Treap tree_;
	gadget::MinHeap heap_;
	std::vector<int> garbage_;

public:

	core(std::string _dir,const int n); //initilaize  _dir is file name, n is the number of the vertices
	~core();
	void read_graph(bool print,bool flag); //read graph from files; print or not print; flag=0 or 1 then the number of vertices begin from 1 or 0 in files

	void ComputeCore(const bool init_idx);  // compute k-core
	std::unordered_set<int>* Insert(const int v1, const int v2);  // coreness maintenace under single-edge insertion
	std::unordered_set<int>* Remove(const int v1, const int v2);  // coreness maintenace under single-edge deletion
	void Check(const std::vector<std::vector<int>>& graph) const;    
	void Keep(const int v, const int K, int & list_t, std::vector<int>& swap);
	void PropagateDismissal(const int K, const int v, std::vector<int>& to_be_clear, std::vector<int>& changed);
	bool before(const int v1, const int v2);

	void core_decompstion(bool print);   // k-core decompstion

	void core_spanning_d(bool flag, ui u, ui v); // static CST recalculate after inserting or deleting (u,v); flag=0 for delete 1 for insert

	void core_spanning(bool print);  // CST initilaize

	void Insert_core_spanning_new(ui u, ui v);  // CST maintenance after inserting (u,v)

	void Delete_core_spanning_new(ui u, ui v);   // CST maintenance after deleting (u,v)

	void Dynamic_maintainence(bool flag, ui u, ui v);   //flag=0/1 for delete and insert  CST maintenance after deleting/inserting (u,v)

	void Dynamic_maintainences(string str, bool flag,bool flag2,bool read_flag);  

	//str:update files; flag=0/1 for dynamic/static; flag2=0/1 for delete/insert; read_flag=0/1 for vertices begin from number 1 or 0

	void Dynamic_maintainences2(string str, bool flag, bool flag2, bool read_flag);
	// stability test
	
	void get_kcores(int k);    // CST for k-core
	long long get_kcores_(int k);  // return time

	void bfs_kcores(int k); //coreness for k-core
	long long bfs_kcores_(int k); // return time

	void connected_kcores(int k); //static method
	long long connected_kcores_(int k); // return time

	void Dynamic_kcores(string str, int flag, bool flag2, bool read_flag,int k1, int k2);   // dynamic k-core 
	//str:update files; flag=0/1 for dynamic/static; flag2=0/1 for delete/insert; read_flag=0/1 for vertices begin from number 1 or 0
	// varying from k1 to k2
private:
	void get_mcd();
};

