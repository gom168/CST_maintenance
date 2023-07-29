#pragma once
#include "Snap.h"
#include "Defines.h"

class MSTruss
{
private:
	std::string dir; // input graph directory
	int n; // nodes of the graph
	int m;  // hyperedges of the graph
	PUNGraph G;
	std::map<std::pair<int, int>, int> eid;
	std::map<int, std::pair<int, int>> rid;
	int *truss_ness;

	PUNGraph H;
	int add_edges;


public:
	MSTruss(std::string dir);
	~MSTruss();
	void read_graph();
	void addEdges(std::string _dir);
	void PeelTruss();
	void compute_pretruss1(int vid, int wid,int truss_copy[]);
	void compute_pretruss2(int vid, int wid,int truss_copy[]);

	void Incremental_truss_maintenance(std::string _dir);
	void IncrementalTraversal(std::multimap<int, std::pair<int, int>,std::greater<int>> m_p,std::map<std::pair<int, int>,int>vis);

};

