#pragma once

#include "Snap.h"
#include "Defines.h"
#include "UnionFind.h"
#include "Utility.h"
#include "Timer.h"

using namespace std;

template <typename T>
inline void hash_combine(std::size_t &seed, const T &val) {
	seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
// auxiliary generic functions to create a hash value using a seed
template <typename T> inline void hash_val(std::size_t &seed, const T &val) {
	hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t &seed, const T &val, const Types &... args) {
	hash_combine(seed, val);
	hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
	std::size_t seed = 0;
	hash_val(seed, args...);
	return seed;
}

struct pairHash {
	template <class T1, class T2>
	std::size_t operator()(const std::pair<T1, T2> &p) const {
		return hash_val(p.first, p.second);
	}
};


class edges
{
private:
	vector<ui>id;
	vector<ui>core_num;
	vector<int>k_num;
	vector<ui>rid;
	ui k_max;   // max edge weight
	ui m;
	ui all_edges;
	unordered_map<PII, ui,pairHash>m_edges;// map edges to id

	unordered_map<ui,PII>rm_edges;// map id to edges.

public:
	edges();
	edges(ui _m, ui _all_edges);
	~edges();

	void Insert(ui u, ui v, ui core_);  // insert(u,v) in Kruskal order  core_ = min(core(u),core(v))

	void Delete(ui u, ui v, ui core_); // delete (u,v) in Kruskal order core_ = min(core(u),core(v))

	void Update(ui u, ui v, ui core_old, ui core_new); // update (u,v) in Kruskal order form edge weight core_old to core_new

	inline void set_m(ui _m)
	{
		m = _m;
	}

	inline void set_edges(ui _all_edges)
	{
		all_edges = _all_edges;
	}


	inline int get_corenum(ui _i)
	{
		return core_num[_i];
	}


	inline PII get_rmedges(ui _i)
	{
		return rm_edges[_i];
	}

	inline ui get_rid(int _i)
	{
		return rid[_i];
	}
	inline ui get_id(int _i)
	{
		return id[_i];
	}

	inline ui get_medges(PII _p)
	{
		return m_edges[_p];
	}

	inline void set_medges(pair<PII, ui>_p)
	{
		m_edges.insert(_p);
	}

	inline void set_rmedges(pair<ui, PII>_p)
	{
		rm_edges.insert(_p);
	}

	inline void set_corenum(ui _u)
	{
		core_num.push_back(_u);
	}

	inline void set_id(ui _id)
	{
		id.push_back(_id);
	}

	inline void set_rid(ui _rid)
	{
		rid.push_back(_rid);
	}

	inline void set_knum(ui _knum)
	{
		k_num.push_back(_knum);
	}

	inline void set_kmax(int _kmax)
	{
		k_max = _kmax;
	}

	inline int get_kmax()
	{
		return k_max;
	}
};

