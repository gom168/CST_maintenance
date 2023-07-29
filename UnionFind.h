#pragma once
#include "Defines.h"

class UnionFind {
private:
	int n; // number of elements
	int *parent;
	int *rank;

public:
	UnionFind(int _n) {
		n = _n;
		parent = rank = nullptr;
	}
	~UnionFind() {
		if (parent != nullptr) {
			delete[] parent;
			parent = nullptr;
		}
		if (rank != nullptr) {
			delete[] rank;
			rank = nullptr;
		}
	}
	int get_n() { return n; }

	void init(int _n = 0) {
		if (_n == 0) _n = n;
		assert(_n <= n);

		if (parent == nullptr) parent = new int[n];
		if (rank == nullptr) rank = new int[n];

		for (int i = 0; i < _n; i++) {
			parent[i] = i;
			rank[i] = 0;
		}
	}

	UnionFind(const UnionFind& other) {
		n = other.n;
		parent = new int[n];
		rank = new int[n];
		for (int i = 0; i < n; i++) {
			parent[i] = other.parent[i];
			rank[i] = other.rank[i];
		}
	}

	void init(int *ids, int _n) {
		assert(_n <= n);

		if (parent == nullptr) parent = new int[n];
		if (rank == nullptr) rank = new int[n];

		for (int i = 0; i < _n; i++) {
			int u = ids[i];
			parent[u] = u;
			rank[u] = 0;
		}
	}

	void add(int u, int v) {
		parent[u] = v;
		rank[u] = 0;
	}

	int UF_find(int u) {
		int res = u;
		while (parent[res] != res) res = parent[res];
		while (parent[u] != res) {
			int tmp = parent[u];
			parent[u] = res;
			u = tmp;
		}
		return res;
	}

	// return the new root of the merged tree
	int UF_union(int u, int v) {
		int tu = UF_find(u);
		int tv = UF_find(v);

		if (tu == tv) return tu;

		int res;
		if (rank[tu] > rank[tv]) {
			res = tu;
			parent[tv] = tu;
		}
		else {
			res = tv;
			parent[tu] = tv;
			if (rank[tu] == rank[tv]) ++rank[tv];
		}
		return res;
	}
};


