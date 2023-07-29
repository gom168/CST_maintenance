#include "core.h"
#define DIRNAME "graphviz"
#include "defs.h"

using namespace std;

core::core(string _dir, const int n_) :tree_(n_), heap_(n_)
{
	dir = _dir;
	n = 0;
	m = 0;
	max_core = 0;
	G = TUNGraph::New();
	spt_num = 0;
	TG = TUNGraph::New();
	head_ = std::vector<int>(n_, -1);
	tail_ = std::vector<int>(n_, -1);
	node_ = std::vector<ListNode>(n_ + 1);
	mcd_ = std::vector<int>(n_, 0);
	deg_ = std::vector<int>(n_, 0);
	rank_ = std::vector<int>(n_, 0);
	root_ = std::vector<int>(n_, n_);
	visited_ = std::vector<bool>(n_, false);
	evicted_ = std::vector<bool>(n_, false);
}

void core::read_graph(bool print,bool flag)
{
	printf("# Start reading graph, Require files  \"corresponding_graph.txt\"\n");

	ifstream ifs;

	if (print)
	{
		ifs.open("3.txt", ifstream::in);
	}
	else
	{
		ifs.open(dir.c_str(), ifstream::in);
	}

	std::string s;
	getline(ifs, s);
	std::stringstream k(s);

	//input number of vertices and edges
	k >> n, k >> m;
	degree = new ui[n + 1];
	memset(degree, 0, sizeof(ui)*(n + 1));
	peel_sequence = new ui[n + 1];
	for (int i = 1; i <= n; i++)
	{
		G->AddNode(i);
		TG->AddNode(i);
	}
	ui edgeid = 0;

	while (getline(ifs, s))
	{
		std::stringstream input(s);
		int x = 0;
		int y = 0;
		input >> x, input >> y;
		if(flag==1)x++, y++;
		if (G->IsEdge(x, y))continue;
		G->AddEdge(x, y);
		degree[x-1]++;
		degree[y-1]++;
	}
	m = G->GetEdges();
	cores = new ui[n + 1];

	/*for (int i = 1; i <= n; i++)
	{
		cout << "degree[i]:" << degree[i] << endl;
	}*/

	ifs.close();
	cout << "read Graph oK" << endl;
}

void core::core_decompstion(bool print)
{
	max_core = 0;
	ui *tem_degree = new ui[n + 1];
	memcpy(tem_degree, degree, sizeof(ui)*(n + 1));
	for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++)
	{
		int v = NI.GetId();
		max_core = max(max_core, degree[v-1]);
	}
	int cnt = 0;

	ui *bin = new ui[n+1];
	
	for (int i = 0; i <= n; i++)
	{
		bin[i] = 0;
	}
	for (int i = 0; i < n; i++)
	{
		++bin[degree[i]];
	}

	ui start = 0;
	for (int i = 0; i <= max_core; i++)
	{
		ui num = bin[i];
		bin[i] = start;
		start += num;
	}

	ui *pos = new ui[n];
	ui *vert = new ui[n];
	for (ui i = 0; i < n; i++)
	{
		ui v = i;
		pos[v] = bin[degree[v]];
		vert[pos[v]] = v;
		++bin[degree[v]];
	}

	for (ui i = max_core; i > 0; i--)
	{
		bin[i] = bin[i - 1];
	}
	bin[0] = 0;

	//ui ctt = 0;
	//memset(cores, 0, sizeof(ui)*(n + 1));
	for (ui i = 0; i < n; i++)
	{
		ui v = vert[i];
		//cout << v << endl;
		cores[v] = degree[v];
		max_core = max(max_core, cores[v]);
		peel_sequence[cnt++] = v;
		TUNGraph::TNodeI NV = G->GetNI(v+1);
		for (int e = 0; e < NV.GetOutDeg(); e++)
		{
			//ctt++;
			ui u = NV.GetOutNId(e);
			u = u - 1;
			if (degree[u] <= degree[v])continue;
			
			ui pu = pos[u];
			ui pw = bin[degree[u]];
			ui w = vert[pw];
			
			if (w!=u)
			{
				pos[u] = pw;
				vert[pu] = w;
				pos[w] = pu;
				vert[pw] = u;
			}
			++bin[degree[u]];
			--degree[u];

		}

	}
	//cout << max_core << endl;
	memcpy(degree, tem_degree, sizeof(ui)*(n + 1));

	delete[] bin;
	delete[] pos;
	delete[] vert;
	delete[] tem_degree;

}

void core::ComputeCore(const bool init_idx)
{
	ui *deg = new ui[n + 1];
	memcpy(deg, degree, sizeof(ui)*(n + 1));
	//auto& deg = degree;
	int max_deg = 0;
	for (int i = 0; i < n; ++i) {
		if (deg[i] > max_deg) {
			max_deg = deg[i];
		}
	}
	std::vector<int> bin(max_deg + 1, 0);
	for (int i = 0; i < n; ++i) {
		++bin[deg[i]];
	}
	int start = 0;
	for (int i = 0; i <= max_deg; ++i) {
		int temp = bin[i];
		bin[i] = start;
		start += temp;
	}
	std::vector<int> vert(n);
	std::vector<int> pos(n);
	for (int i = 0; i < n; ++i) {
		pos[i] = bin[deg[i]];
		vert[pos[i]] = i;
		++bin[deg[i]];
	}
	for (int i = max_deg; i > 0; --i) {
		bin[i] = bin[i - 1];
	}
	bin[0] = 0;
	int k = 0;
	int cnt = 0;
	auto vis = std::vector<bool>(n, false);
	for (int i = 0; i < n; ++i) {
		const int v = vert[i];
		if (deg[v] > k) k = deg[v];
		ASSERT(bin[deg[v]] == i);
		++bin[deg[v]];
		cores[v] = k;
		max_core = max(max_core, cores[v]);
		vis[v] = true;
		int rem = 0;
		peel_sequence[cnt++] = v;
		TUNGraph::TNodeI NV = G->GetNI(v + 1);
		for (int e = 0; e < NV.GetOutDeg(); e++)
		{
			ui u = NV.GetOutNId(e);
			u = u - 1;
			if (vis[u]) continue;
			++rem;
			const int pw = bin[deg[u]];
			const int pu = pos[u];
			if (pw != pu) {
				const int w = vert[pw];
				vert[pu] = w;
				pos[w] = pu;
				vert[pw] = u;
				pos[u] = pw;
			}
			++bin[deg[u]];
			--deg[u];
			if (pos[u] == i + 1) {
				bin[deg[u]] = pos[u];
			}
		}
		if (init_idx) {
			node_[v].rem = rem;
			if (head_[k] == -1) {
				node_[v].prev = node_[v].next = n;
				head_[k] = tail_[k] = v;
			}
			else {
				node_[v].next = n;
				node_[v].prev = tail_[k];
				node_[tail_[k]].next = v;
				tail_[k] = v;
			}
			tree_.Insert(v, false, root_[k]);
		}
	}
	{
		/*			printf("node_[i].rem: ");
					for (int i = 0; i < n_; i++) {
						printf("%d ", node_[i].rem);
					}
					printf("\n");*/
	}
	if (init_idx) {
		for (int v = 0; v < n; ++v) {
			mcd_[v] = 0;
			TUNGraph::TNodeI NV = G->GetNI(v + 1);
			for (int e = 0; e < NV.GetOutDeg(); e++)
			{
				ui u = NV.GetOutNId(e);
				u = u - 1;
				if (cores[u] >= cores[v]) {
					++mcd_[v];
				}
			}
		}
	}
	/*
	for (int i = 0; i < n; i++)
	{
		cout << "i:" << i << ' ' << peel_sequence[i] << endl;
	}*/
	cout << "Compute core is OK" << endl;
	cout << max_core << endl;
}

bool core::before(const int v1, const int v2) {
	return ((cores[v1] == cores[v2] &&
		tree_.Rank(v1) < tree_.Rank(v2)) ||
		cores[v1] < cores[v2]);
}
std::unordered_set<int>* core::Insert(const int v1, const int v2) {
	// insert the edge
	// update mcd
	G->AddEdge(v1 + 1, v2 + 1);
	m++;
	if (cores[v1] <= cores[v2]) ++mcd_[v1];
	if (cores[v2] <= cores[v1]) ++mcd_[v2];
	// the source node and the current core number
	int src = v1;
	const int K = cores[v1] <= cores[v2] ? cores[v1] : cores[v2];
	if ((cores[v1] == cores[v2] &&
		tree_.Rank(v1) > tree_.Rank(v2)) ||
		cores[v1] > cores[v2]) {
		src = v2;
	}
	// update core number
	++node_[src].rem;
	// there is no need to update the core numbers
	if (node_[src].rem <= K) {
		return new std::unordered_set<int>();;
	}
	// preparing the heap
	heap_.Insert(GetRank(src), src);
	//
	std::vector<int> swap;
	// the set of vertices, denoted as A, that doesn't need to be updated
	int list_h = -1, list_t = -1;
	for (int cur = head_[K]; n != cur; ) {
		if (heap_.Empty() || (node_[cur].ext == 0 && node_[cur].rem <= K)) {
			const int start = cur;
			const int end = heap_.Empty() ? tail_[K] : node_[heap_.Top().val].prev;
			// advance the cur pointer
			cur = node_[end].next;
			// remove this sub-list and reinsert it into A
			node_[node_[start].prev].next = node_[end].next;
			node_[node_[end].next].prev = node_[start].prev;
			node_[start].prev = n;
			node_[end].next = n;
			if (-1 == list_h) {
				list_h = start;
				list_t = end;
			}
			else {
				node_[start].prev = list_t;
				node_[list_t].next = start;
				list_t = end;
			}
			continue;
		}
		// update the heap
		// invariant: heap.Top().val == cur
		ASSERT(heap_.Top().val == cur);
		heap_.Delete(heap_.Top().key);
		// deal with cur
		const int next = node_[cur].next;
		const int cur_deg = node_[cur].ext + node_[cur].rem;
		if (likely(cur_deg <= K)) {
			// insert into A
			node_[node_[cur].prev].next = node_[cur].next;
			node_[node_[cur].next].prev = node_[cur].prev;
			if (likely(-1 != list_h)) {
				node_[cur].next = n;
				node_[cur].prev = list_t;
				node_[list_t].next = cur;
				list_t = cur;
			}
			else {
				node_[cur].prev = node_[cur].next = n;
				list_h = list_t = cur;
			}
			node_[cur].rem = cur_deg;
			node_[cur].ext = 0;
			Keep(cur, K, list_t, swap);
		}
		else {
			// cur is temporarily marked as evicted, i.e.,
			// its core number may be updated finally
			evicted_[cur] = true;
			TUNGraph::TNodeI NV = G->GetNI(cur + 1);
			for (int e = 0; e < NV.GetOutDeg(); e++)
			{
				ui u = NV.GetOutNId(e);
				u = u - 1;
				if (cores[u] == cores[cur] && GetRank(u) > rank_[cur]) {
					++node_[u].ext;
					if (!heap_.Contains(rank_[u])) {
						heap_.Insert(rank_[u], u);
					}
				}
			}
		}
		cur = next;
	}
	ASSERT(heap_.Empty());
	head_[K] = list_h;
	tail_[K] = list_t;
	for (const int v : swap) {
		tree_.Delete(v, root_[K]);
		tree_.InsertAfter(v, node_[v].prev, root_[K]);
	}
	// cope with those vertices whose core need to be updated
	std::unordered_set<int>* V_star = new std::unordered_set<int>();
	if (evicted_[src]) {
		auto tail = -1; // tail
		for (auto v = src; n != v; v = node_[v].next) {
			++cores[v];
			V_star->insert(v);
			node_[v].ext = 0;
			tail = v;
			// update mcd
			TUNGraph::TNodeI NV = G->GetNI(v + 1);
			for (int e = 0; e < NV.GetOutDeg(); e++)
			{
				ui u = NV.GetOutNId(e);
				u = u - 1;
				if (evicted_[u]) continue;
				if (K + 1 == cores[u]) {
					++mcd_[u];
				}
				else if (K == cores[u]) {
					--mcd_[v];
				}
			}
			// remove from the current tree
			tree_.Delete(v, root_[K]);
		}
		for (auto v = tail; n != v; v = node_[v].prev) {
			evicted_[v] = false;
			tree_.Insert(v, true, root_[K + 1]);
		}
		// merge list
		if (-1 == head_[K + 1]) {
			head_[K + 1] = src;
			tail_[K + 1] = tail;
		}
		else {
			node_[head_[K + 1]].prev = tail;
			node_[tail].next = head_[K + 1];
			head_[K + 1] = src;
		}
	}
	for (const int v : garbage_) rank_[v] = 0;
	garbage_.clear();

	return V_star;
}
std::unordered_set<int>* core::Remove(const int v1, const int v2)
{
	G->DelEdge(v1 + 1, v2 + 1);
	m--;
	if (cores[v1] <= cores[v2]) --mcd_[v1];
	if (cores[v2] <= cores[v1]) --mcd_[v2];
	// set the root and core number
	const int root = cores[v1] <= cores[v2] ? v1 : v2;
	const int K = cores[root];
	// update rem
	if (cores[v1] == cores[v2]) {
		if (tree_.Rank(v1) > tree_.Rank(v2)) {
			--node_[v2].rem;
		}
		else {
			--node_[v1].rem;
		}
	}
	else {
		--node_[root].rem;
	}
	// update cores
	std::vector<int> to_be_clear;
	std::vector<int> changed;
	if (cores[v1] != cores[v2]) {
		visited_[root] = true;
		deg_[root] = mcd_[root];
		to_be_clear.push_back(root);
		if (deg_[root] < K) {
			PropagateDismissal(K, root, to_be_clear, changed);
		}
	}
	else {
		visited_[v1] = true;
		deg_[v1] = mcd_[v1];
		to_be_clear.push_back(v1);
		if (deg_[v1] < K) {
			PropagateDismissal(K, v1, to_be_clear, changed);
		}
		if (!visited_[v2]) {
			visited_[v2] = true;
			deg_[v2] = mcd_[v2];
			to_be_clear.push_back(v2);
			if (deg_[v2] < K) {
				PropagateDismissal(K, v2, to_be_clear, changed);
			}
		}
	}
	// clear
	for (const int u : to_be_clear) {
		visited_[u] = false;
		deg_[u] = 0;
	}
	std::unordered_set<int>* V_star = new std::unordered_set<int>();
	for (auto x : changed) {
		V_star->insert(x);
	}
	if (!changed.empty()) {
		while (n!= head_[K] && evicted_[head_[K]]) {
			head_[K] = node_[head_[K]].next;
		}
		while (n != tail_[K] && evicted_[tail_[K]]) {
			tail_[K] = node_[tail_[K]].prev;
		}
		if (n == head_[K]) {
			head_[K] = tail_[K] = -1;
		}
		for (const int v : changed) {
			node_[v].rem = 0;
			TUNGraph::TNodeI NV = G->GetNI(v + 1);
			for (int e = 0; e < NV.GetOutDeg(); e++)
			{
				ui u = NV.GetOutNId(e);
				u = u - 1;
				if (cores[u] == K) {
					--mcd_[u];
					if (!evicted_[u] && GetRank(v) > GetRank(u)) {
						--node_[u].rem;
					}
				}
				else if (cores[u] == K - 1 && !evicted_[u]) {
					++mcd_[v];
				}
				if (cores[u] >= K || (evicted_[u] && !visited_[u])) {
					++node_[v].rem;
				}
			}
			visited_[v] = true;
		}
		for (const auto v : changed) {
			evicted_[v] = false;
			visited_[v] = false;
			tree_.Delete(v, root_[K]);
			tree_.Insert(v, false, root_[K - 1]);
			// remove from current list
			node_[node_[v].next].prev = node_[v].prev;
			node_[node_[v].prev].next = node_[v].next;
			node_[v].next = node_[v].prev = n;
			// merge list
			if (-1 == head_[K - 1]) {
				head_[K - 1] = tail_[K - 1] = v;
			}
			else {
				node_[tail_[K - 1]].next = v;
				node_[v].prev = tail_[K - 1];
				tail_[K - 1] = v;
			}
		}
	}
	for (const int g : garbage_) rank_[g] = 0;
	garbage_.clear();
	return V_star;
}
void core::Check(const std::vector<std::vector<int>>& graph) const {
	for (int v = 0; v < n; ++v) {
		int local_mcd = 0;
		for (const auto u : graph[v]) {
			if (cores[u] >= cores[v]) ++local_mcd;
		}
		ASSERT(mcd_[v] == local_mcd);
		ASSERT(!visited_[v]);
		ASSERT(!evicted_[v]);
		ASSERT(rank_[v] == 0);
		ASSERT(deg_[v] == 0);
	}
	std::vector<bool> vis(n, false);
	for (int v = 0; v < n; ++v) {
		if (vis[v]) continue;
		const int K = cores[v];
		int tail = -1;
		ASSERT(-1 != head_[K]);
		for (int tmp = head_[K]; n != tmp; tmp = node_[tmp].next) {
			ASSERT(!vis[tmp]);
			vis[tmp] = true;
			tail = tmp;
			ASSERT(cores[tmp] == K);
			ASSERT(node_[tmp].ext == 0);
			if (n != node_[tmp].next) {
				ASSERT(node_[node_[tmp].next].prev == tmp);
			}
		}
		ASSERT(tail_[K] == tail);
		ASSERT(node_[head_[K]].prev == n);
		ASSERT(node_[tail_[K]].next == n);

		for (int tmp = head_[K], rid = 0; n != tmp; tmp = node_[tmp].next) {
			ASSERT(tree_.Rank(tmp) == ++rid);
		}
		for (int tmp = head_[K]; n != tmp; tmp = node_[tmp].next) {
			int local = 0;
			for (const auto u : graph[tmp]) {
				if (cores[u] > cores[tmp] ||
					(cores[u] == cores[tmp] &&
						tree_.Rank(u) > tree_.Rank(tmp))) {
					++local;
				}
			}
			ASSERT(local == node_[tmp].rem);
			ASSERT(node_[tmp].rem <= K);
		}
	}
	ASSERT(garbage_.empty());
	ASSERT(heap_.Empty());
}

void core::Keep(const int v, const int K, int& list_t, std::vector<int>& swap) {
	// update
	std::queue<int> bfs;
	TUNGraph::TNodeI NV = G->GetNI(v + 1);
	for (int e = 0; e < NV.GetOutDeg(); e++)
	{
		ui u = NV.GetOutNId(e);
		u = u - 1;
		if (cores[u] == cores[v] && evicted_[u]) {
			--node_[u].rem;
			if (node_[u].rem + node_[u].ext <= K) {
				visited_[u] = true;
				bfs.push(u);
			}
		}
	}
	while (!bfs.empty()) {
		const int u = bfs.front(); bfs.pop();
		visited_[u] = false;
		evicted_[u] = false;
		// insert u into the list
		node_[node_[u].prev].next = node_[u].next;
		node_[node_[u].next].prev = node_[u].prev;
		swap.push_back(u);
		node_[list_t].next = u;
		node_[u].next = n;
		node_[u].prev = list_t;
		node_[u].rem += node_[u].ext;
		node_[u].ext = 0;
		// advance the tail of list
		list_t = u;
		// find more vertices to keep
		TUNGraph::TNodeI NV = G->GetNI(u + 1);
		for (int e = 0; e < NV.GetOutDeg(); e++)
		{
			ui w = NV.GetOutNId(e);
			w = w - 1;
			if (cores[w] != cores[u]) continue;
			if (rank_[w] > rank_[v]) {
				--node_[w].ext;
				if (0 == node_[w].ext) {
					heap_.Delete(rank_[w]);
				}
			}
			else if (rank_[w] > rank_[u] && evicted_[w]) {
				--node_[w].ext;
				if (!visited_[w] && node_[w].ext + node_[w].rem <= K) {
					visited_[w] = true;
					bfs.push(w);
				}
			}
			else if (evicted_[w]) {
				--node_[w].rem;
				if (!visited_[w] && node_[w].ext + node_[w].rem <= K) {
					visited_[w] = true;
					bfs.push(w);
				}
			}
		}
	}
}
void core::PropagateDismissal(const int K, const int v,
	std::vector<int>& to_be_clear,
	std::vector<int>& changed) {
	evicted_[v] = true;
	--cores[v];
	changed.push_back(v);
	TUNGraph::TNodeI NV = G->GetNI(v + 1);
	for (int e = 0; e < NV.GetOutDeg(); e++)
	{
		ui u = NV.GetOutNId(e);
		u = u - 1;
		if (K == cores[u]) {
			if (!visited_[u]) {
				deg_[u] = mcd_[u];
				visited_[u] = true;
				to_be_clear.push_back(u);
			}
			--deg_[u];
			if (deg_[u] < K && !evicted_[u]) {
				PropagateDismissal(K, u, to_be_clear, changed);
			}
		}
	}
}

void core::core_spanning(bool print)
{
	Timer timer;
	UnionFind *uf = new UnionFind(n);

	uf->init(n);
	

	bool *vis = new bool[n + 1];
	memset(vis, false, sizeof(bool)*(n+1));
	ui cnt_all_trees = 0;

	edg.set_m(G->GetEdges());
	edg.set_edges(G->GetEdges());

	ui pre_core = -1;

	ui id = 0;

	int *k_num = new int[n];
	memset(k_num, 0, sizeof(int)*(n));
	k_num[0] = -1;

	ui cnttt = 0;

	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		//cout << "peel[" << i-1 << "]:" << u << endl;
		TUNGraph::TNodeI NU = G->GetNI(u + 1);
		for (int e = 0; e < NU.GetOutDeg(); e++)
		{
			ui v = NU.GetOutNId(e);
			v--;
			cnttt++;
			if (vis[v])
			{
				edg.set_medges({ { u,v }, id });
				edg.set_medges({ { v,u }, id });
				edg.set_rmedges({ id,{u,v} });
				edg.set_corenum(cores[u]);

				if (cores[u] != pre_core)
				{
					if (pre_core == -1)
					{
						//edg.k_max = cores[u];
						edg.set_kmax(cores[u]);
						
						k_num[cores[u]] = 0;
						pre_core = cores[u];
					}
					else
					{
						while (pre_core != cores[u])
						{
							pre_core--;
							k_num[pre_core] = id;
						}
						pre_core = cores[u];
					}
				}
				edg.set_rid(id);
				edg.set_id(id);
				id++;

				ui ru = uf->UF_find(u);
				ui rv = uf->UF_find(v);
				if (ru != rv)
				{
					spt.pb(make_pair(make_pair(u, v), cores[u]));
					sptt.push(make_pair(make_pair(u, v), cores[u]));
					spt_num++;
					cnt_all_trees += cores[u];
					uf->UF_union(ru, rv);
					TG->AddEdge(u + 1, v + 1);
				}
			}
		}
		vis[u] = 1;
	}

	for (ui i = 0; i <= edg.get_kmax(); i++)
	{
		if (k_num[i] == 0 && i != edg.get_kmax())
		{
			if (i == 1)
			{
				k_num[i] = m;
			}
			else
			{
				k_num[i] = k_num[i - 1];
			}
		}
		
		edg.set_knum(k_num[i]);
	}

	cnt_tree = cnt_all_trees;
	
	if (print)
	{
		FILE *f = Utility::open_file((string("org_core_spanning_tree.txt")).c_str(), "w");
		fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
		fprintf(f, "weight of the spanning tree is: %u\n", cnt_all_trees);
		fprintf(f, "vertex vertex weight\n");
		map<PII, ui> t;
		for (ui i = 0; i < spt_num; i++)
		{
			auto tem = sptt.front();
			sptt.pop();
			sptt.push(tem);
			t.insert(tem);
			t.insert({ {tem.first.second,tem.first.first},tem.second });
			fprintf(f, "%u\t%u\t%u\n", tem.first.first, tem.first.second,tem.second);
		}

		for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++)
		{
			ui u = EI.GetSrcNId()-1, v = EI.GetDstNId()-1;
			if (t.find({ u,v })!=t.end())continue;
			else
			{
				fprintf(f, "%u\t%u\t%u\n", u, v, 0);
			}
		}

		for (ui i = 0; i < n; i++)
		{
			fprintf(f, "%u\t%u\n", i, cores[i]);
		}

		//for (ui i = spt.size(); i > 0; i--) fprintf(f, "%u\t%u\t%u\n", spt[i - 1].first.first, spt[i - 1].first.second, spt[i - 1].second);
		fclose(f);
	}
	delete[] vis;
	delete[] k_num;

	uf_tree = new UnionFind(*uf);
	delete uf;
}
void core::core_spanning_d(bool flag, ui u, ui v)
{
	if (flag == 0)
	{
		G->DelEdge(u + 1, v + 1);
		degree[u]--;
		degree[v]--;
		m--;
	}
	else
	{
		G->AddEdge(u + 1, v + 1);
		degree[u]++;
		degree[v]++;
		m++;
	}

	core_decompstion(true);
	UnionFind *uf = new UnionFind(n);

	uf->init(n);

	bool *vis = new bool[n + 1];
	memset(vis, false, sizeof(bool)*(n + 1));
	ui cnt_all_trees = 0;
	while (spt.size())spt.pop_back();
	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		TUNGraph::TNodeI NU = G->GetNI(u + 1);
		for (int e = 0; e < NU.GetOutDeg(); e++)
		{
			ui v = NU.GetOutNId(e);
			v--;
			if (vis[v])
			{
				ui ru = uf->UF_find(u);
				ui rv = uf->UF_find(v);
				if (ru != rv)
				{
					spt.pb(make_pair(make_pair(u, v), cores[u]));
					cnt_all_trees += cores[u];
					uf->UF_union(ru, rv);
				}
			}
		}
		vis[u] = 1;
	}
	cnt_tree = cnt_all_trees;


	delete[] vis;
	delete uf;
}


void core::Insert_core_spanning_new(ui u, ui v)
{
	int K = min(cores[u], cores[v]);
	unordered_set<int> *V_ = Insert(u, v);
	//Insertion(u, v);
	edg.Insert(u, v, min(cores[u], cores[v]));
	
	vector<pair<ui,PII>>E_star;
	unordered_map<PII, ui, pairHash>E_find;

	for (unordered_set<int>::iterator it = V_->begin(); it != V_->end(); it++)
	{
		ui w = *it;
		TUNGraph::TNodeI NW = G->GetNI(w + 1);
		for (int e = 0; e < NW.GetOutDeg(); e++)
		{
			ui y = NW.GetOutNId(e) - 1;

			if (cores[y] <= K)
			{
				continue;
			}
			else
			{
				if (E_find.find({ w,y }) == E_find.end())
				{
					E_star.push_back({ (ui)K + 1,{w,y} });
					E_find.insert({ {w,y},(ui)K + 1 });
					E_find.insert({ {y,w},(ui)K + 1 });
					if ((w == u && y == v) || (w == v && y == u))
					{
						continue;
					}
					edg.Update(w, y, K, K + 1);

					
				}
			}
		
		}
	}
	if (E_find.find({ u,v }) == E_find.end())
	{
		E_star.push_back({ (ui)min(cores[u],cores[v]),{u,v} });
		E_find.insert({ {u,v},(ui)min(cores[u],cores[v]) });
		E_find.insert({ {u,v},(ui)min(cores[u],cores[v]) });
	}
	

	queue<pair<PII,ui>>spt_tem;

	UnionFind *uf = new UnionFind(n);

	uf->init(n);

	ui cnt_all_trees = 0;

	ui old_num = sptt.size() + 1;
	ui i = 0;
	spt_num = 0;
	while ((i < E_star.size()) && sptt.size())
	{
		if (spt_num == old_num+1)break;
		auto spt_ = sptt.front();

		ui kwy = E_star[i].first;
		

		if (spt_.second >= kwy)
		{
			
			sptt.pop();
			PII wy = spt_.first;

			ui w = wy.first, y = wy.second;
			if (E_find.find({ w,y }) != E_find.end())continue;
			ui ru = uf->UF_find(w);
			ui rv = uf->UF_find(y);
			if (ru != rv)
			{
				spt_tem.push(make_pair(wy, spt_.second));
				spt_num++;
				cnt_all_trees += spt_.second;
				uf->UF_union(ru, rv);
			}
			else
			{
				if (TG->IsEdge(w + 1, y + 1))TG->DelEdge(w + 1, y + 1);
			}
		}
		else
		{
			PII wy = E_star[i].second;
			i++;
			ui w = wy.first, y = wy.second;

			ui ru = uf->UF_find(w);
			ui rv = uf->UF_find(y);
			if (ru != rv)
			{
				spt_num++;
				spt_tem.push(make_pair(wy, kwy));
				cnt_all_trees += kwy;
				uf->UF_union(ru, rv);
				if (!TG->IsEdge(w + 1, y + 1))TG->AddEdge(w + 1, y + 1);
			}
		}
	}

	while (sptt.size())
	{
		auto spt_ = sptt.front();
		sptt.pop();
		if (spt_num == old_num)continue;
		PII wy = spt_.first;

		ui w = wy.first, y = wy.second;
		if (E_find.find({ w,y }) != E_find.end())continue;
		ui ru = uf->UF_find(w);
		ui rv = uf->UF_find(y);
		if (ru != rv)
		{
			spt_tem.push(make_pair(wy, spt_.second));
			spt_num++;
			cnt_all_trees += spt_.second;
			uf->UF_union(ru, rv);
		}
		else
		{
			if (TG->IsEdge(w + 1, y + 1))TG->DelEdge(w + 1, y + 1);
		}
	}

	while (i < E_star.size())
	{
		if (spt_num == old_num)break;
		ui kwy = E_star[i].first;
		PII wy = E_star[i].second;
		i++;
		ui w = wy.first, y = wy.second;

		ui ru = uf->UF_find(w);
		ui rv = uf->UF_find(y);
		if (ru != rv)
		{
			spt_num++;
			spt_tem.push(make_pair(wy, kwy));
			cnt_all_trees += kwy;
			uf->UF_union(ru, rv);
			if(!TG->IsEdge(w+1,y+1))TG->AddEdge(w + 1, y + 1);
		}
	}


	while (spt_tem.size())
	{
		auto tt = spt_tem.front();
		spt_tem.pop();
		sptt.push(tt);
	}
	cnt_tree = cnt_all_trees;

	delete uf_tree;
	uf_tree = new UnionFind(*uf);

	delete uf;
}

void core::Delete_core_spanning_new(ui u, ui v)
{
	int K = min(cores[u], cores[v]);
	ui begin_id = edg.get_id(edg.get_medges({ u,v }));
	ui end_id = edg.get_id(edg.get_medges({ u,v }));

	unordered_set<int> *V_ = Remove(u, v);
	int flag_tree = 0;

	//Deletion(u, v);
	edg.Delete(u, v, K);

	
	unordered_map<PII, ui,pairHash>E_find;
	unordered_map<ui, ui>V_find;

	for (unordered_set<int>::iterator it = V_->begin(); it != V_->end(); it++)
	{
		V_find.insert({ *it,1 });
	}


	for (unordered_set<int>::iterator it = V_->begin(); it != V_->end(); it++)
	{
		ui w = *it;

		TUNGraph::TNodeI NW = G->GetNI(w + 1);
		for (int e = 0; e < NW.GetOutDeg(); e++)
		{
			ui y = NW.GetOutNId(e) - 1;
			if (cores[y] <= K - 1&&V_find.find(y)==V_find.end())
			{
				continue;
			}
			else
			{
				if (E_find.find({ w,y }) == E_find.end())
				{
					ui ttid = edg.get_id(edg.get_medges({ w,y }));
					begin_id = min(begin_id, ttid);
					E_find.insert({ {w,y},(ui)K - 1 });
					E_find.insert({ {y,w},(ui)K - 1 });

					if ((w == u && y == v) || (w == v && y == u))
					{
						continue;
					}
					edg.Update(w, y, K, K - 1);
					
					ui tt2id = edg.get_id(edg.get_medges({ w,y }));
					end_id = max(end_id, tt2id);
				}

			}
		}
	}
	if (begin_id == end_id)return;
	if (!TG->IsEdge(u, v))end_id = m;
	
	UnionFind *uf = new UnionFind(n);

	uf->init(n);
	ui cnt_all_trees = 0;
	
	queue<pair<PII,ui>>spt_tem;
	spt_num = 0;

	queue<pair<PII,ui>>spt_tem2;

	ui spt_size = sptt.size();
	while (sptt.size())
	{
		auto spt_ = sptt.front();
		sptt.pop();
		if ((spt_.first.first==u&&spt_.first.second==v)|| (spt_.first.first == v && spt_.first.second == u))continue;
		if (E_find.find(spt_.first) == E_find.end())
		{
			PII wy = spt_.first;

			ui w = wy.first, y = wy.second;

			ui ru = uf->UF_find(w);
			ui rv = uf->UF_find(y);
			if (ru != rv)
			{
				spt_tem.push(make_pair(wy, spt_.second));
				spt_num++;
				cnt_all_trees += spt_.second;
				uf->UF_union(ru, rv);
			}
			else
			{
				TG->DelEdge(w + 1, y + 1);
			}
		}
	}
	
	for (ui i = begin_id; i < end_id; i++)
	{	
		ui kwy = edg.get_corenum(i);
		PII wy = edg.get_rmedges(edg.get_rid(i));
		ui w = wy.first, y = wy.second;
		ui ru = uf->UF_find(w);
		ui rv = uf->UF_find(y);
		if (ru != rv)
		{
			spt_num++;
			spt_tem2.push(make_pair(wy, kwy));
			cnt_all_trees += kwy;
			uf->UF_union(ru, rv);
			TG->AddEdge(w + 1, y + 1);
		}
		if (spt_num == spt_size)break;
	}


	
	while (spt_tem.size() && spt_tem2.size())
	{
		auto t1 = spt_tem.front();
		auto t2 = spt_tem2.front();
		if (t1.second >= t2.second)
		{
			spt_tem.pop();
			sptt.push(t1);
		}
		else
		{
			spt_tem2.pop();
			sptt.push(t2);
		}
	}
	
	while (spt_tem.size())
	{
		sptt.push(spt_tem.front());
		spt_tem.pop();
	}

	while (spt_tem2.size())
	{
		sptt.push(spt_tem2.front());
		spt_tem2.pop();
	}
	cnt_tree = cnt_all_trees;
	
	delete uf_tree;
	uf_tree = new UnionFind(*uf);
	delete uf;
}


void core::Dynamic_maintainence(bool flag, ui u, ui v)
{
	core_decompstion(true);
	core_spanning(true);

	if (flag == 0)
	{
		Delete_core_spanning_new(u, v);
	}
	else
	{
		Insert_core_spanning_new(u, v);
	}
}


void core::Dynamic_maintainences(string str, bool flag,bool flag2,bool read_flag)
{
	if (flag == 0)ComputeCore(true);
	else core_decompstion(true);
	core_spanning(true);
	printf("# Start reading graph, Require files  \"update.txt\"\n");

	ifstream ifs;
	ifs.open(str, ifstream::in);

	std::string s;
	Timer timer;
	while (getline(ifs, s))
	{
		std::stringstream input(s);
		int x = 0;
		int y = 0;
		input >> x, input >> y;
		if(read_flag)x++, y++;
		
		if (flag2 == 0)
		{
			if (G->IsEdge(x, y) == 0)continue;

			if (flag == 0)
			{
				Delete_core_spanning_new(x - 1, y - 1);
			}
			else
			{
				core_spanning_d(0, x - 1, y - 1);
			}
		}
		else
		{
			if (G->IsEdge(x, y)==1)continue;
			
			if (flag == 0)
			{
				Insert_core_spanning_new(x - 1, y - 1);
			}
			else
			{
				core_spanning_d(1, x - 1, y - 1);
			}
		}
		
	}

	cout << "here we begin to compute" << endl;

	ifs.close();

	if (flag == 0)
	{
		FILE *f = Utility::open_file((string("core_spanning_tree.txt")).c_str(), "w");
		fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
		fprintf(f, "weight of the spanning tree is: %u\n", cnt_tree);
		fprintf(f, "vertex vertex weight\n");
		for (ui i = 0; i < spt_num; i++)
		{
			auto tem = sptt.front();
			sptt.pop();
			sptt.push(tem);

			fprintf(f, "%u\t%u\t%u\n", tem.first.first, tem.first.second, tem.second);
		}

		for (ui i = 0; i < n; i++)
		{
			fprintf(f, "%u\t%u\n", i, cores[i]);
		}
		fclose(f);
	}
	else
	{
		FILE *f = Utility::open_file((string("D://学习资料//科研资料//datasets//static_core_spanning_tree.txt")).c_str(), "w");
		fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
		fprintf(f, "weight of the spanning tree is: %u\n", cnt_tree);
		fprintf(f, "vertex vertex weight\n");
		for (ui i = spt.size(); i > 0; i--) fprintf(f, "%u\t%u\t%u\n", spt[i - 1].first.first, spt[i - 1].first.second, spt[i - 1].second);
		fclose(f);
	}
}

void core::Dynamic_maintainences2(string str, bool flag, bool flag2, bool read_flag)
{
	if (flag == 0)ComputeCore(true);
	else core_decompstion(true);
	core_spanning(true);
	printf("# Start reading graph, Require files  \"update.txt\"\n");

	FILE *f = Utility::open_file((string("core_spanning_tree_stablity.txt")).c_str(), "w");

	ifstream ifs;
	ifs.open(str, ifstream::in);

	std::string s;
	
	int cnt_all = 0;
	while (cnt_all < 100)
	{
		cnt_all++;
		int single_cnt = 0;
		Timer timer;
		while (single_cnt < 500)
		{
			single_cnt++;
			
			getline(ifs, s);

			std::stringstream input(s);
			int x = 0;
			int y = 0;
			input >> x, input >> y;
			if (read_flag)x++, y++;

			if (flag2 == 0)
			{
				if (G->IsEdge(x, y) == 0)continue;

				if (flag == 0)
				{
					Delete_core_spanning_new(x - 1, y - 1);
				}
				else
				{
					core_spanning_d(0, x - 1, y - 1);
				}
			}
			else
			{
				if (G->IsEdge(x, y) == 1)continue;
				if (flag == 0)
				{
					Insert_core_spanning_new(x - 1, y - 1);
				}
				else
				{
					core_spanning_d(1, x - 1, y - 1);
				}
			}

		}
		fprintf(f, "%s\n",Utility::integer_to_string(timer.elapsed()).c_str());
		printf("%s\n", Utility::integer_to_string(timer.elapsed()).c_str());
		//fprintf(f, "\tTotal processing time of the %d epoch: %s\n\n", cnt_all, Utility::integer_to_string(timer.elapsed()).c_str());
	}
	ifs.close();
	fclose(f);
}


void core::get_kcores(int k)
{
	FILE *f = Utility::open_file((string("k_cores.txt")).c_str(), "w");

	PUNGraph TG;  // 无向图
	TG = TUNGraph::New();

	for (int i = 1; i <= n; i++)
	{
		TG->AddNode(i);
	}
	int cnt = 0;
	while (sptt.size())
	{
		auto x = sptt.front();
		//cout << x.second << endl;
		int weight = x.second;
		auto vu = x.first;
		auto v = vu.first, u = vu.second;
		sptt.pop();
		//cout << v << ' ' << u << endl;
		TG->AddEdge(v + 1, u + 1);
	}
	//Timer timer;
	bool *vis;
	vis = new bool[n];
	memset(vis, 0, sizeof(bool)* n);

	vector<vector<PII>>k_core;
	for (int i = 0; i < n; i++)k_core.push_back(vector<PII>());

	
	cout << "-----------------" << endl;
	queue<int>q;
	for (int i = 0; i < n; i++)
	{
		if (!vis[i]&&cores[i]>=k)
		{
			q.push(i);
			vis[i] = true;
			int weight = k;
			while (q.size())
			{
				int v = q.front();
				q.pop();
				//k_core[i].push_back({ v,weight });
				TUNGraph::TNodeI NV = TG->GetNI(v + 1);
				for (int e = 0; e < NV.GetOutDeg(); e++)
				{
					ui u = NV.GetOutNId(e);
					u = u - 1;
					//cout << u << ' ' << v << endl;
					cnt++;
					if (vis[u]||cores[u]<k)continue;
					q.push(u);
					vis[u] = true;
					//TG->DelEdge(u+1, v+1);
					
				}
			}
		}
	}
	
	//fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
	/*fprintf(f, "all_cnt_is:%d\n", cnt);
	for (int i = 0; i < n; i++)
	{
		if (k_core[i].size() != 0)
		{
			sort(k_core[i].begin(), k_core[i].end());
			fprintf(f, "a %d-core is:",k);
			for (int j = 0; j < k_core[i].size(); j++)
			{
				fprintf(f, "%d ", k_core[i][j]);
			}
			fprintf(f, "\n");
		}
	}*/
	
	//cout << "all_cnt_is:" << cnt << endl;
	//fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
	fclose(f);
	delete vis;
}

long long core::get_kcores_(int k)
{
	int cnt = 0;
	Timer timer;
	bool *vis;
	vis = new bool[n];
	memset(vis, 0, sizeof(bool)* n);

	vector<vector<PII>>k_core;
	for (int i = 0; i < n; i++)k_core.push_back(vector<PII>());
	
	queue<int>q;
	for (int i = 0; i < n; i++)
	{
		if (!vis[i] && cores[i] >= k)
		{
			q.push(i);
			vis[i] = true;
			int weight = k;
			while (q.size())
			{
				int v = q.front();
				q.pop();
				k_core[i].push_back({ v,weight });
				TUNGraph::TNodeI NV = TG->GetNI(v + 1);
				for (int e = 0; e < NV.GetOutDeg(); e++)
				{
					ui u = NV.GetOutNId(e);
					u = u - 1;
					cnt++;
					if (vis[u] || cores[u] < k)continue;
					q.push(u);
					vis[u] = true;

				}
			}
		}
	}
	int ttt = timer.elapsed();
	
	for (int i = 0; i < n; i++)
	{
		if (k_core[i].size() != 0)
		{
			sort(k_core[i].begin(), k_core[i].end());
			printf("a %d-core is:", k);
			for (int j = 0; j < k_core[i].size(); j++)
			{
				printf("%d ", k_core[i][j]);
			}
			printf("\n");
		}
	}
	printf("\n");
	
	
	delete vis;

	return ttt;
}


void core::bfs_kcores(int k)
{
	FILE *f = Utility::open_file((string("bfs_k_cores.txt")).c_str(), "w");
	//Timer timer;
	bool *vis;
	vis = new bool[n];
	memset(vis, 0, sizeof(bool)* n);

	vector<vector<PII>>k_core;
	for (int i = 0; i < n; i++)k_core.push_back(vector<PII>());

	int cnt = 0;
	
	queue<int>q;
	for (int i = 0; i < n; i++)
	{
		if (!vis[i] && cores[i] >= k)
		{
			q.push(i);
			vis[i] = true;
			int weight = k;
			while (q.size())
			{
				int v = q.front();
				q.pop();
				k_core[i].push_back({ v,weight });
				TUNGraph::TNodeI NV = G->GetNI(v + 1);
				for (int e = 0; e < NV.GetOutDeg(); e++)
				{
					ui u = NV.GetOutNId(e);
					u = u - 1;
					cnt++;
					if (vis[u] || cores[u] < k)continue;
					q.push(u);
					vis[u] = true;
				}
			}
		}
	}
	
	//fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
	fprintf(f, "all_cnt_is:%d\n", cnt);
	for (int i = 0; i < n; i++)
	{
		if (k_core[i].size() != 0)
		{
			sort(k_core[i].begin(), k_core[i].end());
			fprintf(f, "a %d-core is:", k);
			for (int j = 0; j < k_core[i].size(); j++)
			{
				fprintf(f, "%d ", k_core[i][j]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
	
	delete vis;
}

long long core::bfs_kcores_(int k)
{
	Timer timer;
	bool *vis;
	vis = new bool[n];
	memset(vis, 0, sizeof(bool)* n);

	vector<vector<PII>>k_core;
	for (int i = 0; i < n; i++)k_core.push_back(vector<PII>());

	int cnt = 0;

	queue<int>q;
	for (int i = 0; i < n; i++)
	{
		if (!vis[i] && cores[i] >= k)
		{
			q.push(i);
			vis[i] = true;
			int weight = k;
			while (q.size())
			{
				int v = q.front();
				q.pop();
				k_core[i].push_back({ v,weight });
				TUNGraph::TNodeI NV = G->GetNI(v + 1);
				for (int e = 0; e < NV.GetOutDeg(); e++)
				{
					ui u = NV.GetOutNId(e);
					u = u - 1;
					cnt++;
					if (vis[u] || cores[u] < k)continue;
					q.push(u);
					vis[u] = true;
				}
			}
		}
	}
	
	int ttt = timer.elapsed();

	for (int i = 0; i < n; i++)
	{
		if (k_core[i].size() != 0)
		{
			sort(k_core[i].begin(), k_core[i].end());
			printf("a %d-core is:", k);
			for (int j = 0; j < k_core[i].size(); j++)
			{
				printf("%d ", k_core[i][j]);
			}
			printf("\n");
		}
	}
	printf("\n");
	delete vis;
	return ttt;
}


void core::connected_kcores(int k)
{
	FILE *f = Utility::open_file((string("connected_k_cores.txt")).c_str(), "w");
	//Timer timer;*/
	bool *vis;
	vis = new bool[n];
	memset(vis, 0, sizeof(bool)* n);

	vector<vector<PII>>k_core;
	for (int i = 0; i < n; i++)k_core.push_back(vector<PII>());
	int cnt = 0;
	queue<int>q;
	for (int i = 0;i < n; i++)
	{
		if (degree[i] < k)q.push(i);
	}

	while (q.size())
	{
		int v = q.front();
		q.pop();
		TUNGraph::TNodeI NV = G->GetNI(v + 1);
		for (int e = 0; e < NV.GetOutDeg(); e++)
		{
			cnt++;
			ui u = NV.GetOutNId(e);
			u = u - 1;
			degree[u] = degree[u] - 1;
			//cout << u << ' ' << degree[u] << endl;
			if (degree[u] == k-1)q.push(u);
			//if(G->IsEdge(u+1,v+1))G->DelEdge(u+1, v+1);
		}
		G->DelNode(v + 1);
	}


	for (int i = 0; i < n; i++)
	{
		if (!G->IsNode(i + 1))
		{
			vis[i] = true;
			continue;
		}
		if (vis[i]||degree[i]<k)continue;

		q.push(i);
		vis[i] = true;
		while (q.size())
		{
			int v = q.front();
			q.pop();
			k_core[i].push_back({ v,k });
			TUNGraph::TNodeI NV = G->GetNI(v + 1);
			for (int e = 0; e < NV.GetOutDeg(); e++)
			{
				ui u = NV.GetOutNId(e);
				u = u - 1;
				cnt++;
				if (vis[u]||degree[u]<k)continue;
				q.push(u);
				vis[u] = true;
			}
		}
	}
	
	//fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
	fprintf(f,"all_cnt_is:%d\n", cnt);
	for (int i = 0; i < n; i++)
	{
		if (k_core[i].size() != 0)
		{
			sort(k_core[i].begin(), k_core[i].end());
			fprintf(f, "a %d-core is:", k);
			for (int j = 0; j < k_core[i].size(); j++)
			{
				fprintf(f, "%d ", k_core[i][j]);
			}
			fprintf(f, "\n");
		}
	}

	fclose(f);

	delete vis;
}

PUNGraph CopyGraph(const PUNGraph& G)
{
	PUNGraph NewG = TUNGraph::New();
	for (TUNGraph::TNodeI NI = G->BegNI(); NI != G->EndNI(); NI++)
	{
		NewG->AddNode(NI.GetId());
	}
	for (TUNGraph::TEdgeI EI = G->BegEI(); EI != G->EndEI(); EI++)
	{
		NewG->AddEdge(EI.GetSrcNId(), EI.GetDstNId());
	}
	return NewG;
}

long long core::connected_kcores_(int k)
{
	Timer timer;
	ui *tem_degree = new ui[n + 1];
	memcpy(tem_degree, degree, sizeof(ui)*(n + 1));
	bool *vis;
	vis = new bool[n];
	memset(vis, 0, sizeof(bool)* n);

	PUNGraph TG = TUNGraph::New();
	TG = CopyGraph(G);
	

	vector<vector<PII>>k_core;
	for (int i = 0; i < n; i++)k_core.push_back(vector<PII>());
	int cnt = 0;
	queue<int>q;
	for (int i = 0; i < n; i++)
	{
		if (tem_degree[i] < k)q.push(i);
	}

	while (q.size())
	{
		int v = q.front();
		q.pop();
		TUNGraph::TNodeI NV = TG->GetNI(v + 1);
		for (int e = 0; e < NV.GetOutDeg(); e++)
		{
			cnt++;
			ui u = NV.GetOutNId(e);
			u = u - 1;
			tem_degree[u] = tem_degree[u] - 1;
			if (tem_degree[u] == k - 1)q.push(u);
		}
		TG->DelNode(v + 1);
	}


	for (int i = 0; i < n; i++)
	{
		if (!TG->IsNode(i + 1))
		{
			vis[i] = true;
			continue;
		}
		if (vis[i] || tem_degree[i] < k)continue;

		q.push(i);
		vis[i] = true;
		while (q.size())
		{
			int v = q.front();
			q.pop();
			k_core[i].push_back({ v,k });
			TUNGraph::TNodeI NV = TG->GetNI(v + 1);
			for (int e = 0; e < NV.GetOutDeg(); e++)
			{
				ui u = NV.GetOutNId(e);
				u = u - 1;
				cnt++;
				if (vis[u] || tem_degree[u] < k)continue;
				q.push(u);
				vis[u] = true;
			}
		}
	}
	delete vis;
	delete tem_degree;
	return timer.elapsed();
}



void core::Dynamic_kcores(string str, int flag, bool flag2, bool read_flag, int k1, int k2)
{
	if (flag == 0||flag==1)ComputeCore(true);
	else core_decompstion(true);
	core_spanning(true);
	printf("# Start reading graph, Require files  \"update.txt\"\n");

	ifstream ifs;
	ifs.open(str, ifstream::in);

	std::string s;
	Timer timer;
	long long alltime = 0;
	int cntt = 0;
	while (getline(ifs, s))
	{
		std::stringstream input(s);
		int x = 0;
		int y = 0;
		input >> x, input >> y;
		if (read_flag)x++, y++;
		int k = edg.get_kmax();
		cntt++;
		if (flag2 == 0)
		{
			if (G->IsEdge(x, y) == 0)continue;

			if (flag == 0)
			{
				Timer timer2;
				Delete_core_spanning_new(x - 1, y - 1);
				alltime += timer2.elapsed();	
				for(int kk =k1; kk<=k2; kk++) alltime += get_kcores_(kk);
			}
			else if(flag == 1)
			{
				Timer timer2;
				Remove(x-1, y-1);
				alltime += timer2.elapsed();
				for (int kk = k1; kk <= k2; kk++) alltime += bfs_kcores_(kk);
			}
		}
		else
		{
			if (G->IsEdge(x, y) == 1)
			{
				G->DelEdge(x, y);
			}
			if (flag == 0)
			{
				Timer timer2;
				Insert_core_spanning_new(x - 1, y - 1);
				alltime += timer2.elapsed();
				for (int kk = k1; kk <= k2; kk++) alltime += get_kcores_(kk);
			}
			else if (flag == 1)
			{
				Timer timer2;
				Insert(x-1, y-1);
				alltime += timer2.elapsed();
				for (int kk = k1; kk <= k2; kk++) alltime += bfs_kcores_(kk);
			}
		}

	}


	ifs.close();

	if (flag == 0)
	{
		FILE *f = Utility::open_file((string("core_spanning_tree.txt")).c_str(), "w");
		fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(alltime).c_str());

		fclose(f);
	}
	else if (flag == 1)
	{
		FILE *f = Utility::open_file((string("bfs_core_spanning_tree.txt")).c_str(), "w");
		fprintf(f, "\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(alltime).c_str());

		fclose(f);
	}
}




