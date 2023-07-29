//#include "MSTruss.h"
//
//using namespace std;
//
//
//
//MSTruss::MSTruss(std::string _dir)
//{
//	dir = string(_dir);
//	n = 0;
//	m = 0;
//	G = TUNGraph::New();
//	H = TUNGraph::New();
//}
//
//void MSTruss::read_graph()
//{
//	printf("# Start reading graph, Reqintre files  \"corresponding_graph.txt\"\n");
//
//	ifstream ifs(dir.c_str(), std::ifstream::in);
//	
//	std::string s;
//	getline(ifs, s);
//	std::stringstream k(s);
//
//	//输入点数和边数,第一行放点和边的数量
//	k >> n, k >> m;
//
//	for (int i = 1; i <= n; i++)
//	{
//		G->AddNode(i);
//	}
//	int edgeid = 0;
//
//	while (getline(ifs, s))
//	{
//		std::stringstream input(s);
//		int x = 0;
//		int y = 0;
//		input >> x, input >> y;
//		G->AddEdge(x, y);
//		eid.insert({ make_pair(x,y), ++edgeid });
//		eid.insert({ make_pair(y,x),edgeid });
//		rid.insert({ edgeid,make_pair(x,y) });
//	}
//
//	ifs.close();
//}
//
//
//void MSTruss::PeelTruss()
//{
//	truss_ness = new int[m + 1];
//	int * supp = new int[m + 1];
//	memset(supp, 0, sizeof(int)*(m + 1));
//
//
//	/*计算出每个边的supp(u,v)*/
//	int max_truss = 0;
//	for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++)
//	{
//		int u = EI.GetSrcNId(), v = EI.GetDstNId();
//		TUNGraph::TNodeI NU = G->GetNI(u), NV = G->GetNI(v);
//		int edgeNum = eid[make_pair(u,v)]; // 这个地方等着看一下是不是从0开始编号
//
//		for (int e = 0; e < NU.GetOutDeg(); e++)
//		{
//			int w = NU.GetOutNId(e);
//			if (G->IsEdge(w, v))supp[edgeNum]++;
//		}
//		max_truss = max(max_truss, supp[edgeNum]);
//	}
//	
//	//int max_truss = 0;
//	vector<int> *bin = new vector<int>[max_truss + 1];
//	for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++)
//	{
//		int u = EI.GetSrcNId(), v = EI.GetDstNId();
//		TUNGraph::TNodeI NU = G->GetNI(u), NV = G->GetNI(v);
//		int edgeNum = eid[make_pair(u, v)];
//		bin[supp[edgeNum]].push_back(edgeNum);
//	}
//
//	vector<bool>del(m + 1, false);
//	
//	for (int i = 0; i <= max_truss; i++)
//	{
//		for (int j = 0; j < bin[i].size(); j++)
//		{
//			int id = bin[i][j];
//			if (del[id] == true)continue;
//			truss_ness[id] = i + 2;
//			del[id] = true;
//			
//			PII xy = rid[id];
//			int x = xy.first, y = xy.second;
//			//int y = EI.GetDstNId(), x = EI.GetSrcNId();
//			TUNGraph::TNodeI NU = G->GetNI(x), NV = G->GetNI(y);
//
//			if (NU.GetDeg() > NV.GetDeg())
//			{
//				TUNGraph::TNodeI tem = NU;
//				NU = NV;
//				NV = tem;
//
//				int t = x;
//				x = y;
//				y = t;
//			}
//
//			for (int e = 0; e < NU.GetOutDeg(); e++)
//			{
//				int w = NU.GetOutNId(e);
//				TUNGraph::TEdgeI wx = G->GetEI(x, w);
//
//				int wx_id = eid[make_pair(x, w)];    
//
//				if (del[wx_id] == true)continue;
//				if (G->IsEdge(w, y))
//				{
//					TUNGraph::TEdgeI wy = G->GetEI(y, w);
//					int wy_id = eid[make_pair(y, w)];
//					if (del[wy_id] == true)continue;
//
//					if (supp[wy_id] > i)
//					{
//						supp[wy_id]--;
//						bin[supp[wy_id]].push_back(wy_id);
//					}
//					
//					if (supp[wx_id] > i)
//					{
//						supp[wx_id]--;
//						bin[supp[wx_id]].push_back(wx_id);
//					}
//				}
//			}
//
//		}
//	}
//
//	for (int i = 1; i <= m; i++)
//	{
//		cout << "i:"<<truss_ness[i] << endl;
//	}
//
//	delete[] bin;
//}
//
//
//void MSTruss::addEdges(std::string _dir)
//{
//	printf("# Start reading graph, Reqintre files  \"corresponding_graph.txt\"\n");
//	add_edges = 0;
//	ifstream ifs(_dir.c_str(), std::ifstream::in);
//
//	std::string s;
//	getline(ifs, s);
//	std::stringstream k(s);
//
//	//输入点数和边数,第一行放点和边的数量
//	int nodes, edges;
//	k >> nodes, k >> edges;
//	add_edges = edges;
//
//	for (int i = 1; i <= nodes; i++)
//	{
//		H->AddNode(i);
//	}
//	//int edgeid = 0;
//
//	while (getline(ifs, s))
//	{
//		std::stringstream input(s);
//		int x = 0;
//		int y = 0;
//		input >> x, input >> y;
//		H->AddEdge(x, y);
//	}
//
//	ifs.close();
//}
//
//void MSTruss::Incremental_truss_maintenance(std::string _dir)
//{
//	addEdges(_dir);
//
//	while (add_edges!=0)
//	{
//		map<PII,int>E_ms;
//		map<PII, int>E_ms1;
//		map<PII, int>E_ms2;
//		queue<int>V_ms;
//
//		int newly_edges = 0;
//
//		while (!H.Empty())
//		{
//			for (TUNGraph::TNodeI NI = H->BegNI(); NI < H->EndNI(); NI++)
//			{
//				int flag = 0;
//				int v = NI.GetId();
//				queue<PII>E_deltav;
//
//				// find E△v Hard to do 
//				for (int u = 0; u < NI.GetOutDeg()-1; u++)
//				{
//					for (int w = u+1; w < NI.GetOutDeg(); w++)
//					{
//						int intd = NI.GetOutNId(u);
//						int wid = NI.GetOutNId(w);
//						if (H->IsEdge(u, w))
//						{
//							if (E_ms.find(make_pair(intd, wid)) != E_ms.end())
//							{
//								flag = 1;
//								break;
//							}
//							else
//							{
//								E_deltav.push(make_pair(intd, wid));
//								E_deltav.push(make_pair(wid, intd));
//							}
//
//							if (E_ms.find(make_pair(intd, v)) != E_ms.end())
//							{
//								flag = 1;
//								break;
//							}
//							else
//							{
//								E_deltav.push(make_pair(intd, v));
//								E_deltav.push(make_pair(v, intd));
//							}
//
//							if (E_ms.find(make_pair(wid, v)) != E_ms.end())
//							{
//								flag = 1;
//								break;
//							}
//							else
//							{
//								E_deltav.push(make_pair(wid, v));
//								E_deltav.push(make_pair(v, wid));
//							}
//						}
//						if (flag == 1)
//						{
//							break;
//						}
//					}
//					if (flag == 1)break;
//				}
//				if (flag == 0)
//				{
//					V_ms.push(v);
//					while (!E_deltav.empty())
//					{
//						PII tem = E_deltav.front();
//						E_deltav.pop();
//						E_ms.insert({ tem, 1 });
//					}
//				}
//
//				
//				
//			}
//		}
//
//		// Ems1 <- Union v\in MS E_G(v)
//		while (V_ms.size())
//		{
//			int vid = V_ms.front();
//			V_ms.pop();
//			TUNGraph::TNodeI v = H->GetNI(vid);
//			for (int e = 0; e < v.GetOutDeg(); e++)
//			{
//				E_ms1.insert({ make_pair(vid, v.GetOutNId(e)) ,1});
//				newly_edges++;
//			}
//		}
//
//		// △E <- △E \ E_ms1 
//		for (auto it = E_ms1.begin(); it != E_ms1.end(); it++)
//		{
//			int vid = it->first.first;
//			int wid = it->first.second;
//			H->DelEdge(vid, wid);
//		}
//
//
//		//for e\ in  \delta E do 
//		for (TUNGraph::TEdgeI EI = H->BegEI(); EI < H->EndEI(); EI++)
//		{
//			int vid = EI.GetSrcNId(), wid = EI.GetDstNId();
//			TUNGraph::TNodeI v = H->GetNI(vid), w = H->GetNI(wid);
//			queue<PII>E_deltae;
//
//			// find E_delta e
//			int flag = 0;
//			for (int e = 0; e < v.GetOutDeg(); e++)
//			{
//				int intd = v.GetOutNId(e);
//				if (H->IsEdge(intd, wid))
//				{
//					if (E_ms.find(make_pair(intd, wid)) == E_ms.end())
//					{
//						E_deltae.push(make_pair(intd, wid));
//						E_deltae.push(make_pair(wid, intd));
//					}
//					else
//					{
//						flag = 1;
//						break;
//					}
//					
//					if (E_ms.find(make_pair(intd, vid)) == E_ms.end())
//					{
//						E_deltae.push(make_pair(intd, vid));
//						E_deltae.push(make_pair(vid, intd));
//					}
//					else
//					{
//						flag = 1;
//						break;
//					}
//				}
//			}
//			if (flag == 0)
//			{
//				E_ms2.insert({make_pair(vid,wid),1});
//				newly_edges++;
//				while (E_deltae.size())
//				{
//					PII tem = E_deltae.front();
//					E_deltae.pop();
//					E_ms.insert({tem,1 });
//				}
//			}
//			
//		}
//		
//		// 
//		while (V_ms.size())
//		{
//			int x = V_ms.front();
//			V_ms.pop();
//			H->DelNode(x);
//			G->AddNode(x);
//		}
//
//		truss_ness = (int *)realloc(truss_ness, sizeof(int)*(m + 1+newly_edges));
//
//
//		// 这个地方感觉应该先拷贝一份
//		int *truss_copy = new int[m + 1 + newly_edges];
//
//		memcpy(truss_copy, truss_ness, sizeof(int)*(m + 1 + newly_edges));
//
//		
//
//		for (auto it = E_ms1.begin(); it != E_ms1.end(); it++)
//		{
//			int vid = it->first.first;
//			int wid = it->first.second;
//			G->AddEdge(vid, wid);
//			m++;
//			eid.insert({ make_pair(vid,wid), m });
//			eid.insert({ make_pair(wid,vid),m });
//			rid.insert({ m,make_pair(vid,wid) });
//		}
//
//		for (auto it = E_ms2.begin(); it != E_ms2.end(); it++)
//		{
//			int vid = it->first.first;
//			int wid = it->first.second;
//
//			G->AddEdge(vid, wid);
//			m++;
//			eid.insert({ make_pair(vid,wid), m });
//			eid.insert({ make_pair(wid,vid),m });
//			rid.insert({ m,make_pair(vid,wid) });
//			//E_ms1.insert(*it);
//		}
//
//		//map<PII, int>m_p;
//		multimap<int, PII,greater<int>>m_p;
//		map<PII, int>vis;
//
//		// compute=pre\tao(e_0)
//		for (auto it = E_ms1.begin(); it != E_ms1.end(); it++)
//		{
//			int vid = it->first.first;
//			int wid = it->first.second;
//			compute_pretruss1(vid, wid, truss_copy);
//		}
//
//		for (auto it = E_ms2.begin(); it != E_ms2.end(); it++)
//		{
//			int vid = it->first.first;
//			int wid = it->first.second;
//			compute_pretruss2(vid, wid, truss_copy);
//
//			E_ms1.insert(*it);
//		}
//
//		for (auto it = E_ms1.begin(); it != E_ms1.end(); it++)
//		{
//			int vid = it->first.first;
//			int wid = it->first.second;
//			
//			TUNGraph::TNodeI v = G->GetNI(vid);
//
//			for (int e = 0; e < v.GetOutDeg(); e++)
//			{
//				int intd = v.GetOutNId(e);
//				if (G->IsEdge(intd, wid))
//				{
//					int vw_truss = truss_ness[eid[{intd, wid}]];
//					int uw_truss = truss_ness[eid[{vid, wid}]];
//
//					int uv_truss= truss_ness[eid[{vid, intd}]];
//					int k = min(vw_truss, uw_truss);
//
//					if (k <= uv_truss)
//					{
//						if (vw_truss == k)
//						{
//							m_p.insert({k,{vid,wid} });
//							vis.insert({ {vid,wid} ,k });
//							//vis.insert({ {wid,vid} ,k });
//						}
//						else if (uw_truss == k)
//						{
//							m_p.insert({k, {intd,wid} });
//							vis.insert({ {intd,wid},k });
//							//vis.insert({ {wid,intd},k });
//						}
//					}
//				}
//			}
//		}
//		IncrementalTraversal(m_p,vis);
//		memcpy(truss_ness, truss_copy, sizeof(int)*(m + 1 + newly_edges));
//
//	}
//	
//}
//
//// E_ms1 pretruss compute
//void MSTruss::compute_pretruss1(int vid, int wid, int truss_copy[])
//{
//	TUNGraph::TNodeI w = G->GetNI(wid);
//	TUNGraph::TNodeI v = G->GetNI(vid);
//
//	for (int e = 0; e < w.GetOutDeg(); e++)
//	{
//		int intd = w.GetOutNId(e);
//		
//		
//	}
//
//
//}
//
//// E_ms2 pretruss compute
//void MSTruss::compute_pretruss2(int intd, int vid, int truss_copy[])
//{
//	int k = 0;
//	TUNGraph::TNodeI u = G->GetNI(intd);
//	TUNGraph::TNodeI v = G->GetNI(vid);
//
//	map<int, int>cntk;
//	priority_queue<int,vector<int>, less<int>>pk;
//
//	for (int e = 0; e < u.GetOutDeg(); e++)
//	{
//		int wid = u.GetOutNId(e);
//		if (G->IsEdge(wid, vid))
//		{
//			int uwtruss = truss_ness[eid[{wid, vid}]];
//			int vwtruss = truss_ness[eid[{intd, wid}]];
//
//			int minw = min(uwtruss, vwtruss);
//
//			if (cntk.find(minw) == cntk.end())
//			{
//				cntk.insert({ minw,1 });
//				pk.push(minw);
//			}
//			else cntk[minw]++;
//		}
//	}
//	int num = 0;
//
//	while (pk.size())
//	{
//		int nowk = pk.top();
//		pk.pop();
//		num += nowk;
//		if (num >= nowk - 2)
//		{
//			truss_copy[eid[{intd, vid}]] = nowk;
//			break;
//		}
//
//	}
//	while (pk.size())pk.pop();
//	
//	
//}
//
//
//void MSTruss::IncrementalTraversal(std::multimap<int, std::pair<int, int>, std::greater<int>> m_p, std::map<PII, int>vis)
//{
//	for (auto it = m_p.begin(); it != m_p.end(); it++)
//	{
//		int k = it->first;
//		auto itrange = m_p.equal_range(k);
//
//		queue<PII>q;
//		map<PII, int>S;
//
//		for (auto itr = itrange.first; itr != itrange.second; itr++)
//		{
//			q.push(itr->second);
//		}
//		
//
//		while (q.size())
//		{
//			PII edge = q.front();
//			q.pop();
//			S.insert({ edge ,0});
//
//			int xid = edge.first, yid = edge.second;
//
//			S.insert({ {yid,xid} ,0});
//			
//			TUNGraph::TNodeI NI = G->GetNI(xid);
//			for (int e = 0; e < NI.GetOutDeg(); e++)
//			{
//				int zid = NI.GetOutNId(e);
//				if (G->IsEdge(yid, zid))
//				{
//					int xztruss = truss_ness[ eid[{xid, zid}]];
//					int yztruss = truss_ness[ eid[{yid, zid}]];
//
//					if (min(xztruss, yztruss) < k)
//					{
//						continue;
//					}
//
//					S[edge]++;
//					S[{yid, xid}]++;
//
//					if (xztruss == k && (vis.find({ xid,zid }) == vis.end()&& vis.find({ zid,xid })==vis.end()))
//					{
//						q.push({ xid,zid });
//						m_p.insert({ k, { xid,zid } });
//						vis.insert({ { xid,zid }, k });
//						//vis.insert({ { zid,xid }, k });
//					}
//
//					if (yztruss == k && (vis.find({ yid,zid }) == vis.end()&& vis.find({ zid,yid }) == vis.end()))
//					{
//						q.push({ yid,zid });
//						m_p.insert({ k, { yid,zid } });
//						vis.insert({ { yid,zid }, k });
//						//vis.insert({ { zid,yid }, k });
//					}
//				}
//			}
//
//
//
//		}
//		
//
//		queue<PII>q_s;
//		for (auto s_c = S.begin(); s_c != S.end(); s_c++)
//		{
//			int tem_x= s_c->first.first,tem_y=s_c->first.second;
//			int tem_k = s_c->second;
//			if (tem_k <= k - 2)
//			{
//				q_s.push({ tem_x,tem_y });
//			}
//
//		}
//
//		while (q_s.size())
//		{
//			PII tem= q_s.front();
//			q_s.pop();
//			int tem_x = tem.first, tem_y = tem.second;
//
//			if (vis.find({ tem_x, tem_y }) != vis.end())
//			{
//				vis[{tem_x, tem_y}] = 0;
//			}
//			else vis[{tem_y, tem_x}] = 0;
//
//			TUNGraph::TNodeI Nx = G->GetNI(tem_x);
//			for (int e = 0; e < Nx.GetOutDeg(); e++)
//			{
//				int zid = Nx.GetOutNId(e);
//				if (G->IsEdge(zid, tem_y))
//				{
//					int xz_truss = truss_ness[eid[{tem_x, zid}]];
//					int yz_truss = truss_ness[eid[{tem_y, zid}]];
//
//					if (min(xz_truss, yz_truss) < k)continue;
//
//					
//					
//					if (yz_truss == k)
//					{
//						int flag = 0;
//						if (vis.find({ tem_y,zid }) == vis.end())flag = 1;
//						else if (vis[{tem_y, zid}] == 0)flag = 1;
//
//						if (vis.find({ zid,tem_y }) == vis.end())flag = 1;
//						else if (vis[{zid, tem_y}] == 0)flag = 1;
//
//						if (flag == 1)continue;
//						else
//						{
//							S[{tem_y, zid}]--;
//							S[{zid, tem_y}]--;
//						}
//
//							
//					}
//
//					if (xz_truss == k)
//					{
//						int flag = 0;
//						if (vis.find({ tem_x,zid }) == vis.end())flag = 1;
//						else if (vis[{tem_x, zid}] == 0)flag = 1;
//
//						if (vis.find({ zid,tem_x }) == vis.end())flag = 1;
//						else if (vis[{zid, tem_x}] == 0)flag = 1;
//
//						if (flag == 1)continue;
//						else
//						{
//							S[{tem_x, zid}]--;
//							S[{zid, tem_x}]--;
//						}
//					}
//				}
//			}
//
//		}
//		for (auto itr = itrange.first; itr != itrange.second; itr++)
//		{
//			PII edge = itr->second;
//			if (vis[edge] != 0)
//			{
//				truss_ness[eid[edge]] = k + 1;
//			}
//		}
//
//
//		it = itrange.second;
//		it--;
//	}
//}