#include "edges.h"


edges::edges()
{

}

edges::~edges()
{

}

edges::edges(ui _m, ui _all_edges)
{
	m = _m;
	all_edges = all_edges;
	k_max = 0;
}


void edges::Insert(ui u, ui v, ui core_)
{
	ui edges_id;
	if (m_edges.find({ u,v }) == m_edges.end())
	{
		edges_id = all_edges;
		m_edges.insert({ {u,v},all_edges });
		m_edges.insert({ { v,u }, all_edges});
		rm_edges.insert({ all_edges,{u,v} });
	}
	else
	{
		edges_id = m_edges[{u, v}];
	}
	all_edges++;
	m++;

	core_num.push_back(core_);
	id.push_back(edges_id);
	rid.push_back(edges_id);

	if (core_ > k_max)
	{
		ui begin_k = k_max+1;
		k_max = core_;

		while (begin_k <= core_)
		{
			if(begin_k!=core_)k_num.push_back(1);
			else k_num.push_back(0);
			begin_k++;
		}
		//k_num[k_max] = 0;
	}

	for (ui i = 1; i < core_; i++)
	{
		int pos = k_num[i];
		if (pos == id[edges_id])continue;

		swap(core_num[pos], core_num[id[edges_id]]);
		ui tt = id[edges_id];
		swap(id[rid[pos]], id[edges_id]);
		swap(rid[pos], rid[tt]);
		edges_id = rid[pos];

		if (core_ > i)
		{
			k_num[i] = pos + 1;
		}
	}

	ui i = core_ + 1;
	while (i <= k_max)
	{
		if (k_num[i] == 0)break;
		if (core_num[k_num[i] - 1] <= i)
		{
			k_num[i] = k_num[i] - 1;
		}
		else break;
	}
	k_max = core_num[0];
}

void edges::Delete(ui u, ui v, ui core_)
{
	ui edges_id;
	edges_id = m_edges[{u, v}];

	ui begin_id=id[edges_id];


	for (ui i = core_; i >= 1; i--)
	{
		if (k_num[i] == m)continue;
		int pos = k_num[i - 1] - 1;
		if (i - 1 == 0)pos = all_edges-1;

		if (id[edges_id] == k_num[i - 1] - 1)
		{
			//k_num[i] = k_num[i-1];
			if (i - 1 != 0)
			{
				k_num[i - 1] --;
			}
			continue;
		}

		swap(core_num[pos], core_num[id[edges_id]]);

		ui tt = id[edges_id];
		swap(id[rid[pos]], id[edges_id]);
		swap(rid[pos], rid[tt]);
		edges_id = rid[pos];

		if (i - 1 != 0)
		{
			k_num[i - 1] --;
		}
	}

	



	ui end_id = id[all_edges - 1];
	PII end_edges = rm_edges[all_edges-1];
	
	rid[id[all_edges-1]] = edges_id;
	id[edges_id] = end_id;

	m_edges[end_edges] = edges_id;
	m_edges[{end_edges.second,end_edges.first}] = edges_id;
	m_edges.erase({ u,v });
	m_edges.erase({ v,u });

	rm_edges[edges_id] = end_edges;
	rm_edges.erase(all_edges - 1);

	id.pop_back();
	core_num.pop_back();
	rid.pop_back();

	m--;
	all_edges--;
	k_max = core_num[0];

}

void edges::Update(ui u, ui v, ui core_old, ui core_new)
{
	int flag = 0;  // 1 为增加1 -1为减少1
	if (core_new == core_old + 1)
	{
		flag = 1;
	}
	else flag = -1;
	
	ui edges_id;
	edges_id = m_edges[{u, v}];

	/*
	cout << "cores:";
	for (ui i = 0; i < m; i++)
	{
		cout << core_num[i] << ' ';
	}
	cout << endl;

	cout << "rid:";
	for (ui i = 0; i < m; i++)
	{
		cout << rid[i] << ' ';
	}
	cout << endl;

	cout << "id:";
	for (ui i = 0; i < m; i++)
	{
		cout << id[i] << ' ';
	}
	cout << endl;

	cout << "k_num:";
	for (ui i = 0; i <= k_max; i++)
	{
		cout << k_num[i] << ' ';
	}
	cout << endl;
	*/
	if (flag == 1)
	{
		ui pos0 = k_num[core_old];
		if (pos0 != id[edges_id])
		{
			swap(core_num[pos0], core_num[id[edges_id]]);
			ui tt = id[edges_id];
			swap(id[rid[pos0]], id[edges_id]);

			
			swap(rid[pos0], rid[tt]);
			edges_id = rid[pos0];
		}


		/* ui pos = k_num[core_new];
		if (core_num[pos]!=core_new)
		{
			k_num[core_new] = edges_id;	
		}*/
		core_num[id[edges_id]]++;
		if (id[edges_id] + 1 >= all_edges)
		{
			k_num[core_old] = -1;
		}
		else
		{
			k_num[core_old]++;
		}
	}
	else
	{
		ui pos0 = k_num[core_old-1]-1;
	
		if (pos0 != id[edges_id])
		{
			swap(core_num[pos0], core_num[id[edges_id]]);

			ui tt = id[edges_id];
			swap(id[rid[pos0]], id[edges_id]);

			swap(rid[pos0], rid[tt]);
			edges_id = rid[pos0];
		}

		ui pos = k_num[core_new];
		k_num[core_new] = id[edges_id];
		/*
		ui tt = 1;
		while (k_num[core_new-tt]==id[edges_id]&&core_new-tt>0)
		{
			if (core_num[k_num[core_new - tt]] != core_new - tt)
			{
				k_num[core_new - tt] = edges_id + 1;
				tt++;
			}
			else
			{
				break;
			}
		}*/
		core_num[id[edges_id]]--;
	
	}
	k_max = core_num[0];
	/*
	cout << "cores:";
	for (ui i = 0; i < m; i++)
	{
		cout << core_num[i] << ' ';
	}
	cout << endl;

	cout << "rid:";
	for (ui i = 0; i < m; i++)
	{
		cout << rid[i] << ' ';
	}
	cout << endl;

	cout << "id:";
	for (ui i = 0; i < m; i++)
	{
		cout << id[i] << ' ';
	}
	cout << endl;

	cout << "k_num:";
	for (ui i = 0; i <= k_max; i++)
	{
		cout << k_num[i] << ' ';
	}
	cout << endl;*/
}



