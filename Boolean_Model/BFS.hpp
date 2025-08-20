#include <queue>

void Print_AL(std::vector<std::vector<int>> G)
{
	for(int i = 0; i < G.size(); i++)
	{
		std::cout << i << "-->";
		for(int j = 0; j < G[i].size(); j++)
		{
			std::cout << G[i][j] << " ";
		}
		std::cout << "\n";
	}
}

std::vector<std::vector<int>> Topology_To_AL(RegNetwork R)
{
	std::vector<std::vector<int>> G;
	int source, target;
	for(int i = 0; i < R.numnodes; i++)
	{
		G.push_back(std::vector<int>());
	}
	for(int i = 0; i < R.numedges; i++)
	{
		source = R.interactions[i][0];
		target = R.interactions[i][1];
		G[source].push_back(target);
	}

	return(G);
}

std::vector<int> BFS(std::vector<std::vector<int>> G, int source)
{
	int u, v;
	std::vector<int> color, d, pi;
	for(int i = 0; i < G.size(); i++)
	{
		color.push_back(0);
		d.push_back(INFINITY);
		pi.push_back(-1);
	}
	color[source] = 1;
	d[source] = 0;
	pi[source] = -1;

	std::queue<int> Q;
	Q.push(source);
	while(Q.size() > 0)
	{
		u = Q.front();
		Q.pop();
		for(int i = 0; i < G[u].size(); i++)
		{
			v = G[u][i];
			if(color[v] == 0)
			{
				color[v] = 1;
				d[v] = d[u] + 1;
				pi[v] = u;
				Q.push(v);
			}
		}
		color[u] = 2;
	}

	return(color);
}

double GRC(RegNetwork R)
{
	std::vector<std::vector<int>> G = Topology_To_AL(R);
	std::vector<int> color;
	std::vector<double> LRC;
	for(int i = 0; i < R.numnodes; i++)
	{
		color = BFS(G, i);
		LRC.push_back(0.0);
		for(int j = 0; j < color.size(); j++)
		{
			if(color[j] > 0 && j != i)
			{
				LRC[i] += 1.0;
			}
		}
		LRC[i] = LRC[i] / ((double) (R.numnodes - 1));
	}
	double max_LRC = -1.0;
	for(int i = 0; i < LRC.size(); i++)
	{
		if(LRC[i] > max_LRC)
		{
			max_LRC = LRC[i];
		}
	}
	double GRC = 0.0;
	for(int i = 0; i < LRC.size(); i++)
	{
		GRC += max_LRC - LRC[i];
	}
	GRC = GRC / ((double) R.numnodes - 1);

	return(GRC);
}