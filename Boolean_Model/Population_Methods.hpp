std::vector<RegNetwork> Generate_Same_Network_Population(RegNetwork R, int popusize)
{
	std::vector<RegNetwork> P;
	for(int i = 0; i < popusize; i++)
	{
		P.push_back(R);
	}

	return(P);
}

std::vector<RegNetwork> Generate_Random_Network_Population(RegNetwork R, int popusize)
{
	std::vector<RegNetwork> P;
	for(int i = 0; i < popusize; i++)
	{
		P.push_back(R);
		for(int j = 0; j < (int) ((std::log(1e7)*R.numedges) / 2.0); j++)
		{
			Mutate_Network(&P[i]);
		}
	}

	return(P);
}

void Simulate_No_Selection(std::vector<RegNetwork> P, int num_generations, int num_steps, std::ofstream *f, std::ofstream *g, std::ofstream *h, std::ofstream *meanfile, int world_rank)
{
	std::ofstream networkfile;

	std::vector<RegNetwork> New_P;
	std::vector<std::vector<int>> F, score, states_list;
	for(int i = 0; i < P.size(); i++)
	{
		F.push_back(std::vector<int>());
		score.push_back(std::vector<int>());
		for(int j = 0; j < num_steps; j++)
		{
			F[i].push_back(0);
			score[i].push_back(0);
		}
	}
	std::vector<int> state;
	double meanf = 0.0;

	for(int numgen = 0; numgen < num_generations; numgen++)
	{
		meanf = 0.0;
		for(int i = 0; i < P.size(); i++)
		{
			states_list.clear();
			for(int j = 0; j < num_steps; j++)
			{
				state = GetRandomState(P[i].numnodes);
				state = TrajectoryAsynUpdate(5000, state, P[i]);
				F[i][j] = CalculateFrustration(state, P[i]);
				score[i][j] = Get_Phenotypic_Score(state, P[i]);
				states_list.push_back(state);
				if(numgen % 50 == 0)
				{
					(*f) << F[i][j] << " ";
					(*g) << score[i][j] << " ";
				}
				meanf += F[i][j];
			}
			if(numgen % 50 == 0)
			{
				for(int j = 0; j < num_steps; j++)
				{
					for(int k = j + 1; k < num_steps; k++)
					{
						(*h) << StateOverlap(states_list[j], states_list[k], P[i]) << " ";
					}
				}
			}
		}
		(*meanfile) << meanf / ((double) (P.size()*num_steps)) << "\n";
		if(numgen % 50 == 0)
		{
			(*f) << "\n";
			(*g) << "\n";
			(*h) << "\n";

			for(int netindex = 0; netindex < P.size(); netindex++)
			{
				networkfile.open("outputfiles/network_files" + std::to_string(world_rank) + "/network_" + std::to_string(numgen) + "_" + std::to_string(netindex) + ".topo");
				Print_Network(P[netindex], &networkfile);
				networkfile.close();
			}
		}
		New_P.clear();
		for(int i = 0; i < P.size(); i++)
		{
			New_P.push_back(P[i]);
		}
		for(int i = 0; i < P.size(); i++)
		{
			P[i] = New_P[(int) (distribution(generator)*New_P.size())];
			if(distribution(generator) < 0.05)
			{
				Mutate_Network(&P[i]);
			}
		}
	}

	return;
}

std::vector<int> Get_Networks_Order(std::vector<std::vector<int>> F)
{
	std::vector<int> index;
	for(int i = 0; i < F.size(); i++)
	{
		index.push_back(i);
	}
	std::vector<double> quantity;
	for(int i = 0; i < F.size(); i++)
	{
		quantity.push_back((double) F[i][std::min_element(F[i].begin(), F[i].end()) - F[i].begin()]);
	}
	std::sort(index.begin(), index.end(), [&quantity] (size_t i1, size_t i2) {return quantity[i1] < quantity[i2];});

	return(index);
}

void Simulate_With_Selection(std::vector<RegNetwork> P, double selection_frac, int num_generations, int num_steps, std::ofstream *f, std::ofstream *g, std::ofstream *h, std::ofstream *meanfile, int world_rank)
{
	std::ofstream networkfile;

	std::vector<RegNetwork> New_P;
	std::vector<std::vector<int>> F, score, states_list;
	for(int i = 0; i < P.size(); i++)
	{
		F.push_back(std::vector<int>());
		score.push_back(std::vector<int>());
		for(int j = 0; j < num_steps; j++)
		{
			F[i].push_back(0);
			score[i].push_back(0);
		}
	}
	std::vector<int> state, index;
	double meanf = 0.0;

	for(int numgen = 0; numgen < num_generations; numgen++)
	{
		for(int i = 0; i < P.size(); i++)
		{
			states_list.clear();
			for(int j = 0; j < num_steps; j++)
			{
				state = GetRandomState(P[i].numnodes);
				state = TrajectoryAsynUpdate(5000, state, P[i]);
				F[i][j] = CalculateFrustration(state, P[i]);
				score[i][j] = Get_Phenotypic_Score(state, P[i]);
				states_list.push_back(state);
				if(numgen % 50 == 0)
				{
					(*f) << F[i][j] << " ";
					(*g) << score[i][j] << " ";
				}
				meanf += F[i][j];
			}
			if(numgen % 50 == 0)
			{
				for(int j = 0; j < num_steps; j++)
				{
					for(int k = j + 1; k < num_steps; k++)
					{
						(*h) << StateOverlap(states_list[j], states_list[k], P[i]) << " ";
					}
				}
			}
		}
		(*meanfile) << meanf / ((double) P.size()*num_steps) << "\n";
		if(numgen % 50 == 0)
		{
			(*f) << "\n";
			(*g) << "\n";
			(*h) << "\n";

			for(int netindex = 0; netindex < P.size(); netindex++)
			{
				networkfile.open("outputfiles/network_files" + std::to_string(world_rank) + "/network_" + std::to_string(numgen) + "_" + std::to_string(netindex) + ".topo");
				Print_Network(P[netindex], &networkfile);
				networkfile.close();
			}
		}
		New_P.clear();
		index = Get_Networks_Order(F);
		for(int i = 0; i < (int) (selection_frac*P.size()); i++)
		{
			New_P.push_back(P[index[i]]);
		}
		for(int i = 0; i < P.size(); i++)
		{
			P[i] = New_P[(int) (distribution(generator)*New_P.size())];
			if(distribution(generator) < 0.05)
			{
				Mutate_Network(&P[i]);
			}
		}
	}

	return;
}
