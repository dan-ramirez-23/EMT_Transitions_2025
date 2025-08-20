#include <random>

std::mt19937 generator;
std::uniform_real_distribution<> distribution{0.0, 1.0};

void InitializeRandomNumGen(long seed)
{
	generator = std::mt19937(seed);
}

long BinaryToDecimal(std::vector<int> state)
{
	long dec_state = 0;
	for(int i = 0; i < state.size(); i++)
	{
		dec_state += state[i]*((long) std::pow(2, state.size() - i - 1));
	}

	return(dec_state);
}

std::vector<int> GetRandomState(int numnodes)
{
	std::vector<int> state;
	for(int i = 0; i < numnodes; i++)
	{
		if(distribution(generator) < 0.5)
		{
			state.push_back(0);
		}
		else
		{
			state.push_back(1);
		}
	}

	return(state);
}

std::vector<int> DecimalToBinary(long DecimalState, int numnodes)
{
    std::vector<int> state;
    if(DecimalState == 0)
    {
    	state.push_back(0);
    }
    else
    {
        while(DecimalState != 0)
        {
            state.push_back(DecimalState % 2);
            DecimalState = DecimalState / 2;
        }
    }

    if(state.size() > numnodes)
    {
    	std::cout << "Mistakes have been made in converting to binary for this network. Basically, you suck." << "\n";
    	return(std::vector<int>());
    }
    while(state.size() < numnodes)
    {
        state.push_back(0);
    }
    std::reverse(state.begin(), state.end());

    return(state);
}

int CalculateFrustration(std::vector<int> state, RegNetwork R)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "Invalid network state. Stand there in your wrongness and be wrong and get used to it." << "\n";
		return(-1);
	}

	int F = 0;

	int source, target, interactiontype, fedge;
	int s_i, s_j, J_ij;
	for(int i = 0; i < R.numedges; i++)
	{
		source = R.interactions[i][0];
		target = R.interactions[i][1];
		interactiontype = R.interactions[i][2];
		s_i = state[source];
		s_j = state[target];
		J_ij = interactiontype;
		if(s_i == 0)
		{
			s_i = -1;
		}
		if(s_j == 0)
		{
			s_j = -1;
		}
		if(J_ij == 2)
		{
			J_ij = -1;
		}
		fedge = s_i*s_j*J_ij;
		if(fedge < 0)
		{
			F += 1;
		}
	}

	return(F);
}

std::vector<int> TrajectoryAsynUpdateC(int numsteps, std::vector<int> state, RegNetwork R)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "Invalid network state. Run fast. Run far." << "\n";
		return(std::vector<int>());
	}

	int errflag = 0;

	int nodetoupdate, input;

	std::vector<int> nextstate;

	int source, target, interactiontype;
	for(int step = 0; step < numsteps; step++)
	{
		nodetoupdate = (int) (distribution(generator)*R.numnodes);
		if (R.Clamp[nodetoupdate] == 0)
		  {
		    int input1=0;
		    if (R.target_edges[nodetoupdate][0] != -1)
		      {
			for (int i=0;i < R.target_edges[nodetoupdate].size();i++)
			  {
			    int edge=R.target_edges[nodetoupdate][i];
			    int J=3 - 2*R.interactions[edge][2];
			    source = R.interactions[edge][0];
			    int s=2*state[source]-1;
			    input1 += J*s;
			  }
		      }
		    if(input1 > 0)
		      {
			state[nodetoupdate] = 1;
		      }
		    else if(input1 < 0)
		      {
			state[nodetoupdate] = 0;
		      }
		    else
		      {
			state[nodetoupdate] = state[nodetoupdate];
		      }
		  }
	}

	return(state);
}

std::vector<int> TrajectoryAsynUpdate(int numsteps, std::vector<int> state, RegNetwork R)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "Invalid network state. Run fast. Run far." << "\n";
		return(std::vector<int>());
	}

	int errflag = 0;

	int nodetoupdate, input;

	std::vector<int> nextstate;
	for(int i = 0; i < state.size(); i++)
	{
		nextstate.push_back(state[i]);
	}
	int source, target, interactiontype;
	for(int step = 0; step < numsteps; step++)
	{
		nodetoupdate = (int) (distribution(generator)*R.numnodes);
		/*
		input = 0;
		for(int i = 0; i < R.numedges; i++)
		{
			source = R.interactions[i][0];
			target = R.interactions[i][1];
			if(target != nodetoupdate)
			{
				continue;
			}
			interactiontype = R.interactions[i][2];
			if(interactiontype == 1)
			{
				if(nextstate[source] == 0)
				{
					input += (1)*(-1);
				}
				else if(nextstate[source] == 1)
				{
					input += (1)*(1);
				}
				else
				{
					errflag = 1;
				}
			}
			else if(interactiontype == 2)
			{
				if(nextstate[source] == 0)
				{
					input += (-1)*(-1);
				}
				else if(nextstate[source] == 1)
				{
					input += (-1)*(1);
				}
				else
				{
					errflag = 1;
				}
			}
			else
			{
				errflag = 1;
			}
		}
		*/
		int input1=0;
		if (R.target_edges[nodetoupdate][0] != -1)
		  {
		    for (int i=0;i < R.target_edges[nodetoupdate].size();i++)
		      {
			int edge=R.target_edges[nodetoupdate][i];
			int J=3 - 2*R.interactions[edge][2];
			source = R.interactions[edge][0];
			int s=2*state[source]-1;
			input1 += J*s;
		      }
		  }
		/*
		if (input1 != input)
		  {
		    std::cout << "BUG: " << nodetoupdate << " " << input << " " << input1 << "\n";
		    errflag=1;
		  }
		*/
		    
		
		if(errflag == 1)
		{
			std::cout << "Wicked weird stuff going on here. Take some time and think about what you have done." << "\n";
			return(std::vector<int>());
		}
		if(input1 > 0)
		{
			state[nodetoupdate] = 1;
		}
		else if(input1 < 0)
		{
			state[nodetoupdate] = 0;
		}
		else
		{
			state[nodetoupdate] = state[nodetoupdate];
		}
		//	nextstate[nodetoupdate] = state[nodetoupdate];
	}

	return(state);
}

double StateOverlap(std::vector<int> state1, std::vector<int> state2, RegNetwork R)
{
	if(state1.size() != R.numnodes || state2.size() != R.numnodes)
	{
		std::cout << "Invalid network state. This calculation will now come to a sudden arboreal stop." << "\n";
		return(-INFINITY);
	}

	double Q = 0.0;

	int s_a, s_b;
	for(int i = 0; i < state1.size(); i++)
	{
		s_a = state1[i];
		s_b = state2[i];
		if(s_a == 0)
		{
			s_a = -1;
		}
		if(s_b == 0)
		{
			s_b = -1;
		}
		Q += s_a*s_b;
	}
	Q = Q / (double) state1.size();

	return(Q);
}

std::vector<int> NextStateSynUpdate(std::vector<int> state, RegNetwork R)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "Invalid network state. Run fast. Run far." << "\n";
		return(std::vector<int>());
	}

	int errflag = 0;

	std::vector<int> nextstate, input;
	for(int i = 0; i < state.size(); i++)
	{
		nextstate.push_back(state[i]);
		input.push_back(0);
	}
	int source, target, interactiontype;
	for(int i = 0; i < R.numedges; i++)
	{
		source = R.interactions[i][0];
		target = R.interactions[i][1];
		interactiontype = R.interactions[i][2];
		if(interactiontype == 1)
		{
			if(state[source] == 0)
			{
				input[target] += (1)*(-1);
			}
			else if(state[source] == 1)
			{
				input[target] += (1)*(1);
			}
			else
			{
				errflag = 1;
			}
		}
		else if(interactiontype == 2)
		{
			if(state[source] == 0)
			{
				input[target] += (-1)*(-1);
			}
			else if(state[source] == 1)
			{
				input[target] += (-1)*(1);
			}
			else
			{
				errflag = 1;
			}
		}
		else
		{
			errflag = 1;
		}
	}
	for(int i = 0; i < R.numnodes; i++)
	{
		if(input[i] > 0)
		{
			nextstate[i] = 1;
		}
		else if(input[i] < 0)
		{
			nextstate[i] = 0;
		}
		else
		{
			nextstate[i] = state[i];
		}
	}

	if(errflag == 1)
	{
		std::cout << "Wicked weird stuff going on here. Take some time and think about what you have done." << "\n";
		return(std::vector<int>());
	}

	return(nextstate);
}

bool IsStableStateC(std::vector<int> state, RegNetwork R)
{
	std::vector<int> nextstate = NextStateSynUpdate(state, R);
	if(nextstate.size() == 0)
	{
		return(false);
	}
	for(int i = 0; i < state.size(); i++)
	{
	  if((nextstate[i] != state[i]) & (R.Clamp[i]==0))
	    {
	      return(false);
	    }
	}

	return(true);
}

bool IsStableState(std::vector<int> state, RegNetwork R)
{
	std::vector<int> nextstate = NextStateSynUpdate(state, R);
	if(nextstate.size() == 0)
	{
		return(false);
	}
	for(int i = 0; i < state.size(); i++)
	{
		if(nextstate[i] != state[i])
		{
			return(false);
		}
	}

	return(true);
}

int Get_Phenotypic_Score(std::vector<int> state, RegNetwork R)
{
	if(R.phenotype_pos.size() == 0 || R.phenotype_neg.size() == 0)
	{
		std::cout << "I pity the fool who tries to calculate a phenotypic score for this network. I do. I do." << "\n";
		return(-50);
	}
	if(R.phenotype_pos.size() != R.phenotype_neg.size())
	{
		std::cout << "Something is seriosuly broken over here. Go back to genesis, or to the part where the network files were read." << "\n";
		return(-50);
	}

	int score = 0;
	for(int i = 0; i < R.phenotype_pos.size(); i++)
	{
		if(state[R.phenotype_pos[i]] == 1)
		{
			score += 1;
		}
		else if(state[R.phenotype_pos[i]] == 0)
		{
			score -= 1;
		}

		if(state[R.phenotype_neg[i]] == 1)
		{
			score -= 1;
		}
		else if(state[R.phenotype_neg[i]] == 0)
		{
			score += 1;
		}
	}

	return(score);
}

void OutputJij(RegNetwork R)
{
  // Prints out Jij matrix in sparse form
  std::ofstream f("networkfiles/Jij");
  for (int i=0;i<R.numedges;i++)
    {
      f << R.interactions[i][0] << " " << R.interactions[i][1] << " " << 3 - 2*R.interactions[i][2] << "\n";
    }
  f.close();
}
      
      
