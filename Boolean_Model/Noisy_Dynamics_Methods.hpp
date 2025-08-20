double CalculateEnergy(std::vector<int> state, RegNetwork R)
{
	double energy = 0.0;
	int J_ij, s_i, s_j;
	for(int i = 0; i < R.numedges; i++)
	{
		s_i = state[R.interactions[i][0]];
		s_j = state[R.interactions[i][1]];
		J_ij = R.interactions[i][2];

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

		energy += J_ij*s_i*s_j;
	}

	return(-energy);
}

double CalculateDEnergyH(std::vector<int> state, RegNetwork R,int node)
{
  int nedges=R.node_edges[node].size();
  int J_ij, s_i, s_j;
  double energy=0;
  // std::cout << "Enter DEng " << node << " nedges: " << nedges <<  "\n";
  for (int edge=0;edge < nedges;edge++)
  {
    int i=R.node_edges[node][edge];
    s_i = state[R.interactions[i][0]];
    s_j = state[R.interactions[i][1]];
    J_ij = R.interactions[i][2];
    s_i = 2*s_i - 1; //0 -> -1
    s_j = 2*s_j - 1;
    J_ij = 3 - 2*J_ij; //2-> -1
    energy += J_ij*s_i*s_j;
    //  std::cout << i << " " << s_i << " " <<  s_j << " " << J_ij << " " << energy << "\n";
  }
  energy += 2*R.H[node]*(2*state[node]-1);
  return(2*energy);
}

double CalculateDEnergy(std::vector<int> state, RegNetwork R,int node)
{
  int nedges=R.node_edges[node].size();
  int J_ij, s_i, s_j;
  double energy=0;
  // std::cout << "Enter DEng " << node << " nedges: " << nedges <<  "\n";
  for (int edge=0;edge < nedges;edge++)
  {
    int i=R.node_edges[node][edge];
    s_i = state[R.interactions[i][0]];
    s_j = state[R.interactions[i][1]];
    J_ij = R.interactions[i][2];
    s_i = 2*s_i - 1; //0 -> -1
    s_j = 2*s_j - 1;
    J_ij = 3 - 2*J_ij; //2-> -1
    energy += J_ij*s_i*s_j;
    //  std::cout << i << " " << s_i << " " <<  s_j << " " << J_ij << " " << energy << "\n";
  }
  return(2*energy);
}

std::vector<int> Metropolis5(double T, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *g)
//  Metropolis with CLAMPING AND PARALLEL UPDATE and output of Pheno Score to g
{
  //std::cout << "Enter Metropolis3 \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(state);
	}

  /*
    std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */
  // Create array of nodes to flip or not
  std::vector<int> flipnodes;
  std::vector<double> vDelE;
  for (int ii=0; ii<R.numnodes; ii++) {
    flipnodes.push_back(0);
    vDelE.push_back(0.0);
  }
    
	// numsteps=3;
	int nodetoupdate, flag = 0;
	// double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	int sc;
	int tf=0;
	for(int t = 0; t < numsteps; t++)
	{
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
		    nodetoupdate = i;
		    vDelE[i]=0.0;
		    flipnodes[i]=0;
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				// double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
		    if (R.Clamp[nodetoupdate]  == 0)
		      {
			double DeltaE1 = CalculateDEnergy(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
		        flag = 0;
		        vDelE[i]=DeltaE1;
		        if(DeltaE1 < 0.0)
			  {
				flag = 1;
			  }
		        else
			  {
		            if(distribution(generator) < std::exp(-DeltaE1 / T))
			      {
					flag = 1;
		              }
			  }
		       flipnodes[i]=flag;
		     }
		}
		for (int i=0;i<R.numnodes;i++) {
		  if (flipnodes[i]==1) {
		    state[i]=1-state[i];
		  }
		}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		
		
		sc=Get_Phenotypic_Score(state, R);
		if (sc <= -8) {
		  if (tf == 0) {
		    tf = t + 20;
		    //		    std::cout << "tf Set: " << tf << "\n";
		  }
		  else if (t > tf) {
		    return(state);
		  }
		}
		else {
		  tf = 0;	
		}
		(*g) << t << " " <<  sc <<  " ";
		for (int k=0;k<R.numnodes;k++) {
		  (*g) << state[k]+1 << " ";
		}
		(*g) << "\n";
	}
	return(state);
}

std::vector<int> Metropolis4a(double T, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *g)
//  Metropolis with CLAMPING and output of every spin flip and Pheno Score to g
{
  // std::cout << "Enter Metropolis4 \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(state);
	}

  /*
    std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */

	// numsteps=3;
	int nodetoupdate, flag = 0;
	// double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	int sc;
	int tf=0;
	for(int t = 0; t < numsteps; t++)
	{
	  std::cout << t << " ";
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				// double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
			if (R.Clamp[nodetoupdate]  == 0)
			  {
			double DeltaE1 = CalculateDEnergy(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
			flag = 0;
			if(DeltaE1 < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE1 / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = 1-state[nodetoupdate];
				sc=Get_Phenotypic_Score(state, R);
				(*g) << nodetoupdate << " " << sc << " " << t << "\n";
			}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		}
		}
		sc=Get_Phenotypic_Score(state, R);
		if (sc <= -8) {
		  if (tf == 0) {
		    tf = t + 20;
		    //		    std::cout << "tf Set: " << tf << "\n";
		  }
		  else if (t > tf) {
		    return(state);
		  }
		}
		else {
		  tf = 0;	
		}
	}
	return(state);
}


std::vector<int> Metropolis4(double T, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *g)
//  Metropolis with CLAMPING and output of Pheno Score to g
{
  // std::cout << "Enter Metropolis4 \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(state);
	}

  /*
    std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */

	// numsteps=3;
	int nodetoupdate, flag = 0;
	// double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	int sc;
	int tf=0;
	for(int t = 0; t < numsteps; t++)
	{
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				// double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
			if (R.Clamp[nodetoupdate]  == 0)
			  {
			double DeltaE1 = CalculateDEnergy(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
			flag = 0;
			if(DeltaE1 < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE1 / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = 1-state[nodetoupdate];
			}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		}
		}
		sc=Get_Phenotypic_Score(state, R);
		if (sc <= -8) {
		  if (tf == 0) {
		    tf = t + 20;
		    //		    std::cout << "tf Set: " << tf << "\n";
		  }
		  else if (t > tf) {
		    return(state);
		  }
		}
		else {
		  tf = 0;	
		}
		(*g) << t << " " <<  sc <<  " ";
		for (int k=0;k<R.numnodes;k++) {
		  (*g) << state[k]+1 << " ";
		}
		(*g) << "\n";
	}
	return(state);
}

std::vector<int> Metropolis3(double T, int numsteps, std::vector<int> state, RegNetwork R)
//  Metropolis with CLAMPING
{
  //std::cout << "Enter Metropolis3 \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(state);
	}

  /*
    std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */

	// numsteps=3;
	int nodetoupdate, flag = 0;
	// double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	for(int t = 0; t < numsteps; t++)
	{
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				// double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
			if (R.Clamp[nodetoupdate]  == 0)
			  {
			double DeltaE1 = CalculateDEnergy(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
			flag = 0;
			if(DeltaE1 < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE1 / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = 1-state[nodetoupdate];
			}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		}
		}
	}
	return(state);
}

std::vector<int> Metropolis2(double T, int numsteps, std::vector<int> state, RegNetwork R)
//   Metropolis WITH H Field
{
  //std::cout << "Enter Metropolis2 \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(state);
	}

  /*
    std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */

	// numsteps=3;
	int nodetoupdate, flag = 0;
	// double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	for(int t = 0; t < numsteps; t++)
	{
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				// double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
			
			double DeltaE1 = CalculateDEnergyH(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
			flag = 0;
			if(DeltaE1 < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE1 / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = 1-state[nodetoupdate];
			}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		}
	}
	return(state);
}

std::vector<int> Metropolis1(double T, int numsteps, std::vector<int> state, RegNetwork R)
{
  //std::cout << "Enter Metropolis1 \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(state);
	}

  /*
    std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */

	// numsteps=3;
	int nodetoupdate, flag = 0;
	// double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	for(int t = 0; t < numsteps; t++)
	{
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				// double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
			
			double DeltaE1 = CalculateDEnergy(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
			flag = 0;
			if(DeltaE1 < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE1 / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = 1-state[nodetoupdate];
			}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		}
	}
	return(state);
}


		



double Get_Lifetime(double T, std::vector<int> state, RegNetwork R)
{
  int time=0;
  int f=100;
  const int dt=1;
  while (f > 10)
    {
      state=Metropolis1(T, dt, state, R);
      f=CalculateFrustration(state, R);
      time += dt;
    }
  return(time);
}

double Get_Lifetime_Limit(double T, std::vector<int> state, RegNetwork R,int limit)
{
  int time=0;
  int f=100;
  const int dt=1;
  while (f > 10 && time<limit)
    {
      state=Metropolis1(T, dt, state, R);
      f=CalculateFrustration(state, R);
      time += dt;
    }
  return(time);
}
    
  
  
    
 


void Metropolis(double T, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *f, std::ofstream *g)
{
  std::cout << "Enter Metropolis \n";
  if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return;
	}

  /*
	std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}
  */

	// numsteps=3;
	int nodetoupdate, flag = 0;
	//	double DeltaE;
	//	std::cout << R.numnodes << " " << state.size() << "\n";
	for(int t = 0; t < numsteps; t++)
	{
	  /*
	  for(int i=0; i<state.size();i++)
	    {
	      std::cout << state[i] << " ";
	      if (abs(2*state[i]-1) != 1)
		{
		  std::cout << "State BUG " << state[i] << "\n";
		  return;
		}
	    }
	  std::cout << "\n";
	  */
		for(int i = 0; i < R.numnodes; i++)
		{
		  //	  std::cout << t << "flip " << i << "\n";
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			/*
				for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
				//	double Eold=CalculateEnergy(state, R);
				//	double Enew=CalculateEnergy(newstate, R);
			
				// DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
				//		std::cout << "Shubham "  << Eold << " " << Enew << " " << DeltaE << "\n";
				*/
			double DeltaE1 = CalculateDEnergy(state, R, nodetoupdate);
			//	std::cout << "Deltas " << DeltaE1 << " " <<  DeltaE << "\n";
			/*
			if (DeltaE != DeltaE1)
		  {
		    std::cout << "BUG " << DeltaE << " " << DeltaE1 << " " << nodetoupdate <<  "\n";
		    } 
			*/ 
		
			flag = 0;
			if(DeltaE1 < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE1 / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = 1-state[nodetoupdate];
			}
			//			std::cout << "Del NewE " << DeltaE1 << " " <<  CalculateEnergy(state,R) << "\n";
		}


		

		if(t > numsteps / 2 && t % 50 == 0)
		{
			(*f) << CalculateFrustration(state, R) << " ";
			(*g) << Get_Phenotypic_Score(state, R) << " ";
		}
	}
	(*f) << "\n";
	(*g) << "\n";

	return;
}

std::vector<int> Noisy_Asyn_Dynamics(double eta, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *f, std::ofstream *g)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(std::vector<int>());
	}

	int errflag = 0;

	int nodetoupdate, input;

	int source, target, interactiontype;
	for(int runcount = 0; runcount < numsteps; runcount++)
	{
		for(int step = 0; step < R.numnodes; step++)
		{
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
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
					if(state[source] == 0)
					{
						input += (1)*(-1);
					}
					else if(state[source] == 1)
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
					if(state[source] == 0)
					{
						input += (-1)*(-1);
					}
					else if(state[source] == 1)
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
			if(errflag == 1)
			{
				std::cout << "Wicked weird stuff going on here. Take some time and think about what you have done." << "\n";
				return(std::vector<int>());
			}
			if(input > 0)
			{
				state[nodetoupdate] = 1;
			}
			else if(input < 0)
			{
				state[nodetoupdate] = 0;
			}
			else
			{
				state[nodetoupdate] = state[nodetoupdate];
			}

			if(distribution(generator) < eta)
			{
				state[nodetoupdate] = (int) (!state[nodetoupdate]);
			}
		}

		if(runcount > numsteps / 2 && runcount % 50 == 0)
		{
			(*f) << CalculateFrustration(state, R) << " ";
			(*g) << Get_Phenotypic_Score(state, R) << " ";
		}
	}
	(*f) << "\n";
	(*g) << "\n";

	return(state);
}
