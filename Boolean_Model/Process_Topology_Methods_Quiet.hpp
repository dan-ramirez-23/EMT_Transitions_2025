RegNetwork ReadTopologyFile(std::string filename1, std::string filename2, std::string filename3)
{
	RegNetwork R;
	R.numnodes = 0;
	R.numedges = 0;
	R.nodeIDs = std::vector<std::string>();
	R.interactions = std::vector<std::vector<int> >();
	R.phenotype_pos = std::vector<int>();
	R.phenotype_neg = std::vector<int>();
	// An array that associates edges with nodes
	R.node_edges = std::vector<std::vector<int>>();
	R.target_edges = std::vector<std::vector<int>>();
	std::string nodeName;
	int nodeID, numlines;

	std::ifstream f(filename1);
	if(!static_cast<bool>(f))
	{
		std::cout << "Error reading the .ids file." << "\n";
		return(R);
	}
	while(f >> nodeName >> nodeID)
	{
		R.nodeIDs.push_back(nodeName);
		R.numnodes += 1;
	}
	f.close();

	f.open(filename1);
	numlines = 0;
	while(getline(f, nodeName))
	{
		numlines += 1;
	}
	f.close();

	for (int i=0;i<R.numnodes;i++)
	{
	  R.node_edges.push_back({-1});
	  R.target_edges.push_back({-1});
	}

	if(numlines != R.nodeIDs.size())
	{
		std::cout << "The format of the .ids file is ridonc, and frankly, completely wrong." << "\n";
		return(R);
	}

	std::string source, target;
	int sourceID, targetID, interactionType;

	f.open(filename2);
	if(!static_cast<bool>(f))
	{
		std::cout << "Error reading the .topo file." << "\n";
		return(R);
	}
	getline(f, source);
	while(f >> source >> target >> interactionType)
	{
		if(std::find(R.nodeIDs.begin(), R.nodeIDs.end(), source) == R.nodeIDs.end())
		{
			std::cout << "One of the network nodes is missing from the .ids file." << "\n";
			return(R);
		}
		if(std::find(R.nodeIDs.begin(), R.nodeIDs.end(), target) == R.nodeIDs.end())
		{
			std::cout << "One of the network nodes is missing from the .ids file." << "\n";
			return(R);
		}
		sourceID = std::find(R.nodeIDs.begin(), R.nodeIDs.end(), source) - R.nodeIDs.begin();
		targetID = std::find(R.nodeIDs.begin(), R.nodeIDs.end(), target) - R.nodeIDs.begin();
		R.interactions.push_back({0, 0, 0});
		R.interactions[R.numedges][0] = sourceID;
		R.interactions[R.numedges][1] = targetID;
		R.interactions[R.numedges][2] = interactionType;
		//  Update node_edges
		if (targetID != sourceID) //NO SELF-LINKS
		  {
		    if (R.node_edges[sourceID][0] == -1)
		      {
			R.node_edges[sourceID][0]=R.numedges;
		      }
		    else
		      {
			R.node_edges[sourceID].push_back(R.numedges);
		      }
		    /*
		    std::cout << sourceID << "\t";
		    for (int i=0;i<R.node_edges[sourceID].size();i++)
		      {
			std::cout << R.node_edges[sourceID][i] << "\t";
		      }
		    std::cout << "\n";
		    */
		    if (R.node_edges[targetID][0] == -1)
		      {
			R.node_edges[targetID][0]=R.numedges;
		      }
		    else
		      {
			R.node_edges[targetID].push_back(R.numedges);
		      }
		    //    std::cout << targetID << "\t";
		    /*   for (int i=0;i<R.node_edges[targetID].size();i++)
		      {
			std::cout << R.node_edges[targetID][i] << "\t";
		      }
		    std::cout << "\n";
		    */
		  }
		// Update target_edges
		  {
		    if (R.target_edges[targetID][0] == -1)
		      {
			R.target_edges[targetID][0]=R.numedges;
		      }
		    else
		      {
			R.target_edges[targetID].push_back(R.numedges);
		      }
		    //    std::cout << "S " << targetID << "\t";
		    /*  for (int i=0;i<R.target_edges[targetID].size();i++)
		      {
			std::cout << R.target_edges[targetID][i] << "\t";
		      }
	       	    std::cout << "\n";
		    */
		  }
		R.numedges += 1;
	}
	f.close();

	f.open(filename2);
	numlines = 0;
	while(getline(f, source))
	{
		numlines += 1;
	}
	f.close();

	if(numlines - 1 != R.interactions.size())
	{
		std::cout << "The format of the .topo file is ridonc, and frankly, completely wrong." << "\n";
		return(R);
	}

	int nodeid_pos, nodeid_neg;
	f.open(filename3);
	if(!static_cast<bool>(f))
	{
		std::cout << "Phenotype file could not be read for this network. Moving on." << "\n";
	}
	else
	{
		while(f >> nodeid_pos >> nodeid_neg)
		{
			R.phenotype_pos.push_back(nodeid_pos);
			R.phenotype_neg.push_back(nodeid_neg);
		}
	}
	f.close();

	if(R.numnodes == R.nodeIDs.size() && R.numedges == R.interactions.size())
	{
		R.isvalid = true;
	}
	int count=0;
	for (int node=0;node<R.numnodes;node++)
	  {
	    count += R.node_edges[node].size();
	  }
	//	std::cout << "From node_edges: " << count << " Numedges: " << R.numedges << "\n";
	count=0;
	for (int node=0;node<R.numnodes;node++)
	  {
	    if (R.target_edges[node][0]!= -1)
	      {
		count += R.target_edges[node].size();
	      }
	  }
	//    std::cout << "From target_edges: " << count << " Numedges: " << R.numedges << "\n";

	return(R);
}

void Print_Network(RegNetwork R, std::ofstream *f)
{
	(*f) << "Source\tTarget\tType\n";
	for(int i = 0; i < R.interactions.size(); i++)
	{
		(*f) << R.nodeIDs[R.interactions[i][0]] << "\t";
		(*f) << R.nodeIDs[R.interactions[i][1]] << "\t";
		(*f) << R.interactions[i][2] << "\n";
	}

	return;
}
