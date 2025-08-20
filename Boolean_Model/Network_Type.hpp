class RegNetwork
{
	public:
		bool isvalid = false;
		int numnodes = 0, numedges = 0;
		std::vector<std::string> nodeIDs;
		std::vector< std::vector<int> > interactions;
		std::vector<int> phenotype_pos, phenotype_neg;
                std::vector<std::vector<int>> node_edges;
                std::vector<std::vector<int>> target_edges;
                std::vector<double> H;
                std::vector<int> Clamp;
};
