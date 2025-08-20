// This version also prints out the "initial" state, to correlate its fate
//       with the init. genome
// This version prints out the phenotypic score during the clamping
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <vector>
//#include <mpi.h>
#include "Network_Type.hpp"
#include "Process_Topology_Methods_Quiet.hpp"
#include "Boolean_Modeling_Methods.hpp"
#include "Noisy_Dynamics_Methods.hpp"

void Quench(double Temp, int Time0, std::string filename, int network_id, int trials, std::vector<int> Clamps)
{
  //RegNetwork R = 	ReadTopologyFile(filename + "_" + std::to_string(network_id) + ".ids", filename + "_" + std::to_string(network_id) + ".topo", filename + "_" + std::to_string(network_id) + ".phs");
  //  std::cout <<  filename << '\n';
      	RegNetwork R = 	ReadTopologyFile(filename + std::to_string(network_id) + "network.ids", filename  + std::to_string(network_id) + "network.topo", filename + std::to_string(network_id) + "network.phs");
	if(!R.isvalid)
	{
		std::cout << "Network input failed. Your input files disrespected Tuco Salamanca." << "\n";
		return;
	}
	

	std::vector<int> state,state0;
	std::string str="outputfiles/RUN" + std::to_string(network_id) + "/TRANS10_" + std::to_string(Temp) + "_" + std::to_string(Time0);
	for (int i=0; i<Clamps.size(); i++)
	  {
	    str=str + "_" + std::to_string(Clamps[i]);
	  }
	str=str + ".dat";
	std::ofstream f(str);
	//	std::cout << "OutputFile: " << str << "\n";
	
	str="outputfiles/RUN" + std::to_string(network_id) + "/TRANS10_" + std::to_string(Temp) + "_" + std::to_string(Time0);
	for (int i=0; i<Clamps.size(); i++)
	  {
	    str=str + "_" + std::to_string(Clamps[i]);
	  }
	str=str + "_score.dat";
	std::ofstream g(str);
	//	std::cout << "Score OutputFile: " << str << "\n";
	if(!static_cast<bool>(f))
	{
		std::cout << "The output file specified is AWOL. Go figure." << "\n";
		return;
	}
	

	int fr0=100;
	int fr,sc0;
	const int dt=10;
	

	//	std::cout << "Generating Clamp List" << "\n";
	// numbers correspond to node number (starting at 1!!!!)
	// pos=clamp on;  neg=clamp off
      
	  for (int i=0; i<R.numnodes; i++)
	    {
	      R.Clamp.push_back(0);
	      //	      std::cout << i << " " << R.H[i] << "\n";
	    }
	  //  std::cout << "Number of Clamps: " << Clamps.size() << "\n";
	  for (int i=0; i<Clamps.size(); i++)
	    {
	      R.Clamp[std::abs(Clamps[i])-1]=1;
	      //	      std::cout << "Clamp i: " << abs(Clamps[i])-1 << "\n";
	    }

	 
	  //	  std::cout << "trials: " << trials << "\n";
	  int succeed=0;
	for (int i = 0; i < trials; i++)
	{
	  sc0 = -12;
	  // start from low frustration state
	  while (sc0 < 5)
	  {
	    fr0 = 100;
	    state0 = GetRandomState(R.numnodes);
       	    while (fr0 > 10)
	    {
	      state0=Metropolis1(Temp, 5, state0, R);
	      fr0=CalculateFrustration(state0, R);
	    }
	   /* Only Look at One Type of Init. Cond. */
	    sc0=Get_Phenotypic_Score(state0, R);
	  }
	  // Apply Clamping
	  for (int j=0; j<Clamps.size(); j++)
	    {
	      if (Clamps[j] < 0)
		{
		  state0[std::abs(Clamps[j])-1]=0;
		}
	      else
		{
		  state0[Clamps[j]-1]=1;
		}
	    }
	  //  (f) << fr0 << "  " << sc0 <<  " ";
	  /*
	    for (int i=0;i<R.numnodes;i++)
	    {
	      // std::cout << state0[i] << " ";
	      (f) << state0[i] << " ";
	    }
	  */
	  // std::cout << "\n";




	  //state = Metropolis4(Temp, Time0, state0, R, &g);
	  state = Metropolis4(Temp, Time0, state0, R, &g);
	  
	  
	  /*	  for (int i=0;i<R.numnodes;i++)
	    {
	      std::cout << state[i] << " ";
	    }
	    std::cout << "\n"; */
	  fr=CalculateFrustration(state, R);
	  int sc=Get_Phenotypic_Score(state, R);


	  // Turn off Clamping - Quench at T=0.5

	  state = Metropolis1(0.5,10, state, R);
	  int frf=CalculateFrustration(state, R);
	  int scf=Get_Phenotypic_Score(state, R);
          
	  if (scf<-5 && sc0 >5)
	    {
	      succeed++;
	    }

	  //	  std::cout << i << " Before : " << fr0 << "  " << sc0 << " After: " << fr << "  " << sc << " Final: " << frf << " " << scf << "\n";

	  //  (f)  << fr << "  " << sc << " " <<  frf << " " << scf <<  "\n";	 	
	}

	f.close();
	std::cout << (double)succeed/(double)trials << "\n";

		

	return;
}

int main(int argc, char *argv[])
{
//	MPI_Init(NULL, NULL);
        int world_size;
//      MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int world_rank = 0;
//      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	std::vector<int> Clamps;

	InitializeRandomNumGen(std::time(NULL) + world_rank);

	int NClamp=argc - 6;

       	std::string filename = argv[1];
	int network_id = std::stoi(argv[2]);
	double Temp = std::stod(argv[3]);
	int Time0 = std::stoi(argv[4]);
	int trials = std::stoi(argv[5]);
	for (int i=0;i<NClamp;i++)
	  {
	    Clamps.push_back(std::stoi(argv[6+i]));
	  }


	//	std::cout <<  filename << '\n';

	//	Sample_States(T, Tid, "inputfiles/" + filename, network_id, world_rank, 500);
	Quench(Temp, Time0, "networkfiles/" + filename, network_id, trials, Clamps);

//	MPI_Finalize();
	
	return(0);
}
