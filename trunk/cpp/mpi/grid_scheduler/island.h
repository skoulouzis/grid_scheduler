/***************************************************************************
 *   Copyright (C) 2007 by S. Koulouzis, T. Wood   *
 *   skoulouz@science.uva.nl   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef ISLAND_H
#define ISLAND_H

/**
This is a schduling implementation of the island model. Each islend holds sub-set of the population, in which it applyes a GA. This class has a migration operation to enable for mixing of populations/solution

	@author S. Koulouzis, T. Wood <skoulouz@science.uva.nl>
*/

#include "definitions.h"
#include "geneticalgorithm.h"

class Island{
public:
	/**
	An emty constructor, that will create an island.
	*/
	Island(int populationSize,int popInitMethod,int fitnessType, int myrank);
	
	~Island();	
	
	
//	/**
//	This method initializes the population in this island. In the parallel version the population is a sub-set of the global population (a //sub-set of solutions to evaluate). The batchSize sets the size of an individual schedule, that is used in initializing the population, and //running the GA, while the mipsRatings and miRequirements are ther to be used by the GA
//	@param  int
//	@param  int
//	@param ResourceList
//	@param vector<double>
//	*/
	//void initialize(ResourceList resourceList,vector<double> miRequirements);
	
	//	/**
	// set the mpis rating of resources (a vector of doubles) for population initialization, and for the // fitness function in the GA
	// @param mipsRatings 
	// */
// 	void setResourceMipsRatings(vector<double> mipsRatings);
	
	// /**
	// returns the mpis rating of resources seen by this island
	// @return vector<double> 
	// */
// 	vector<double>  getResourceMipsRatings();
	
	
	/**
	set the mi rating of tasks (a vector of doubles) for population initialization, and for the fitness function in the GA
	@param mipsRatings 
	*/
	void setTaskMiRatings(vector<double> miRatings);
	
	/**
	returns the mpis rating of tasks seen by this island
	@return vector<double> 
	*/
	vector<double>  getTaskMiRatings();
	
	
	/**
	* After a sertant amount of iterations, the island is migrating some of its
	* population to its neighbor. For MPI this sould be an array for esyer
	* communication.
	* @param vector<Schedule> 
	* @return vector<Schedule>
	*/
	vector<Schedule> migrate(vector<Schedule> aPopulation);
	
	
//	/**
//	* An MPI recive operation for reciving population from neighbor islands (acoarding
//	* to the MPI topology)
//	* @return vector<Schedule>
//	*/
// 	vector<Schedule> receiveImigrants();
	
	/**
	This method starts running a genetic algorithm, based on the population created by the initialize method.After the GA is done it will return a "good" population, optimized by execution time, of size n. Make sure that the initialize method is called before calling this method.
	@param int 
	@return vector<Schedule>
	*/
	Schedule runGA(ResourceList resourceList,vector<double> miRequirements);
	
	/**
	Out puts data to a file
	*/
	void outputData();
	
	/**
	returns the number of iterations it took the GA to return a schedule
	@return int
	*/
	int getIterations();
		
	/**
	returns the fitness of a schedule specified by the fitnessType
	@param Schedule
	@param int
	@return double
	*/
	double getFitness(Schedule aSchedule,int fitnessType);
	
	/**
	sets the rank of this isleand. Only for island model use, with migration (MPI)
	@param int
	*/
	void setMyRank(int rank);
	
	/**
	retunts this iseland's rank. Only for island model use, with migration (MPI)
	@return int
	*/
	int getMyRank();
	
	/**	
	Sets the communicator of this island, so migration and stopping conditions may be used.
	@param MPI_Comm
	*/
	void setCommunicator(MPI_Comm communicator);
	
	/**
	Sets the left and right neighbors for this iseland, in order to know where to migrate
	@param int
	@param int 
	*/
	void setNeighbors(int left,int right);

	/**
	returns the initilization time in sec
	@return double 
	*/
	double getT_init();
	
	/**
	returns the computation time in sec
	@return double 
	*/
	double getT_comp();

	/**
	returns the communication time in sec
	@return double 
	*/
	double getT_comm();

	/**
	returns the comutation time of the GA in sec
	@return *double 
	*/
	double* getT_compGA();
	
	/**
	returns the fitness evaluation time in sec
	@return double 
	*/
	double getT_compFitnessEval();
	
private:
	/**
	Initializes the population randomly 
	*/
	void randomInit();
	
	/**
	Returns a vector of doubles representing the mips ratings of the
	resources listed in resourceList.
	@return vector<double>
	*/
	vector<double> getResourceMipsRatings();

	/**
	Returns a vectorof doubles representing the current loads of the
	resources listed in resourceList...
	@return vector<double>
	*/
	vector<double> getResourceLoads();
	
	
	/**
	Allocate tasks to resources depending on task requirements	
	Tasks are either assigned to task by shortest task to fastest resource
	or if 'option' is set to '1' longest taskto fastest resource... 
	@param int 
	*/
	void taskMi2fastest(int option);

	/**
	Starts counting time 
	*/
	void timeStart();

	/**
	Stops counting time,  could it is communication, computation etc.
	@param int 
	*/
	void timeEnd(int what);
	
private:
	//int populationSize;
	int batchSize;
	int populationSize;int popInitMethod;int fitnessType;
	vector<double> resourceMipsRatings;
	vector<double> resourceLoads;
	ResourceList resourceList;
	vector<double> taskMiRatings;
	vector<Schedule> population;
	GeneticAlgorithm *ga;
	int randomIndex;
	int numOfResources;
	ofstream output;
	int width;
	int iterations;
	double flowTime, waitTime;
	double executionTime;
	//MPI stuff
	int myrank;
	MPI_Comm communicator;
	int myLeft,myRight;
	double sigma;
	double maxdif;
	double t_init;
	double t_comp;
	double t_comm;
	double t_compFitnessEval;
	double t_compGA[5];
	double td;
	long sec,usec;
	struct timeval tp;
};

#endif
