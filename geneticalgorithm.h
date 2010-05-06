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
#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

/**
This is a genetic algorithm imlementation. It provides all the operations of a GA like selection, mutation etc

	@author S. Koulouzis, T. Wood <skoulouz@science.uva.nl>
*/

#include "definitions.h"
#include <math.h>
#include <time.h>

class GeneticAlgorithm{

public:
	/**
	Constructor of the GA 
	@param vector<Schedule>
	*/
	GeneticAlgorithm(vector<Schedule> aPopulation);
	
	GeneticAlgorithm();
	
	~GeneticAlgorithm();
	
	/**
	Selects 2 parents for the first step of the GA, taken from the population. The selection might be random,or according to some heuristic. If the selection is a tournament selection then this operation should select a sub set of the population, and return the 2 best. The different kinds of selection are specified according the kindOfSelection int
	Note this operation reruns only 2 parents.	
	@return vector<Schedule> 
	@param kindOfSelection
	*/
	vector<Schedule> parentsSelection(int kindOfSelection);
	
	
	/**
	It takes a vector of schedules, and eliminates the "weak" ones, by calculating their fitness, which is selected acoarding to the @param int. This method should keep the population size constant
	@param candidates
	@param int
	*/
	void naturalSelection(vector<Schedule> candidates,int kindOfFitness);
	
	/**
	It takes candidates, and evaluates their fitness. So better individuals get higher chance of not being eliminted. Chances proportional to fitness. Assign to each individual a part of the roulette wheel,spin the wheel n times to select n individuals. This method should keep the population size constant
	@param candidates
	@param int
	*/
	void rouletteWheelSelection(vector<Schedule> candidates,int kindOfFitness);
		
	/**
	* This is a corossover between 2 parrents. The crossover might be random, cyclic,
	* or any other kind. The kind of crossover is determined by kindOfXover It returns 2 new offsprings
	* @return vector<Schedule> 
	* @param  parents
	* @param  kindOfXover
	*/
	vector<Schedule> crossover (vector<Schedule> parents, int kindOfXover);
	

	/**
	* Takes an individual as an input, mutates it, and reurns it. The mutation migth
	* be swaping geens from a random, schedule/individual or random change from a
	* resource list. Alter each gene independently with a probability pm 
	* pm is called the mutation rate, typically between 1/pop_size and 1/ chromosome_length
	* @return Schedule
	* @param  aSchedule
	* @param  rateOfMutation
	*/
	Schedule mutation(Schedule aSchedule,double rateOfMutation);
	
	
	Schedule mutationWithHeuristic(Schedule aSchedule,double rateOfMutation);
	
	/**
	* it reurns the fitness of an individual, acoarding what needs to be measured
	* @return double
	* @param  aSchedule
	* @param  kindOfFitness
	*/
	double fitnessFunction (Schedule aSchedule, int kindOfFitness);
	
	/**
	returns the average fitness of the GA's population, according to a specified fitness function
	@param  double
	@return double
	*/
	double getAveragefitness(int kindOfFitness);
	
	/**
	set the initial population for this GA. 
	@param vector<Schedule>
	*/
	void setPopulation(vector<Schedule> aPopulation);
	
	/**
	get The current population
	@return vector<Schedule> 
	*/
	vector<Schedule> getPopulation();
	
	/**
	returns the number of resources seen by this GA
	@return int
	*/
	int getNumOfResources();
	
	/**
	sets the number of resources to this GA. It is used in mutation, to swap a geen( resource id ) in a schedule,with a random resource. 
	@param int
	*/
	void setNumOfResources(int num);
	
	/**
	set the mpis rating of resources (a vector of doubles) for use in the fitness function 
	@param mipsRatings 
	*/
	void setResourceMipsRatings(vector<double> mipsRatings);
	
	/**
	returns the mpis rating of resources seen by this GA
	@return vector<double> 
	*/
	vector<double>  getResourceMipsRatings();
	
	
	/**
	set the mi rating of tasks (a vector of doubles) for use in the fitness function 
	@param mipsRatings 
	*/
	void setTaskMiRatings(vector<double> miRatings);
	
	/**
	returns the mpis rating of tasks seen by this GA
	@return vector<double> 
	*/
	vector<double>  getTaskMiRatings();
	
	
		/**
	set the mi rating of tasks (a vector of doubles) for use in the fitness function 
	@param mipsRatings 
	*/
	void setResourcesLoad(vector<double> resourcesLoad);
	
	/**
	returns the mpis rating of tasks seen by this GA
	@return vector<double> 
	*/
	vector<double>  getResourcesLoad();
	
	/**
	returns n best schedules, ordered by their fitness, which is specified by the kinOfFitness parameter.
	@param int
	@param int
	@return vector<Schedule>
	*/
	vector<Schedule> getBestSchedules(int n,int kinOfFitness);
	
	
	/**
	returns the n best schedule's indexs, ordered by their fitness, which is specified by the kinOfFitness parameter.
	@param int
	@param int
	@return vector<int>
	*/
	vector<int> getIndexOfBestSchedules(int n,int kinOfFitness);
	
	/**
	removes n individuals from the GA's population
	@param int
	*/
	void removeIndividuals(int n);
	
	/**
	Adds an individual to tha GA's population
	@param Schedule
	*/
	void addIndividual(Schedule aSchedule);
	
	/**
	removes an individual from the GA's population according to its index
	@param int
	 */
	void removeIndividualsByIndex(int index);
	
private:
	/**
	Random selection. It takes 2 random int and chooses the parent according to those
	@return vector<Schedule> 
	*/
	vector<Schedule> randomParentsSelection();
	
	/**
	1-point crossover. 
	Choose a random point on the two parents
	Split parents at this crossover point
	Create children by exchanging tails
	@param vector<Schedule>
	@return vector<Schedule>
	*/
	vector<Schedule> onePointscrossover(vector<Schedule> parents );
	
	/**
	calculates the total execution time of a schedule
	@param Schedule
	@return double
	*/
	double executionTimeFunction(Schedule aSchedule);
	
	/**
	calculates the total flow time of a schedule
	@param Schedule
	@return double
	*/
	double flowTimeFunction(Schedule aSchedule);
	
	/**
	calculates the total wait time of a schedule, taking the load of resources 
	@param Schedule
	@return double
	*/
	double getWaitTimeFunction(Schedule aSchedule);
	
	/**
	calculates the resource Utilization of a schedule
	@param Schedule
	@return double
	*/
	double resourceUtilizationFunction(Schedule aSchedule);		
	
private:
// 	Schedule parents[2]; 
	vector<Schedule> population;
	vector<Schedule> orderedPopulation;
	int randomIndex;
	double randomNum;
	vector<double> resourceMipsRatings;
	vector<double> taskMiRatings;
	vector<double> resourcesLoad;
// 	vector<double> executionTimes;
	double executionTime, waitTime;
	double mi;
	double mips;
	double resLoad;
// 	vector<double> flowTimes;
	
	/**
	This is used in mutation, where we need to swap on resource ID to an other. So the random number is from 0 to numOfResources
	*/
	int numOfResources;
};

#endif
