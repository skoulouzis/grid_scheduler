/***************************************************************************
 *   Copyright (C) 2007 by S. Koulouzis, T. Wood   *
 *   skoulouz@science.uva.nl, twood@science.uva.nl  *
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
#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "definitions.h"
//#include <fstream>

class Grid {
	
public:
	/**
	constructor initialises grid environment 
	creates numResources resources and sets their various attributes
	@param double
	*/
	Grid(double numResources);
	
	
	/** Grid destructor */
	~Grid();
	
//	/**
//	returns a vector containing the mips rating of each resource
//	*/
//	vector<double> getMipsRatings();
	
	
	/**
	returns ResourceList to scheduler
	@return ResourceList
	*/
	ResourceList getResourceList();
	
	
	/**
	return resource id of resource at index in vector	
	@param int
	@return int 
	*/
	int getId(int index);
	
	
	/**
	returns mips rating of rsource with id
	@param int
	@return double 
	*/
	double getComputationalPotential(int id);
	
	
	/**
	recieves schedule from scheduler and attempts 
	to allocate tasks to resources specified therein...
	
	also recieves vector of of values to represent the mi requirmenst of each
	task which willbe added to the specified resources load upon sucessful submission	
	
	returns int vector representing task allocation results	
	@param Schedule
	@param vector<double>
	@return vector<double>
	*/
	void executeSchedule(Schedule schedule, vector<double> miReq);

	
	/**
	get the current load of resource with id
	@param int
	@return double
	*/
	double getLoad(int id);
	
	
	/**
	get the current availibilityof this resource
	@param int
	@return bool
	*/
// 	double getAvailable(int id);
	bool getReliability(int id);
	
	
//	/**
//	set the availibility of  resource with id
//	@param int
//	@param	double
//	*/
// 	void setAvailable(int id, double available);

	/**
	set the reliability of  resource with id
	@param int
	@param	double
	*/
	void setReliability(int id, double available);
	
	
	/**
	Updates resource load data. Simulates the passing of one second, at the moment
	by subtracting mips rating from load. The iteartions parameter makes up for the time the scheduler took. The more iterations it took, the more mi are executed.
	@param	iteartions
	 */
	void updateResourceLoadData(int iteartions);
	
	
	/**
	Appends a line of data to the grid output file detailing resource utilzation
	load, mips, waitime, reliability, requests,
	 */
	void appendResourceUtilData();
	
private:
	ResourceList  gridResourceList;
	ofstream output;
	ofstream resourceOutput[1000];
	int width;
	int count;
	double totalMips;
	double totalLoad;
	double resourceUtiliziation;
	double thisMI[1000];//^
	double thisMIExecuted[1000];//^
};
#endif

