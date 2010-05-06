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

#ifndef STDDEFINITIONS_H
#define STDDEFINITIONS_H

#include "mpi.h"
#include <vector>
#include <iostream> 
#include <string>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <math.h> 
#include <map>
#include <iomanip>
#include <sys/time.h>

#define MAX_INT numeric_limits<int>::max()
#define MIN_INT numeric_limits<int>::min()

#define MAX_DOUBLE numeric_limits<double>::max()
#define MIN_DOUBLE numeric_limits<double>::min()
/**
This header file includes all the nessary definition for the SESP, like population, resource, task, etc, to make life easy

@author S. Koulouzis, T. Wood <skoulouz@science.uva.nl>
*/


using namespace std;

/**
  * struct Task
  * A task abstratction. The task is asumed to be atomic (not part of a job), and
  * independent of any other task
  */
struct Task {
	// /**
	// The ID (not sure if it's necessary)
	// */
// 	int id;
	/**
	the task's reqirements in milions of instuctions (mi), it will be used to calculate fitness of a schedule 
	*/
	double mi;
	/**
	for future development
	*/
// 	int dependsOn;		
	/**
	for future development
	*/
// 	int ProvidesTo;
};


/**
  * struct Recource
  * This a resource abstarction. It can be assigned with tasks, which are comleted
  * acoarding to the resource's computational power, and the task's requrements
  */
struct Gresource {
	int id;
	/**
	The computational power of this resource, in mpis, it will be used to calculate fitness of a schedule
	*/
	double mips;
	
	/**
	The number of tasks assigned in this resource. Not sure if it's necessary as the load can provide the nessary information

	/**
	keeps count of the umber of times this resource has been requested by a scheduler
	in a schedule.. (regardless of submission sucess)..)
	*/
	int requests;
	
	/**
	The total load of this resource, in mi (sum of mi requirements of all tasks allocated to this resource)
	*/
	double load;

	/**
	Double value 0 < 1 used to determine the likelyhood that the submission of a task is sucessful.
	e.g where the available value is 0.1 there is a 10% probability that task submission will be sucesful
	whereas a resource with available value 1 will always accept tasks, regardless of current workload.
   	*/
	double reliability;
// 	vector<int> taskQueue; //The tasks assigned in this recource
};

/**
This is an encoding of a schedule/individual. It is just a vector of ints that encodes resurces in the content, and tasks as the index. eg schedule[0,3,6,4] means that resurce 0 is assigned with task 0, resource 3 with task 1 etc.
*/
typedef vector<int> Schedule; 

typedef vector<Gresource> ResourceList;

typedef vector<Task> TaskList;

//return type from grid after submitting a schedule
typedef vector<int> ResultVec;

#endif
