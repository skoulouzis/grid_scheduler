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
#include "island.h"

Island::Island(int populationSize,int popInitMethod,int fitnessType,int myrank){

	
	this->populationSize = populationSize;
	this->popInitMethod= popInitMethod;
	this->fitnessType = fitnessType;
	this->myrank = myrank;

	t_init=0.0;
	t_comp=0.0;
	t_comm=0.0;
	t_compFitnessEval=0.0;
	td=0.0;
	for(int i=0;i<5;i++){
		t_compGA[i]=0.0;
	}
	long sec,usec;

	ga = new GeneticAlgorithm();
	width = 20;
	char fileName[40];
	
	sprintf(fileName,"GAData_%d.txt",myrank);
	
	output.open(fileName);
	output << left                      // Left-justify in each field
		<< "," << "#Iterations"      // Then, repeatedly set the width
		<< "," << "Flow Time"     // and write some data
		<< "," << "Execution Time"
		<< "," << "Wait Time"
		<< endl;
}


void Island::randomInit(){
	Schedule aSchedule;
	
	for(int i=0;i<populationSize;i++){
		aSchedule.clear();
		for(int j=0;j<batchSize;j++){
			aSchedule.push_back((rand()%numOfResources)+0);
		}
		population.push_back(aSchedule);
	}
}


/**
Allocate tasks to resources depending on task requirements

Tasks are either assigned to task by shortest task to fastest resource
or if 'option' is set to '1' longest taskto fastest resource... 
*/
void Island::taskMi2fastest(int option) {
	
	//temporary containers
	Schedule s;
	vector<int> resources;
	vector<int>::iterator iterator;
	//vector<double> mips;
		
	//first sort miRatings in ascending order
	//i.e. shortest task to fastest resource
	sort(taskMiRatings.begin(), taskMiRatings.end());
	
	//if option is '1' then reverse taskMiRatings
	//so they are in descending order
	//i.e longest task to fastest resource
	if(option == 1) {
		reverse(taskMiRatings.begin(), taskMiRatings.end());
	}
	
	for(int i=0;i<populationSize;i++){
		
		//clear temporary containers
		s.clear();
		resources.clear();		
		
		for(int j=0;j<batchSize;j++){
			//pick a ramdom resource ID
			int resID = (rand()%numOfResources)+0;
			//put the resourceID in temp resources container
			resources.push_back(resID);
			//then get resource mips using resID and put it in temp mips container
			//mips.push_back(resourceMipsRatings.at(resID));
		}
		//now we should have 2 containers, one with resource ID's 
		//and one with their mips ratings - and the batch of tasks as well
		
		for(int k = 0; k < batchSize; k++) {
			
			//need to find longest job and fastest resource
			//double longestTask = -1;//init value
			int fastestResource = -1;
			int index;//the index of the fastest resource in resources
			double fastestMips = 0;
			
			for(int l = 0; l < resources.size(); l++) {
				
				//get the mips rating of resource at 'l'
				double temp = resourceMipsRatings.at(resources.at(l));
				
				//check if it is the fastest, and update if it is
				if(temp > fastestMips) {
					fastestMips = temp;
					//set fastest resource
					fastestResource = resources.at(l);
					index = l;//and index in resources
				}
			}
			//Righto, now we have found the fastest resource in resources
			//so push it into the schedule
			
			s.push_back(fastestResource);
			//and also remove the resources
			iterator = resources.begin() + index;
			resources.erase(iterator);
			
			//and now do it all again
		}
		population.push_back(s);
	}
}


void Island::setTaskMiRatings(vector<double> miRatings){
	this->taskMiRatings = miRatings;
	this->batchSize = miRatings.size();
}

// vector<double>  Island::getResourceMipsRatings(){
// 	return this->resourceMipsRatings;
// }
vector<double> Island::getResourceMipsRatings(){
	
	vector<double> mipsRatings;
	//got through the resourceList, and push the resource's mpis into the vector
	for(int i = 0; i < resourceList.size(); i++) {
		Gresource res = resourceList.at(i);
		mipsRatings.push_back(res.mips);
	}
	numOfResources = resourceList.size();
	return mipsRatings;
}

vector<double> Island::getResourceLoads() {
	
	vector<double> resourceLoads;
	//got through the resourceList, and push the resource's loads into the vector
	for(int i = 0; i < resourceList.size(); i++) {
		Gresource res = resourceList.at(i);
		resourceLoads.push_back(res.load);
	}
	return resourceLoads;
}

// vector<double>  Island::getTaskMiRatings(){
// 	return this->taskMiRatings;
// }

Schedule Island::runGA(ResourceList resourceList,vector<double> miRequirements){
// 	double startwtime;
// 	if (myrank == 0){
// 		startwtime = MPI_Wtime();
// 	}
	timeStart();

	this->resourceList = resourceList;
	// 	this->taskMiRatings = miRequirements;
	this->populationSize = populationSize;
	//This also sets the numOfResources 
	this->resourceMipsRatings = this->getResourceMipsRatings();
	this->resourceLoads = this->getResourceLoads();
	population.clear();
			
// 	setResourceMipsRatings(mipsRatings);
	setTaskMiRatings(miRequirements);
	//Initilize the population acoarding to a heuristic
	
	switch (popInitMethod){
	case 0:
		{
			randomInit();
		}
		break;
	case 1:
		{
			taskMi2fastest(0);
		}
		break;
	case 2:
		{
			taskMi2fastest(1);
		}
		break;
	default:
		break;
	}
	
	ga->setPopulation(this->population);
	ga->setNumOfResources(this->numOfResources);
	ga->setResourceMipsRatings(this->resourceMipsRatings);
	ga->setTaskMiRatings(this->taskMiRatings);
	ga->setResourcesLoad(this->resourceLoads);
	
	/**end of init */
	

	
	
	int n = 1;
	
	vector<Schedule> parents;	//will hold the selscted parents from the population
	vector<Schedule> offsprings;	//will hold the offsprings after crossover
	vector<Schedule> candidates;	//will hold the "family", parents and offsprings
	Schedule individual;		//used for mutation
	double fitness,averageFitness,mutationRate,bestMips,warstMips,randomNum,mutationMax,bestFitness,tempBestFitness;
	maxdif = MAX_DOUBLE;
	sigma = 0;
	vector<Schedule> outgoingImmigrants;
	vector<Schedule> incomingImmigrants;
	vector<Schedule> population;
	
	bool end=false;
	int index1,index2,noChange,migrationNoChange;
	iterations = 0;
	noChange = 0;
	migrationNoChange = 0;
	timeEnd(99);
	timeStart();
	
	//find best and warst schedules, (unrealistic), to normalize the fitness to [0,1], so it may be applyed in the mutationRate
	bestMips = *max_element(resourceMipsRatings.begin( ), resourceMipsRatings.end( ));
	warstMips = *min_element(resourceMipsRatings.begin( ), resourceMipsRatings.end( ));
		
	//get the index of the best resource
	for(int i=0;i<resourceMipsRatings.size();i++){
		if(resourceMipsRatings.at(i)==bestMips){
			index1 = i;
			break;
		}
	}
	
	//get the index of the warst resource
	for(int i=0;i<resourceMipsRatings.size();i++){
		if(resourceMipsRatings.at(i)==warstMips){
			index2 = i;
			break;
		}
	}
	
	//create an unrealistic schedule, all the tasks in the best resource
	Schedule aSchedule;
	for(int j=0;j<batchSize;j++){
		aSchedule.push_back(resourceList[index1].id);
	}
	//get the best fitness
	double min = ga->fitnessFunction(aSchedule,fitnessType);
	//create an unrealistic schedule, all the tasks in the warst resource
	aSchedule.clear();
	for(int j=0;j<batchSize;j++){
		aSchedule.push_back(resourceList[index2].id);
	}
		
	//get the warst fitness
	double max = ga->fitnessFunction(aSchedule,fitnessType);
	mutationMax = max + 1.0;
	timeEnd(1);	
	//start the GA loop
	while(!end){
		timeStart();
		//take the best fitness before appling a GA
// 		bestFitness = ga->fitnessFunction(ga->getBestSchedules(1,fitnessType).at(0),fitnessType);
		//take the average fitness before appling a GA
		averageFitness = ga->getAveragefitness(fitnessType);
		timeEnd(2);
		timeStart();
		candidates.clear();
		//select parents. 0 is random selection
		parents = ga->parentsSelection(0);
		timeEnd(3);
		timeStart();
		//apply crossover. 0 is one point crossover
		offsprings = ga->crossover(parents,0);
		//make the "family"
		candidates.push_back(parents[0]);
		candidates.push_back(parents[1]);
		candidates.push_back(offsprings[0]);
		candidates.push_back(offsprings[1]);
		timeEnd(4);
		timeStart();
		//mutate them accoarding to their fitness
		for(int i=0;i<candidates.size();i++){
			randomNum = double(random()) / (double(RAND_MAX) + 1.0);
			fitness = ga->fitnessFunction(candidates[i],fitnessType);
// 			mutationRate = ((fitness/max) * randomNum );
			mutationRate = ((fitness) + randomNum )/mutationMax;
			individual.clear();
			individual = candidates[i];
			candidates[i] = ga->mutation(individual,mutationRate);			
		}
		timeEnd(5);
		timeStart();
		ga->naturalSelection(candidates,fitnessType);
		timeEnd(6);
		// End of computation now migrate

		if(iterations%30==0){
			int immigrantsSize = ga->getPopulation().size()/2;//(random()%6)+1; //1;
		// 	MPI_Bcast ( &sigma, 1, MPI_DOUBLE, 0,communicator );
			timeStart();
			outgoingImmigrants.clear();
			incomingImmigrants.clear();
			population.clear();
			population=ga->getPopulation();
// 			population=ga->getBestSchedules(immigrantsSize,fitnessType);
						
			for(int i=0;i<immigrantsSize;i++){
				outgoingImmigrants.push_back(population.at(i));
			}		
						
			ga->removeIndividuals(immigrantsSize);
// 			vector<int> indexes =  ga->getIndexOfBestSchedules(immigrantsSize,fitnessType);
// 			for(int i=0;i< indexes.size();i++){
// 				ga->removeIndividualsByIndex(indexes.at(i));
// 			}
			timeEnd(1);
			timeStart();
			incomingImmigrants = migrate(outgoingImmigrants);
			timeEnd(0);
			timeStart();
			for(int i=0;i<outgoingImmigrants.size();i++){
				ga->addIndividual(outgoingImmigrants.at(i));
			}
			timeEnd(1);
		}
		timeStart();
// 		tempBestFitness = ga->fitnessFunction(ga->getBestSchedules(1,fitnessType).at(0),fitnessType);
		double tempDelta = (averageFitness - ga->getAveragefitness(fitnessType));
		double delta=0.0;	
		 
		if(tempDelta > delta){
			delta = tempDelta;
		}
// 		if(myrank==0){
// 			cout << " delta:" << delta << ":averageFitness:" << averageFitness << endl;
// 		}
		timeEnd(1);
		timeStart();
		MPI_Allreduce( &delta, &maxdif, 1, MPI_DOUBLE, MPI_MAX, communicator );
		timeEnd(0);
		timeStart();
		if(maxdif<sigma){
			end = true;
		}
		
		if(iterations%100==0){
			sigma = sigma +  (iterations/9990.1);
		}
		
// 		if(myrank==0 && iterations%10==0){
// 			cout << " delta:" << delta << ":averageFitness:" << averageFitness <<":signa:"<< sigma<< endl;
// 		}

		/** end stopping cond stuff */
		
		iterations++;
		if(iterations>3000){
			end = true;
		}
		timeEnd(1);
	}//end while
// 	timeStart();
	
	vector<Schedule> aPopulation;
	aPopulation.clear();

	aPopulation = ga->getBestSchedules(n,fitnessType); 
	
// 	flowTime = ga->fitnessFunction(aPopulation.at(0),1);
// 	executionTime = ga->fitnessFunction(aPopulation.at(0),0);
// 	waitTime = ga->fitnessFunction(aPopulation.at(0),2);
	
	Schedule s = aPopulation.at(0);
	
// 	outputData();
	
// 	timeEnd(2);
// 	if(myrank==0){
// 		cout << MPI_Wtime() - startwtime << endl; 
// 	}

	return s;
}

int Island::getIterations(){
	return iterations;
}

void Island::outputData(){
	output << left                      // Left-justify in each field
		<< "," << this->iterations      // Then, repeatedly set the width
		<< "," << flowTime     // and write some data
		<< "," << executionTime
		<< ", "<< waitTime
		<< endl;
}


double Island::getFitness(Schedule aSchedule,int fitnessType){
	return ga->fitnessFunction(aSchedule,fitnessType);
}

void Island::setCommunicator(MPI_Comm communicator){
	this->communicator = communicator;
}

void Island::setNeighbors(int left,int right){
	this->myLeft = left;
	this->myRight = right;
}

void Island::setMyRank(int rank){
	this->myrank = rank;
}

vector<Schedule> Island::migrate(vector<Schedule> aPopulation){
	Schedule aSchedule;
	int counter=0;
	vector<Schedule> incomingPopulation;
	MPI_Status status; 
	int allocateSize = aPopulation.size() * aPopulation.at(0).size();
	int populationOut[allocateSize];
	int populationIn[allocateSize];

// 	//put the popolation in the buffer
	for(int i=0;i<aPopulation.size();i++){
		aSchedule.clear();
		aSchedule = aPopulation.at(i);
		for(int j=0;j<aSchedule.size();j++){
			populationOut[counter]=aSchedule.at(j);
			counter++;
// 			cout << " populationOut["<< counter << "]="<<  populationOut[counter] << "aSchedule.at(" << i << ")="<< aSchedule.at(j) << endl;
		}
	}	
	
	MPI_Sendrecv(populationOut, allocateSize, MPI_INT, myLeft, 123, populationIn, allocateSize, MPI_INT, myRight, 123, communicator, &status);
 	
 	//extract incoming population
	counter =0;
	Schedule tempSchedule(aPopulation.at(0).size());
	for(int i=0; i< allocateSize;i++){
		if(counter == (batchSize) ){			
			counter = 0;
			incomingPopulation.push_back(tempSchedule);
		}
		//cout << " populationIn["<< i << "]="<<  populationIn[i] << endl;
		tempSchedule.at(counter) = populationIn[i];
		counter++;
		//add the last schedule
		if(i == (allocateSize) - 1){
			incomingPopulation.push_back(tempSchedule);
		}
	}
	
// 	if(myrank==1){
// 		cout << "--------------" << endl;
// 		for(int i=0;i<incomingPopulation.size();i++){
// 			cout << "schedule " << i << endl;
// 			for(int j=0;j<incomingPopulation.at(0).size();j++){
// 				cout << "	"<< incomingPopulation.at(i).at(j) << endl;
// 			}
// 		}
//  	}

	
	return 	incomingPopulation;
}

// vector<Schedule> Island::receiveImigrants(){
// 	int *populationIn;
// 	MPI_Status status; 
// 	
// 	MPI_Recv(populationIn, 4 * sizeof(int), MPI_INT, myRight,114, communicator, &status);
// 	
// }
void Island::timeStart(){

	if(myrank==0)  
	{
		if(gettimeofday(&tp,NULL)==0){
		/* Record the first point of time */
		sec=tp.tv_sec; usec=tp.tv_usec; 
	
		}
		else
		{
		perror("Can not get Wall-Clock time!\n");
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE); 
		}
	}
}

void Island::timeEnd(int what){
	if(myrank==0)
	{
		if(gettimeofday(&tp,NULL)==0)
			{
			sec=tp.tv_sec-sec; 
			usec=tp.tv_usec-usec; 
			/* Get the second point of time */
			td=sec+usec/1.0e6;
			/* Calculate the time difference and print it out */
			switch(what){
				case 0:
					t_comm = t_comm+td;
				break;
				case 1:
					t_comp = t_comp+td;
				break;
				case 2:
					{
					t_comp = t_comp+td;
					t_compFitnessEval = t_compFitnessEval + td;
					}
				break;
				case 3:
					{
					t_comp = t_comp+td;
					t_compGA[0] = t_compGA[0] +td;
					}
				break;
				case 4:
					{
					t_comp = t_comp+td;
					t_compGA[1] = t_compGA[1] +td;
					}
					
				break;			
				case 5:
					{
					t_comp = t_comp+td;
					t_compGA[2] = t_compGA[2] +td;
					}
				break;
				case 6:
					{
					t_comp = t_comp+td;
					t_compGA[3] = t_compGA[3] +td;
					}
					
				break;
				default :
					t_init =t_init + td;
				break;
			}
			
			}
		else
			{
			perror("Can not get Wall-Clock time!\n");
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE); 
			}
	}
	
}


double Island::getT_init(){
	return t_init;
}
	
double Island::getT_comp(){
	return t_comp;
}

double Island::getT_comm(){
	return t_comm;
}

double* Island::getT_compGA(){
	return t_compGA;
}

double Island::getT_compFitnessEval(){
	return t_compFitnessEval;
}

Island::~Island(){
}


