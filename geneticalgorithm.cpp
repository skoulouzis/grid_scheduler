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
#include "geneticalgorithm.h"

GeneticAlgorithm::GeneticAlgorithm(vector<Schedule> aPopulation){
	population = aPopulation;
}

GeneticAlgorithm::GeneticAlgorithm(){

}

vector<Schedule> GeneticAlgorithm::parentsSelection(int kindOfSelection){
	vector<Schedule> parents;
	switch (kindOfSelection){
	case 0:
		{
			parents =  randomParentsSelection();
		} break;
		
	default:
		break;
		
	}
	return parents;
}


//Should it be alowed to select the same parent?
vector<Schedule> GeneticAlgorithm::randomParentsSelection(){
	vector<Schedule> parents;
	vector<Schedule>::iterator iterator;
	for(int i=0;i<2;i++){
		randomIndex = (random()%population.size())+0;
		parents.push_back(population[randomIndex]);
		iterator = population.begin() + randomIndex;
		population.erase(iterator);
	}
// 	cout << "randomSelection population.size() = "<< population.size() << endl;
// 	cout << "randomSelection parents[0].size() = "<< parents[0].size() << endl;
// 	cout << "randomSelection parents[1].size() = "<< parents[1].size() << endl;
	return parents;
}

vector<Schedule> GeneticAlgorithm::crossover(vector<Schedule> parents, int kindOfXover){
	vector<Schedule> offsprings;
// 	cout << "crossover parents[0].size() = "<< parents[0].size() << endl;
	switch (kindOfXover){
	case 0:
		{
			offsprings = onePointscrossover(parents);
		}
	default:
		break;
	}
// 	cout << "crossover offsprings.size() = "<< offsprings.size() << endl;
	return offsprings;
}

vector<Schedule> GeneticAlgorithm::onePointscrossover(vector<Schedule> parents){
	vector<Schedule> offsprings;
	Schedule parent1 = parents[0];
	Schedule parent2 = parents[1];
	Schedule offspring1,offspring2;
	
	randomIndex =  (random()%(parents[0].size()-1))+1;
	
// 	cout << "onePointscrossover randomIndex = "<< randomIndex << endl;
	
// 	cout << "onePointscrossover parent1.size() = "<< parent1.size() << endl;
// 	cout << "onePointscrossover parent2.size() = "<< parent2.size() << endl;
	for(int i=0;i<randomIndex;i++){
		offspring1.push_back(parent1.at(i));
		offspring2.push_back(parent2.at(i));
	}
	for(int i=randomIndex;i<parent1.size();i++){
		offspring1.push_back(parent2.at(i));
		offspring2.push_back(parent1.at(i));
	}
	
// 	cout << "parent1:" << endl;
// 	for(int i=0;i<parent1.size();i++){
// 		cout << "	"<<  parent1.at(i)<< endl;
// 		
// 	}
// 	cout << "parent2:" << endl;
// 	for(int i=0;i<parent1.size();i++){
// 		cout << "	"<<  parent2.at(i)<< endl;
// 		
// 	}
// 	cout << "offspring1:" << endl;
// 	for(int i=0;i<offspring1.size();i++){
// 		cout << "	"<<  offspring1.at(i)<< endl;
// 		
// 	}
// 	cout << "offspring2:" << endl;
// 	for(int i=0;i<offspring2.size();i++){
// 		cout << "	"<<  offspring2.at(i)<< endl;
// 		
// 	}
	
	offsprings.push_back(offspring1);
	offsprings.push_back(offspring2);
	return offsprings;
}

Schedule  GeneticAlgorithm::mutation(Schedule aSchedule, double rateOfMutation){
	
// 	cout << " fitness before=" << this->fitnessFunction(aSchedule,0) << endl;

	int size = aSchedule.size();
	randomIndex = round(size * rateOfMutation);	
	
	vector<int> randomNumbers;
	for(int i=0;i<randomIndex;i++){
		randomNumbers.push_back(i);
	}
	
	random_shuffle(randomNumbers.begin( ), randomNumbers.end( ));

	for(int i=0;i<randomNumbers.size();i++){
		aSchedule[randomNumbers.at(i)] =  random()%numOfResources+0;
	}
	
// 	cout << " fitness after=" << this->fitnessFunction(aSchedule,0) << endl;
	return aSchedule;
}


Schedule  GeneticAlgorithm::mutationWithHeuristic(Schedule aSchedule, double rateOfMutation){
	int size = aSchedule.size();
	randomIndex = round(size * rateOfMutation);	
	
// 	cout << "fitness before " << this->fitnessFunction(aSchedule,2) << endl;
	vector<int> randomNumbers;
	for(int i=0;i<randomIndex;i++){
		randomNumbers.push_back(i);
	}
	
	random_shuffle(randomNumbers.begin(), randomNumbers.end());

	for(int i=0;i<randomNumbers.size();i++){
		randomIndex = random()%numOfResources+0;
		while(this->resourcesLoad.at(randomIndex)!=0){
			randomIndex = random()%numOfResources+0;
		}

		aSchedule[randomNumbers.at(i)] = randomIndex;
	}
// 	cout << "fitness after " << this->fitnessFunction(aSchedule,2) << endl;
	
	return aSchedule;
}

double GeneticAlgorithm::fitnessFunction(Schedule aSchedule,int kinOfFitness){
	double executionTime;
	switch (kinOfFitness){
	case 0:
		{
			executionTime = executionTimeFunction(aSchedule);
		}
		break;
	case 1:
		{
			executionTime = flowTimeFunction(aSchedule);
		}
		break;
	case 2:
		{
			executionTime = getWaitTimeFunction(aSchedule);
		}
		break;
	default:
		break;
	}
	return executionTime;
}


/**This function returns the execution time only */
double GeneticAlgorithm::executionTimeFunction(Schedule aSchedule){
// 	executionTimes.clear();
	executionTime=0;
	mi=0;
	mips=0;

	for(int i=0;i<aSchedule.size();i++){
		mi = taskMiRatings.at(i);
		mips = resourceMipsRatings[aSchedule.at(i)];
		
		double temp = mi / mips;
		
		if(temp > executionTime) {
			executionTime = temp;
		}
	}
	
	return executionTime;
}

double GeneticAlgorithm::flowTimeFunction(Schedule aSchedule){
// 	flowTimes.clear();
	mi=0;
	mips=0;
	resLoad=0;
	double flowTime=0;
	for(int i=0;i<aSchedule.size();i++){
		mi = taskMiRatings.at(i);
		mips = resourceMipsRatings.at(aSchedule.at(i));
		resLoad = resourcesLoad.at(i);
// 		cout <<" flowTimeFunction" << resLoad << endl;
		double executionTime = mi / mips;
		double waitTime = resLoad / mips;
		
		//executionTime = executionTime * 0.2;
		
		if( (executionTime + waitTime) > flowTime ){
			flowTime = executionTime + waitTime;			
// 			if(aSchedule.at(i) == 7 && resLoad > 0){
// 				cout << executionTime << "+" << waitTime << "=" << flowTime << endl;
// 			}
		}
	}
	return flowTime;
}

/**this function returns the waittime only, for a schedule */
double GeneticAlgorithm::getWaitTimeFunction(Schedule aSchedule){
// 	flowTimes.clear();
	
	mips=0;
	resLoad=0;
	waitTime = 0;
	
	for(int i=0;i<aSchedule.size();i++){

		mips = resourceMipsRatings[aSchedule.at(i)];
		resLoad = resourcesLoad.at(i);
		
		double temp =  resLoad/mips;
		
// 		cout <<" getWaitTimeFunction" << resLoad << endl;
		if(temp > waitTime) {
			waitTime = temp;
		}

	}
// 	cout <<"getWaitTimeFunction: waitTime = " << waitTime << endl;
	return waitTime;
}


void GeneticAlgorithm::naturalSelection(vector<Schedule> candidates,int kindOfFitness){
	double fitness;
	double max;
	int index;	
	vector<Schedule>::iterator iterator;
		
// 	cout << "  ----------------------- " << endl;
	for(int i=0;i<2;i++){
		fitness=MAX_INT;
		max=MIN_INT;
// 		cout << "  ----------------------- " << endl;
		for(int j=0;j<candidates.size();j++){
			fitness = fitnessFunction(candidates.at(j),kindOfFitness);
// 			cout << " fitness " << fitness << " " << j << endl;
			if(fitness > max ){
				max = fitness;
				index = j;
			}
		}
		iterator = candidates.begin() + index;
		candidates.erase(iterator);
// 		cout << index << " max="<< max << " candidates.size "<< candidates.size() << endl;
	}
	
// 	cout << index << " max="<< max << " candidates.size "<< candidates.size() << endl;
	for(int i=0;i<candidates.size();i++){
		population.push_back(candidates.at(i));
	}
}

void rouletteWheelSelection(vector<Schedule> candidates,int kindOfFitness){
	double totalFitness;
}

double GeneticAlgorithm::getAveragefitness(int kindOfFitness){
	double averagefitness=0;
	vector<double> avfitness;
	for(int i=0;i<population.size();i++){
		averagefitness = averagefitness + (fitnessFunction(population.at(i),kindOfFitness));
		avfitness.push_back(fitnessFunction(population.at(i),kindOfFitness));
	}
	return averagefitness/population.size();
}

vector<Schedule> GeneticAlgorithm::getPopulation(){
	return this->population;
}

void GeneticAlgorithm::setPopulation(vector<Schedule> aPopulation){
	this->population = aPopulation;
}

void GeneticAlgorithm::setNumOfResources(int num){
	this->numOfResources = num;
}

void GeneticAlgorithm::setResourceMipsRatings(vector<double> mipsRatings){
	this->resourceMipsRatings = mipsRatings;
}

vector<double>  GeneticAlgorithm::getResourceMipsRatings(){
	return this->resourceMipsRatings;
}

void GeneticAlgorithm::setTaskMiRatings(vector<double> miRatings){
	this->taskMiRatings = miRatings;
}

vector<double>  GeneticAlgorithm::getTaskMiRatings(){
	return this->taskMiRatings;
}

void GeneticAlgorithm::setResourcesLoad(vector<double> resourcesLoad){
	this->resourcesLoad = resourcesLoad;
}

vector<double>  GeneticAlgorithm::getResourcesLoad(){
	return this->resourcesLoad;
}

vector<Schedule> GeneticAlgorithm::getBestSchedules(int n, int kinOfFitness){
	Schedule aSchedule;
	orderedPopulation.clear();
	vector<pair<double, int> > indexedFitness(population.size());	
	indexedFitness.clear();
	int index=0;
	
	for(int i=0;i<this->population.size();i++){
		aSchedule = population.at(i);
		indexedFitness.push_back(pair<double, int>(this->fitnessFunction(aSchedule,kinOfFitness), i));
	}
	
// 	cout << " before "<< endl;
// 	for (vector<pair<double, int> > ::iterator i = mirror.begin() ; i != mirror.end() ; ++i){
// 		cout << i->first << ", " << i->second << endl;
// 	}
	
	sort(indexedFitness.begin(), indexedFitness.end());
	
// 	cout << " GeneticAlgorithm: getBestSchedules will return "<< n << " schedules "<< endl;
	
	for (vector<pair<double, int> > ::iterator i = indexedFitness.begin() ; i != indexedFitness.end() ; ++i){
		index++;
		orderedPopulation.push_back( population.at(i->second) );
// 		cout << i->first << ", " << i->second << endl;
		if(index >= n) {
			break;
		}
	}	
	
	return orderedPopulation;
}

vector<int> GeneticAlgorithm::getIndexOfBestSchedules(int n,int kinOfFitness){
	Schedule aSchedule;
	vector<pair<double, int> > indexedFitness(population.size());	
	indexedFitness.clear();
	int index=0;
	vector<int> indexes;
	
	for(int i=0;i<this->population.size();i++){
		aSchedule = population.at(i);
		indexedFitness.push_back(pair<double, int>(this->fitnessFunction(aSchedule,kinOfFitness), i));
	}
	
	sort(indexedFitness.begin(), indexedFitness.end());
		
	for (vector<pair<double, int> > ::iterator i = indexedFitness.begin() ; i != indexedFitness.end() ; ++i){
		index++;
		indexes.push_back( i->second ); 
// 		cout << i->first << ", " << i->second << endl;
		if(index >= n) {
			break;
		}
	}	
	
	return indexes;
}


void GeneticAlgorithm::removeIndividualsByIndex(int index){
	vector<Schedule>::iterator iterator;
	iterator = population.begin() + index;
	population.erase(iterator);
}

void GeneticAlgorithm::removeIndividuals(int n){
	vector<Schedule>::iterator iterator;
	for(int i=0;i<n;i++){
		iterator = population.begin() + i;
		population.erase(iterator);
	}
}

void GeneticAlgorithm::addIndividual(Schedule aSchedule){
	this->population.push_back(aSchedule);
}

GeneticAlgorithm::~GeneticAlgorithm(){
}
