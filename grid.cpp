#include "grid.h"

/**
constructor initialises grid environment 
creates numResources resources and sets their various attributes
*/
Grid::Grid(double numResources) {
	//temp resource object
	Gresource res;
	//create output file
	output.open("grideData.txt");
	char fileName[30];
	totalMips = 0;
	totalLoad = 0;
	resourceUtiliziation = 0; 
	ios_base::fmtflags flags = output.flags();
	width = 20;
// 	string tab= "	";
	
	output << left                      // Left-justify in each field
		<< "," << "#Tota Load"      // Then, repeatedly set the width
		<< "," << "Total Mips "     // and write some data
		<< "," << "Total Load/Mips"
		<< "," << "Resource Utilization"

		<< endl;
// 	output.flags(flags);
	//create numResource resources
	for(int i=0; i<numResources; i++) {
		sprintf(fileName,"R%d_Data.txt",i);
		resourceOutput[i].open(fileName);
		resourceOutput[i] << left
			<< "," << "#LOAD"
			<< "," << "Mips"
			<< "," << "Load/Mips"
			<< "," << "Requests"
			<< "," << "Mis @ this step"
			<< "," << "MIs Executed @ this step"
			<< endl;
		
		res.id = i;
// 		res.mips = ((rand()%10)+1) * 1.2; 
		res.mips = double(rand()) / double(RAND_MAX) * 10; 
		
		totalMips = totalMips + res.mips;
		res.load = 0;//temporary value 

		totalLoad = totalLoad + res.load;
		res.reliability = double(rand()) / double(RAND_MAX);//temporary value
		res.requests = 0;
		//ResourceList.push_back(res);
		gridResourceList.push_back(res);
// 		cout << res.id << " with mips: <<" << res.mips << " requests: "<< res.requests<< endl;
		
		if(res.load>0){
			resourceUtiliziation++;
		}
	}
	
	appendResourceUtilData();
	
}//end Grid constructor


Grid::~Grid(){
	output.close();
	for(int i=0;i<this->gridResourceList.size();i++){
		resourceOutput[i].close();
	}
}


/**
returns ResourceList to scheduler
*/
ResourceList Grid::getResourceList() {
	
	return gridResourceList;
}//end getResourceList


/**
return resource id of resource at index in vector	
*/
int Grid::getId(int index) {
	
	Gresource res = gridResourceList.at(index);//ResurceList[index];
	
	return res.id;
}//end getID


/**
returns mips rating of rsource with id
*/
double Grid::getComputationalPotential(int id) {
	
	Gresource res = gridResourceList.at(id);//ResurceList[id];
	
	return res.mips;
}//end getComputationalPotential


/**
recieves schedule from scheduler and attempts 
to allocate tasks to resources specified therein...

returns int vector representing task allocation results
*/
void Grid::executeSchedule(Schedule schedule, vector<double> miReq) {
	
	
// 	cout << "GRID:" << "  -Recieved schedule with " << schedule.size() << " tasks and " << "miReq vector with " << miReq.size() << " values" << endl;
	
	//vector<int> results;
	
	//iterate through vectors
	for(int i= 0; i < schedule.size(); i++) {
		
		//get resourceID for each task in schedule
		int resourceID = schedule.at(i);
		
		//get the mi requirement of each task in schedule
		double taskMI =  miReq.at(i);
		
		//attempt to submit each task, storing result in result vector		
		
		//results.push_back();
		Gresource res = gridResourceList.at(resourceID);//ResourceList[resourceID];
		res.requests++;
	

		res.load = res.load + taskMI;
		thisMI[resourceID] = taskMI;
		//return updated Gresource object to resourceList 
		gridResourceList.at(resourceID) = res;
	}
// 	cout << "\n  -done submitting tasks\n" << endl;
// 	appendResourceUtilData();//^
	//return results;
}//end executeSchedule


/**
get the current load of resource with id
*/
double Grid::getLoad(int id) {
	
	Gresource res = gridResourceList.at(id);//ResurceList[id];
	
	return res.load;
}//end getLoad


/**
get the current availibility of resource with id
*/
bool Grid::getReliability(int id) {
	
	Gresource res = gridResourceList.at(id);//ResourceList[id];
	
	return res.reliability;
	
}//end getAvailable


/**
Updates resource load data. Simulates the passing of one second, at the moment
by subtracting mips rating from load
*/
void Grid::updateResourceLoadData(int iteartions){
	for(int i = 0; i < gridResourceList.size(); i++) {
	
		//get each resource and update it's load if >0;
		Gresource res = gridResourceList.at(i);
		
		//if res.load greater than zero, update load
// 		if(res.load > 0) { ;
		thisMIExecuted[i] = (res.mips * (iteartions/250.0));//^
		res.load = res.load - thisMIExecuted[i];
// 			if new load is less than zero set it to zero
			if(res.load < 0) {
				res.load = 0;	
			}
			
// 		}
		gridResourceList.at(i) = res;//^ put the resource back
	}
// 	appendResourceUtilData();//^
}


/**
set the availibility of  resource with id
*/
void Grid::setReliability(int id, double reliability) {
	
	Gresource res = gridResourceList.at(id);//ResurceList[id];
	
	res.reliability = reliability;
	
}//end setAvailable

/**
Appends a line of data to the grid output file detailing resource utilzation
load, mips, waitime, reliability, requests,
*/
void Grid::appendResourceUtilData(){
// 	if (output.is_open()) {
// 		
// 		output << " " <<"\n";
// 	}else cout << "Unable to open file";
	
	resourceUtiliziation = 0;
	Gresource res;
	for(int i=0;i<this->gridResourceList.size();i++){
		res = gridResourceList.at(i);
		resourceOutput[i] << left                      // Left-justify in each field
			<< "," << getResourceList().at(i).load       // and write some data
			<< "," << getResourceList().at(i).mips
			<< "," << (getResourceList().at(i).load/getResourceList().at(i).mips)
			<< "," << getResourceList().at(i).requests
			<< "," << thisMI[i] 
			<< "," << thisMIExecuted[i]
			<< endl;
		totalLoad = totalLoad + getResourceList().at(i).load;
		if(res.load>0){
			resourceUtiliziation++;
		}
	}
	output << left                      // Left-justify in each field
		<< "," << totalLoad      // Then, repeatedly set the width
		<< "," << totalMips       // and write some data
		<< "," << (totalLoad/totalMips)
		<< "," << (resourceUtiliziation/gridResourceList.size())*100.0
		<< endl;
	
}
