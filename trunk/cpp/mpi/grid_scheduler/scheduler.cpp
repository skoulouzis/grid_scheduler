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

// cluster/workspace/SESP/debug/src/notThis; mpirun.globus -n 2 /cluster/workspace/SESP/debug/src/sesp -r 100 -p 200 -b 6

#include "mpi.h"
#include "definitions.h"
#include "island.h"
#include "grid.h"
#include <sys/time.h>


static int parseargs(int argc, char *argv[]);
int numResources;
int batchSize;
int populationSize,localPopulationSize;
int popInitMethod;//the way the population is initialised
int fitnessType;//the fitness type used to evaluate candidates
int maxTasks;
int numTasks;
double startwtime = 0.0;
double endwtime=0.0;
double total_time=0.0;
// int arguments[6];
Grid *g;
Island *island;
ResourceList gridResourceList;
//MPI stuff
int numprocs, myrank;
//MPI data type describing the resource
MPI_Datatype Gresource_type;
/*The arrays to fill and pass as arguments to MPI_Type_struct */
int blocklens[5];
MPI_Aint displs[5];
MPI_Datatype oldtypes[5];
MPI_Status status;

Gresource ress[99999];
int *schdulesIn;
int *schdulesOut;
Schedule mySchedule;
vector<Schedule> returnedSchedules;
static int nprocsX;
static int nprocsY;
int dims[2], periods[2], coords[2], neighbour[2];
static MPI_Comm ring; 
static int myLeft;         /* rank of left mpi process in grid/line topology       */
static int myRight;        /* rank of right mpi process in grid/line topology      */
static int myUp;		/* rank of up mpi process in grid/line topology      */
static int myDown;	/* rank of down mpi process in grid/line topology      */

/*Performance measures*/
static int reporter;
static double t_init=0.0;
static double t_comp=0.0;
static double t_compRunGA=0.0;
static double t_comm=0.0;
static	double td=0.0;
static	long sec,usec;
static	struct timeval tp;


/**
Simulates the retrieval of a batch of jobs from a pool of jobs
submitted by users. In fact the pool does not exist and instead 
jobs are created on the fly.
*/
vector<double> getBatch(int batchSize) {
	
	vector<double> batch;
	double mi;
	for(int i = 0; i < batchSize; i++) {
		mi = double(random()) / double(RAND_MAX) * 500;
		batch.push_back(mi);
		numTasks++;
	}

// 	cout << "Main @" << myrank << ":\n  -created batch with " << batch.size() << " tasks " << " total tasks: " << numTasks << endl;
	
	return batch;
}

void refreshGridResourceList(){
	if(myrank==0){
		gridResourceList.clear();
		gridResourceList = g->getResourceList();
	}
}

//why? - for output stupid!!!
int getGaIterations(){
	int iters = island->getIterations();
	return iters;
}


void initResourceStruct(){
	/*Putting each element of struct params into the new datatype... */
	//resource ID
	blocklens[0] = 1;
	displs[0] = 0;
	oldtypes[0] = MPI_INT;
	
	//resource mips
	blocklens[1] = 1;
	displs[1] =  sizeof(int);
	oldtypes[1] = MPI_DOUBLE;
	
	//resource requests	
	blocklens[2] = 1;
	displs[2] = sizeof(int) + sizeof(double);
	oldtypes[2] = MPI_DOUBLE;
	
	//resource load
	blocklens[3] = 1;
	displs[3] = sizeof(int) + sizeof(double) + sizeof(double);
	oldtypes[3] = MPI_DOUBLE;
	
	//reliability
	blocklens[4] = 1;
	displs[4] = sizeof(int) + sizeof(double)+ sizeof(double)+ sizeof(double);
	oldtypes[4] = MPI_DOUBLE;
	
	/*finally, create the datatype and commit changes!*/
	MPI_Type_struct( 5, blocklens, displs, oldtypes, &Gresource_type );
	MPI_Type_commit( &Gresource_type );
}

static void timeStart(){

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

static void timeEnd(int what){
	if(myrank==0)
	{
		if(gettimeofday(&tp,NULL)==0){
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
					t_compRunGA=t_compRunGA+td;
				break;
				default :
					t_init =t_init + td;
				break;
			}
			
		}else{
			perror("Can not get Wall-Clock time!\n");
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE); 
			}
	}
	
}

void packAndBroadcastResourceList(){
// 	ress = new Gresource[numResources];
	if(myrank==0){
		for(int i=0;i<numResources;i++){
			ress[i].id =  gridResourceList.at(i).id;
			ress[i].mips = gridResourceList.at(i).mips;
			ress[i].requests = gridResourceList.at(i).requests;
			ress[i].load = gridResourceList.at(i).load;
			ress[i].reliability =  gridResourceList.at(i).reliability;
		}
// 		cout <<"rank:"  << myrank << " ress[0].mips = " << ress[0].mips << endl;
// 		MPI_Send(ress1, numResources , Gresource_type, 1, 123, grid);	
	}
	MPI_Bcast( &ress, numResources, Gresource_type, 0,MPI_COMM_WORLD );
}

void unpackResourceList(){
	if(myrank!=0){
	gridResourceList.clear();
		for(int i=0;i<numResources;i++){
			gridResourceList.push_back(ress[i]);
		}
	}
// 	cout <<"rank:"  << myrank << " gridResourceList.at(0).mips  = " << gridResourceList.at(0).mips  << endl;
}

void gatherSchdules(){
	//allocat memory for the send buffer, and unpack schdules returned from the GA
	schdulesOut = (int *)malloc( (mySchedule.size()) * sizeof(int));
	for(int i=0;i<mySchedule.size();i++){
		schdulesOut[i] = mySchedule.at(i);
	}
	if(myrank==0){
		schdulesIn = (int *)malloc( (numprocs * mySchedule.size()) * sizeof(int));
		returnedSchedules.clear();
	}
 	MPI_Gather(schdulesOut, mySchedule.size(), MPI_INT, schdulesIn, mySchedule.size() , MPI_INT, 0, MPI_COMM_WORLD);
	
	if(myrank==0){
		int counter=0;
		Schedule tempSchedule(batchSize);
		for(int i=0; i< numprocs * batchSize;i++){
			if(counter == (batchSize) ){			
				counter = 0;
				returnedSchedules.push_back(tempSchedule);
			}
			tempSchedule.at(counter) = schdulesIn[i];
			counter++;
			//add the last schedule
			if(i == (numprocs * batchSize) - 1){
				returnedSchedules.push_back(tempSchedule);
			}
		}
// 		cout <<  endl;
// 		for(int i=0;i<returnedSchedules.size();i++){
// 			tempSchedule.clear();
// 			tempSchedule = returnedSchedules.at(i);
// 			cout << " schedule vector " <<  i << endl;
// 			for(int j=0;j<tempSchedule.size();j++){
// 				cout << " 	" << tempSchedule.at(j) << endl;
// 			}
// 		}
	}	
	free(schdulesOut);
	if(myrank==0){
		free(schdulesIn);
	}
}

Schedule getBestSchedules(int kinOfFitness){
	if(myrank==0){
		double fitness=MAX_INT;
		int index;
		Schedule tempSchedule;
		for(int i=0;i<returnedSchedules.size();i++){
			tempSchedule.clear();
			tempSchedule = returnedSchedules.at(i);
// 			fitness = island->getFitness(tempSchedule,kinOfFitness);
			if(island->getFitness(tempSchedule,kinOfFitness) < fitness){
				fitness = island->getFitness(tempSchedule,kinOfFitness);
				index = i;
			}
		}
// 		cout << "best fitness " << fitness << endl;
		return returnedSchedules.at(index);
	}
}

void runSimulation(){
	numTasks = 0;
	if(myrank==0){
// 		cout << "Main: \n  -starting simulation with " << numprocs << " processors " << " local population is set to "<< localPopulationSize<<endl;
	}

	//first do some set up stuff
	
// 	srand(static_cast<unsigned int>(clock()) * myrank);
	
	//rank 0 creates the grid and broadcasts the resource list to the other
	//nodes
	if(myrank==0){
		//create a new grid object. 
		g= new Grid(numResources);
		//get the resource list
		gridResourceList = g->getResourceList();
	}

	packAndBroadcastResourceList();
	
	unpackResourceList();


	//create an island passing some constants for use later
	island = new Island(localPopulationSize, popInitMethod, fitnessType,myrank);

	island->setCommunicator(ring);	
		
	island->setNeighbors(myLeft,myRight);
// 	island->setMyRank(myrank);
			
	//then start the simulation n tings
	//get time
// 	if(gettimeofday(&tp,NULL)==0 && myrank==0){
// 		sec=tp.tv_sec; usec=tp.tv_usec;
// 	}
	MPI_Barrier(ring);
	
	
	//Start Timing
	if (myrank == 0){
		startwtime = MPI_Wtime();
	}
		
	//the main control loop
	while (numTasks < maxTasks) {
	
		timeStart();	
		//get a batch of jobs
		vector<double> batch = getBatch(batchSize); //all			
		//refresh resource list
		refreshGridResourceList(); // rank 0		
		//get resource list and run scheduler 
				
		mySchedule.clear();
		mySchedule = island->runGA(gridResourceList, batch);//all
		timeEnd(2);
		
		timeStart();
		//recive schdules and submit schedule to the grid 
		gatherSchdules();
		timeEnd(0);

		timeStart();
		if(myrank==0){
			g->executeSchedule(getBestSchedules(fitnessType), batch);//rank 0
// 			update resource load data (a bit like a rtimestep)
			g->updateResourceLoadData(getGaIterations()); // rank 0
// // 			write resource data to file
			g->appendResourceUtilData(); //rank 0
		}
		timeEnd(1);
		
		timeStart();
		packAndBroadcastResourceList();
		timeEnd(0);
		
		timeStart();
		unpackResourceList();	
		timeEnd(1);
	}

	//get time again and calculate execution time
	if(myrank==0){
		if(gettimeofday(&tp,NULL)==0){

// 			sec=tp.tv_sec-sec; usec=tp.tv_usec-usec;
// 
// 			td=sec+usec/1.0e6;
			endwtime = MPI_Wtime();
			total_time = endwtime - startwtime;
			
			cout << numprocs << " " << total_time << endl;

// 			cout << "Sceduler time measures"<< endl;
// 			cout << numprocs <<"	total_time:" << total_time << " t_init:"<< t_init << " t_comp:"<< t_comp << " t_comm:" << t_comm << " run ga:"<< t_compRunGA <<" sum:"<< (t_init+t_comp+t_comm+t_compRunGA) <<endl;
// 			cout << "Island time measures"<< endl;
// 			cout << " 	T_init:" << island->getT_init() << " T_comp:"<< island->getT_comp() <<" T_comm:" <<island->getT_comm() << " sum:"<<(island->getT_init() + island->getT_comp() + island->getT_comm())<< endl;
// 			cout << "GATimes"<< endl;
// 			cout << " 	selection:" << island->getT_compGA()[0] << " Xover:"<< island->getT_compGA()[1] <<" mutation:" <<island->getT_compGA()[2] << "  natural Selection:" << island->getT_compGA()[3]  << " fitness eval:"<< island->getT_compFitnessEval()<<" sum:"<< (island->getT_compGA()[0]+island->getT_compGA()[1]+island->getT_compGA()[2]+island->getT_compGA()[3]+island->getT_compFitnessEval()) << endl;
			
// 			cout << " Iterartons:" << getGaIterations() << " LocalPop:"<< localPopulationSize <<endl;
		}
		g->~Grid();
	}
}


static int parseargs(int argc, char *argv[]) {
	
	int check_options = 0;
	char opt;
	
	/* Check options */
	while (EOF != (opt =  getopt(argc, argv, "r:t:p:b:i:f:") )) {  
		switch (opt) {
			
		case 'r':
			//get argument from argv, I think
			numResources = atoi(optarg);
			//check that it's a valid argument
			if(numResources > 1) check_options += 1;
// 			cout << "numResources set to: " << numResources << endl;
// 			arguments[0] = numResources;
			break;
		case 't':
			//get argument from argv, I think
			maxTasks = atoi(optarg);
			//check that it's a valid argument
			if(maxTasks > 1) check_options += 2;
// 			cout << "maxTasks set to: " << maxTasks << endl;
// 			arguments[1] = maxTasks;
			break;
		case 'p':
			//get argument from argv, I think
			populationSize = atoi(optarg);
			//check that it's a valid argument
			if(populationSize > 1) check_options += 4;
// 			cout << "populationSize set to: " << populationSize << endl;
// 			arguments[2] = populationSize;
			break;
			
		case 'b':
			//get argument from argv, I think
			batchSize = atoi(optarg);
			//check that it's a valid argument
			if(batchSize > 1) check_options += 8;
// 			cout << "batchSize set to: " << batchSize << endl;
// 			arguments[3] = batchSize;
			break;
		case 'i':
			//get argument from argv, I think
			popInitMethod = atoi(optarg);
			//check that it's a valid argument
			check_options += 15;
// 			cout << "popInitMeth set to: " << popInitMethod << endl;
// 			arguments[4] =popInitMethod;
			break;
		case 'f':
			//get argument from argv, I think
			fitnessType = atoi(optarg);
			//check that it's a valid argument
			check_options += 30;
// 			cout << "fitness type set to: " << fitnessType << endl;
// 			arguments[5] = fitnessType;
			break;
			
		default:
// 			cout << "Usage: -r <num_resources>" << endl;
			exit(EXIT_FAILURE);
			break;
		}
	}
	
    	/* Check if all needed options have been set */
	if (check_options != 60) {
		//If required arguments are not supplied print usage and quit
// 		cout << "Usage: -r <num_resources> -p <populationSize> -b <batchSize> -i <initialisation(int)> -f <fitnessType(int)>" << endl;
    		//exit(EXIT_FAILURE);
	}
	return 0;
}//end parseargs

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	if(myrank==0){
		if(parseargs(argc, argv)) { 
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
			exit(EXIT_FAILURE); 
		}				
	}
	
	// Brodcast arguments to the nodes
	MPI_Bcast ( &populationSize, 1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast ( &batchSize, 1, MPI_INT, 0,MPI_COMM_WORLD );
	MPI_Bcast ( &numResources, 1, MPI_INT, 0,MPI_COMM_WORLD );
	MPI_Bcast ( &popInitMethod, 1, MPI_INT, 0,MPI_COMM_WORLD );
	MPI_Bcast ( &fitnessType, 1, MPI_INT, 0,MPI_COMM_WORLD );
	MPI_Bcast ( &maxTasks, 1, MPI_INT, 0,MPI_COMM_WORLD );

	
// 	if (myrank==1){
// 		localPopulationSize = (populationSize/numprocs) + (populationSize%numprocs) ;
// 	}
// 	else{
// 		localPopulationSize = populationSize/numprocs;
// 	}
	if (myrank<populationSize%numprocs){
		localPopulationSize = (populationSize/numprocs) +1;
	}else{
		localPopulationSize = populationSize/numprocs;
	}
		
	initResourceStruct();
		
	dims[0] = numprocs;

	//0 for line topology, 1 for ring topology
	periods[0] = 1; 


	//cretae a ring totpology to impruve prrformnce 
	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &ring);
	
	//get new rank
	MPI_Comm_rank(ring, &myrank);

	
	//get coordinates 
	MPI_Cart_get(ring,1,dims,periods,coords);

	
	int lcoord[dims[0]];
	int rcoord[dims[0]];

	lcoord[0] = coords[0]-1;
	rcoord[0] = coords[0]+1;

	//left right nodes
	MPI_Cart_rank(ring,lcoord,&myLeft);
	MPI_Cart_rank(ring,rcoord,&myRight);

	
// 	myRight = (myrank + 1) % numprocs;
// 	myLeft = myrank - 1;
// 	if (myLeft < 0)
// 		myLeft = numprocs - 1;
	
	srandom(31416+10*(myrank));

// 	ring = MPI_COMM_WORLD;
		
	MPI_Barrier(ring);
	
	runSimulation();
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}

