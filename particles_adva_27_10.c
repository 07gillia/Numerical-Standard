// Compile with
//
// g++ -O3 particles.c -o particles
//
// Run it with
//
// ./particles
//
// Modelling a set number of particles in a unit cube using periodic bounding with forces applied using Lennard-Jones potential.
//
// Alex Gillies



#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <vector>



const int NumberOfParticles = 10;
// a variable to store the number of particles in the simulation

const double a = pow(10,-5);
const double s = pow(10,-5);
// set the constants for use in the force equations

double x[NumberOfParticles][3];
// the coordinates of each paricle

double v[NumberOfParticles][3];
// the velocity of each particle

double timeStepSize = pow(10,-4); // 0.0001
// set the timesteps
// the best timestep for a = s = 0.1 is pow(10,-4)
// the relative best timestep for a = s = pow(10,-5) is pow(10,-4) at distance 0.00001

double current_shortest = 1;
// variable for the timestep calculation

const double extra_dimentions[27][3]  = {{0,0,0},{0,0,1},{0,0,-1},{0,1,0},{0,1,1},{0,1,-1},{0,-1,0},{0,-1,1},{0,-1,-1},{1,0,0},{1,0,1},{1,0,-1},{1,1,0},{1,1,1},{1,1,-1},{1,-1,0},{1,-1,1},{1,-1,-1},{-1,0,0},{-1,0,1},{-1,0,-1},{-1,1,0},{-1,1,1},{-1,1,-1},{-1,-1,0},{-1,-1,1},{-1,-1,-1}};
// a static array that will store the change in dimentions for each of the wrap around particles
// the number of wrap around particles is 26 plus the center box {0,0,0}

std::vector<int> verletList[NumberOfParticles];
std::vector<int> verletListMirror[NumberOfParticles];
// the data structure to store the local particles and their mirror version

const int updateIndex = 10000;
// how often to update the verlet lists

double cutOff = 2.5 * a;
// the cutoff distance

double cutOffSqrt = cutOff * cutOff;
// a variable stored such that there is no need to use sqrt very single time



void setUp(int NumberOfParticles) {
  srand (time(NULL));
  // set the seed of the random number generator
  for (int i = 0; i<NumberOfParticles; i++){
    x[i][0] = (double)rand()/(double)RAND_MAX;
    x[i][1] = (double)rand()/(double)RAND_MAX;
    x[i][2] = (double)rand()/(double)RAND_MAX;

    v[i][0] = 0.0;
    v[i][1] = 0.0;
    v[i][2] = 0.0;
  }
}



double getTimeStepSize(double distance){
  double timeStepSize;
  if (distance > 0.00015){
    timeStepSize = pow(distance,2.7469)*pow(10,10);
  }
  else{
    timeStepSize = distance;
  }
  return timeStepSize;
}



void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "csv_3/" << NumberOfParticles << "/27/result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;

  for (int i=0; i<NumberOfParticles; i++) {
    out << x[i][0]
    << ","
    << x[i][1]
    << ","
    << x[i][2]
    << std::endl;
  }
}



void updateBody(int NumberOfParticles, int current_timestep) {
  current_shortest = 1000;
  // counter for use in calculateion of timesteps

  if (current_timestep % updateIndex == 0){
    // if the current timestep is an update timestep
    for(int i = 0; i<NumberOfParticles; i++){
    //iterate through every particle

      double force[3];
      // a place to store the directional forces
      force[0] = 0.0;
      force[1] = 0.0;
      force[2] = 0.0;

      verletList[i].clear();
      verletListMirror[i].clear();
      // clear the datastructures so they may be updated

      for (int j=0; j<NumberOfParticles; j++) {
        // iterate through all other particles

        if(i != j){
          // don't do anything to apply a force on itself
          for(int k = 0; k < 27; k++){
            // iterate through the mirrors

            double normalX = x[i][0]-(x[j][0] + extra_dimentions[k][0]);
            double normalY = x[i][1]-(x[j][1] + extra_dimentions[k][1]);
            double normalZ = x[i][2]-(x[j][2] + extra_dimentions[k][2]);
            // calculate the distance between the current particle and the mirror particle in the box 'k'

            double distance = (normalX) * (normalX) + (normalY) * (normalY) + (normalZ) * (normalZ);
            // calculate the non square rooted euclidean distance 

            if (distance < cutOffSqrt){
              // if the mirror particle is close enough
              verletList[i].push_back(j);
              verletListMirror[i].push_back(k);
              // store the particle and its translation index
            }
          }
        }
      }
    }
  }
  // this whole if statement is only called if the timestep is 

  for (int i = 0; i < NumberOfParticles; i++){
    double force[3];
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;

    for (int j = 0; j < verletList[i].size(); j++){
      // iterate through the length of list of nearby particles

      double normalX = x[i][0]-(x[verletList[i][j]][0] + extra_dimentions[verletListMirror[i][j]][0]);
      double normalY = x[i][1]-(x[verletList[i][j]][1] + extra_dimentions[verletListMirror[i][j]][1]);
      double normalZ = x[i][2]-(x[verletList[i][j]][2] + extra_dimentions[verletListMirror[i][j]][2]);
      // calculate the distance between the current particles and the jth item in the list from the mirror stored in the jth index of the mirror list

      double distance = sqrt((normalX) * (normalX) + (normalY) * (normalY) + (normalZ) * (normalZ));
      // calculate the euclidian distance between two particles

      if(distance < current_shortest){
        // keep track of the shortest distance between 
        current_shortest = distance;
      }

      double euclidianForce = 4 * a * (((12 * pow(s,12)) / pow(distance,13)) - ((6 * pow(s,6))/pow(distance,7)));
      // calculate the force between the two particles

      force[0] += (normalX)/distance * euclidianForce;
      force[1] += (normalY)/distance * euclidianForce;
      force[2] += (normalZ)/distance * euclidianForce;
      // calculate the force in each direction (x,y,z)

      v[i][0] = v[i][0] + timeStepSize * force[0];
      v[i][1] = v[i][1] + timeStepSize * force[1];
      v[i][2] = v[i][2] + timeStepSize * force[2];
      // calculate the velocity given the force and timestepsize
    }
  }
  // this is always done

  for (int i = 0; i < NumberOfParticles; i++){

    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];
    // calculate the new coordinates of the particle

    if(x[i][0]<=0){
      x[i][0] = 1;
    }
    if(x[i][0]>1){
      x[i][0] = 0;
    }
    if(x[i][1]<=0){
      x[i][1] = 1;
    }
    if(x[i][1]>1){
      x[i][1] = 0;
    }
    if(x[i][2]<=0){
      x[i][2] = 1;
    }
    if(x[i][2]>1){
      x[i][2] = 0;
    }    
    // anything that goes below 0 is wrapped around and anything above 1 is wrapped around
  }
  //printf("%f\n", current_shortest);
  timeStepSize = getTimeStepSize(current_shortest);
}



int main() {

  clock_t begin = clock();

  setUp(NumberOfParticles);
  //printCSVFile(0);
  
  const int timeSteps        = 200000;
  const int plotEveryKthStep = 100;
  for (int i=0; i<timeSteps; i++) {

    printf("%d\r", (i * 100/timeSteps));
    // percentage of how far through we are

    updateBody(NumberOfParticles, i);
    if (i%plotEveryKthStep==0) {
      //printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
    }
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::stringstream filename;
  filename << "csv_3/" << NumberOfParticles << "/27/result.txt";
  std::ofstream out( filename.str().c_str() );

  out << elapsed_secs << std::endl;

  printf("%d\n", NumberOfParticles);
  printf("%f\n", elapsed_secs);

  return 0;
}
