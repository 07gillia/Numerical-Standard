// Translate this file with
//
// g++ -O3 spaceboddies.c -o spaceboddies
//
// Run it with
//
// ./spaceboddies
//
// Open Paraview (www.paraview.org) and do the following:
// - Select File/open and select all the results files. Press the Apply button.
// - Click into the left visualisation window (usually titled Layout #1).
// - Click the result-* item in the window Pipeline Browser. Ensure that your Layout #1 and the item result-* is marked.
// - Select Filters/Alphabetical/TableToPoints. Press Apply button.
// - Switch the representation (on top) from Surface into Points.
// - Press the play button and adopt colours and point sizes.
// - For some Paraview versions, you have to mark your TableToPoints item (usually called TableToPoints1) and explicitly select that X Column is x, Y Column is y, and Z Column is z.
// - What is pretty cool is the Filter TemporalParticlesToPathlines. If you set Mask Points to 1, you see a part of the trajactory.
//
// (C) 2015 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <cmath>

const int N = 2;
// set the number of points
double x[N][3];
// set the size of the array to store the coordinates of the points
double v[N][3];
// set the size of the array to store the coordinate of the points
double dimentionsArray[27][3] = {{0,0,0},{0,1,0},{0,-1,0},{-1,0,0},{1,0,0},{0,0,-1},{0,0,1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{1,0,-1},{-1,0,-1},{1,0,1},{-1,0,1},{0,1,-1},{0,-1,-1},{0,1,1},{0,-1,1},{-1,-1,1},{1,-1,1},{-1,-1,-1},{1,-1,-1},{-1,1,1},{1,1,1},{-1,1,-1},{1,1,-1}};


void setUp(int NumberOfParticles) {
  srand (10);
  // seed the random number generator
  //for (int s=0; s < NumberOfParticles; s++)
  // iterate through all the particles
  {
    x[0][0] = 0.3; //(double)rand()/(double)RAND_MAX;
    x[0][1] = 0.5; //(double)rand()/(double)RAND_MAX;
    x[0][2] = 0.5; //(double)rand()/(double)RAND_MAX;
    // create random X,Y,Z coordinates for each point

    
    x[1][0] = 0.7; //(double)rand()/(double)RAND_MAX;
    x[1][1] = 0.5; //(double)rand()/(double)RAND_MAX;
    x[1][2] = 0.5; //(double)rand()/(double)RAND_MAX;
    // create random X,Y,Z coordinates for each point
    

    v[0][0] = 0.0;
    v[0][1] = 0.0;
    v[0][2] = 0.0;
    // set the inital velocities

    
    v[1][0] = 0.0;
    v[1][1] = 0.0;
    v[1][2] = 0.0;
    // set the inital velocities
      
  }
}


void printCSVFile(int counter) {
  // output the csv files
  std::stringstream filename;
  filename << "csv/result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;

  for (int i=0; i<N; i++) {
    out << x[i][0]
    << ","
    << x[i][1]
    << ","
    << x[i][2]
    << std::endl;
  }
}


double useFormula(double distance)
{
  double a = pow(10,-5);
  double s = pow(10,-5);
  return (4 * a * (((12 * pow(s,12))/pow(distance,13)) - ((6 * pow(s,6)) / pow(distance,7))));
}


void updateBody(int NumberOfParticles) {

  for(int r=0; r<NumberOfParticles; r++){
  // iterate through all the particles

    double force[3];
    // set the size of the array to store the forces
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;
    // set the forces for each particle

    for (int i=0; i<NumberOfParticles; i++) {
      // iterate through each of the particles
      if(i != r){
        // this ensures that the particles doesn't calculate the force it applies on itself
        for (int axis=0; i<27;i++){
          // this is for all possible axis of the simulation
        
        const double distance = sqrt(
          (x[r][0]-(x[i][0])+dimentionsArray[axis][0]) * (x[r][0]-(x[i][0]+dimentionsArray[axis][0])) +
          (x[r][1]-(x[i][1])+dimentionsArray[axis][1]) * (x[r][1]-(x[i][1]+dimentionsArray[axis][1])) +
          (x[r][2]-(x[i][2])+dimentionsArray[axis][2]) * (x[r][2]-(x[i][2]+dimentionsArray[axis][2]))
          ); // find the distance between the two points using pythag

        force[0] += ((x[r][0]-x[i][0])/distance) * useFormula(distance);
        force[1] += ((x[r][1]-x[i][1])/distance) * useFormula(distance);
        force[2] += ((x[r][2]-x[i][2])/distance) * useFormula(distance);
        //calculate the force in each direction
      }
    }
  }

    const double timeStepSize = pow(10,20);
    // the length of the timestep

    v[r][0] = v[r][0] + timeStepSize * force[0];
    v[r][1] = v[r][1] + timeStepSize * force[1];
    v[r][2] = v[r][2] + timeStepSize * force[2];
    // calculates the velocity of each point

    x[r][0] = (x[r][0] + timeStepSize * v[r][0]);
    x[r][1] = (x[r][1] + timeStepSize * v[r][1]);
    x[r][2] = (x[r][2] + timeStepSize * v[r][2]);
    // calculating the coordinates of the current point
  }
}


int main() {

  setUp(N);
  printCSVFile(0);

  const int timeSteps        = 200000;
  const int plotEveryKthStep = 100;
  for (int i=0; i<timeSteps; i++) {
    updateBody(N);
    if (i%plotEveryKthStep==0) {
      printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}