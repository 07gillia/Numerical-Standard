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



const int NumberOfParticles = 2;
// a variable to store the number of particles in the simulation
const double a = 0.1;
const double s = 0.1;
// set the constants for use in the force equations
double x[NumberOfParticles][3];
// the coordinates of each paricle
double v[NumberOfParticles][3];
// the velocity of each particle



void setUp(int NumberOfParticles) {
  srand (time(NULL));
  // set the seed of the random number generator
  /*
  for (int i = 0; i<NumberOfParticles; i++){
    x[i][0] = (double)rand()/(double)RAND_MAX;
    x[i][1] = (double)rand()/(double)RAND_MAX;
    x[i][2] = (double)rand()/(double)RAND_MAX;

    v[i][0] = 0.0;
    v[i][1] = 0.0;
    v[i][2] = 0.0;
  }
  */
  
  x[0][0] = 0.1;
  x[0][1] = 0.5;
  x[0][2] = 0.5;

  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  
  x[1][0] = 0.9;
  x[1][1] = 0.5;
  x[1][2] = 0.5;

  v[1][0] = 0.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;
  
}



void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "new_csv/result-" << counter <<  ".csv";
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



void updateBody(int NumberOfParticles) {
  for(int i = 0; i<NumberOfParticles; i++){
    //iterate through every particle

    double force[3];
    // a place to store the directional forces
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    for (int j=0; j<NumberOfParticles; j++) {
      // iterate through every particle that can effect the current particle
      if(i != j){
        // don't work out the force applied on itself

        double normalX = x[i][0]-x[j][0];
        double normalY = x[i][1]-x[j][1];
        double normalZ = x[i][2]-x[j][2];

        double distance = sqrt((normalX) * (normalX) + (normalY) * (normalY) + (normalZ) * (normalZ));
        // calculate the euclidian distance between two particles

        double euclidianForce = 4 * a * (((12 * pow(s,12)) / pow(distance,13)) - ((6 * pow(s,6))/pow(distance,7)));

        force[0] += (normalX)/distance * euclidianForce;
        force[1] += (normalY)/distance * euclidianForce;
        force[2] += (normalZ)/distance * euclidianForce;
        // calculate the force in each direction (x,y,z)

        //printf("X distance: %E, Total Distance: %E, Euclidian Force: %E, Total: %E\n", (x[i][0]-x[j][0]), distance, euclidianForce, (x[i][0]-x[j][0])/distance * euclidianForce);
        //printf("Y distance: %E, Total Distance: %E, Euclidian Force: %E, Total: %E\n", (x[i][1]-x[j][1]), distance, euclidianForce, (x[i][1]-x[j][1])/distance * euclidianForce);
        //printf("Z distance: %E, Total Distance: %E, Euclidian Force: %E, Total: %E\n", (x[i][2]-x[j][2]), distance, euclidianForce, (x[i][2]-x[j][2])/distance * euclidianForce);

        double otherX;
        double otherY;
        double otherZ;

        if(normalX>0){
          otherX = -1 + normalX;
        }
        else if(normalX<0){
          otherX = 1 + normalX;
        }
        else{
          otherX = 0;
        }
        // find the distance in the X direction
        if(normalY>0){
          otherY = -1 + normalY;
        }
        else if(normalY<0){
          otherY = 1 + normalY;
        }
        else{
          otherY = 0;
        }
        // find the distance in the X direction
        if(normalZ>0){
          otherZ = -1 + normalZ;
        }
        else if(normalZ<0){
          otherZ = 1 + normalZ;
        }
        else{
          otherZ = 0;
        }
        // find the distance in the X direction
        // find the distance in each axis in wrap around
        //printf("%E %E %E\n", otherX, otherY, otherZ);

        distance = sqrt(otherX * otherX + otherY * otherY + otherZ * otherZ);
        // calculate the euclidian distance between two particles in wrap around
        //printf("%E\n", distance);

        double otherEuclidianForce = 4 * a * (((12 * pow(s,12)) / pow(distance,13)) - ((6 * pow(s,6))/pow(distance,7)));
        // using the wrap around distance use the wrap around force

        force[0] += otherX/distance * otherEuclidianForce;
        force[1] += otherY/distance * otherEuclidianForce;
        force[2] += otherZ/distance * otherEuclidianForce;
        // split the forces into their x,y and z components
        
        //printf("Other X distance: %E, Total Distance: %E, Euclidian Force: %E, Total: %E\n", otherX, distance, otherEuclidianForce, otherX/distance * otherEuclidianForce);
        //printf("Other Y distance: %E, Total Distance: %E, Euclidian Force: %E, Total: %E\n", otherY, distance, otherEuclidianForce, otherY/distance * otherEuclidianForce);
        //printf("Other Z distance: %E, Total Distance: %E, Euclidian Force: %E, Total: %E\n", otherZ, distance, otherEuclidianForce, otherZ/distance * otherEuclidianForce);
        //printf("\n");
      }
    }
    const double timeStepSize = pow(10,-4);
    // set the timesteps

    v[i][0] = v[i][0] + timeStepSize * force[0];
    v[i][1] = v[i][1] + timeStepSize * force[1];
    v[i][2] = v[i][2] + timeStepSize * force[2];
    // calculate the velocity given the force and timestepsize

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
}



int main() {

  setUp(NumberOfParticles);
  printCSVFile(0);
  
  const int timeSteps        = 2000000;
  const int plotEveryKthStep = 100;
  for (int i=0; i<timeSteps; i++) {
    updateBody(NumberOfParticles);
    if (i%plotEveryKthStep==0) {
      printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}