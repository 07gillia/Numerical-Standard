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


double x[3][3];
double v[3][3];
double mass[3];


void setUp() {
  x[0][0] = 0.4;
  x[0][1] = 0.0;
  x[0][2] = 0.0;

  x[1][0] = 2.0;
  x[1][1] = 0.0;
  x[1][2] = 0.0;

  x[2][0] = 3.0;
  x[2][1] = 0.0;
  x[2][2] = 0.0;

  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = 0.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;

  v[2][0] = 0.0;
  v[2][1] = 1.0;
  v[2][2] = 0.0;

  mass[0] = 0.2;
  mass[1] = 1.0;
  mass[2] = 0.01;
}


void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;

  for (int i=0; i<3; i++) {
    out << x[i][0]
        << ","
        << x[i][1]
        << ","
        << x[i][2]
        << std::endl;
  }
}



void updateBody() {
  double force[3];
  force[0] = 0.0;
  force[1] = 0.0;
  force[2] = 0.0;

  for (int i=0; i<2; i++) {
    const double distance = sqrt(
      (x[2][0]-x[i][0]) * (x[2][0]-x[i][0]) +
      (x[2][1]-x[i][1]) * (x[2][1]-x[i][1]) +
      (x[2][2]-x[i][2]) * (x[2][2]-x[i][2])
    );
    force[0] += (x[i][0]-x[2][0]) * mass[i]*mass[2] / distance / distance / distance ;
    force[1] += (x[i][1]-x[2][1]) * mass[i]*mass[2] / distance / distance / distance ;
    force[2] += (x[i][2]-x[2][2]) * mass[i]*mass[2] / distance / distance / distance ;
  }

  const double timeStepSize = 0.0001;

  x[2][0] = x[2][0] + timeStepSize * v[2][0];
  x[2][1] = x[2][1] + timeStepSize * v[2][1];
  x[2][2] = x[2][2] + timeStepSize * v[2][2];

  v[2][0] = v[2][0] + timeStepSize * force[0] / mass[2];
  v[2][1] = v[2][1] + timeStepSize * force[1] / mass[2];
  v[2][2] = v[2][2] + timeStepSize * force[2] / mass[2];
}


int main() {

  setUp();
  printCSVFile(0);

  const int timeSteps        = 2000000;
  const int plotEveryKthStep = 100;
  for (int i=0; i<timeSteps; i++) {
    updateBody();
    if (i%plotEveryKthStep==0) {
      printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}