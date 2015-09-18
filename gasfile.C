#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include <fstream>

#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "ComponentConstant.hh"
#include "MediumMagboltz.hh"
#include "TrackHeed.hh"
#include "SolidBox.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  //GARFIELD++	
  Garfield::MediumMagboltz* mediumMagboltz = new Garfield::MediumMagboltz();
  mediumMagboltz->SetComposition("cf4", 100.);
  mediumMagboltz->SetTemperature(293.15);
  mediumMagboltz->SetPressure(30.);
  mediumMagboltz->SetFieldGrid(400,500,2,false);
  mediumMagboltz->GenerateGasTable(10,true);//# of collisions (x1e7), verbose
  mediumMagboltz->WriteGasFile("CF4.gas");

}
