#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1D.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include <fstream>

#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "ComponentConstant.hh"
#include "ComponentAnalyticField.hh"
#include "MediumMagboltz.hh"
#include "TrackHeed.hh"
#include "SolidBox.hh"
#include "AvalancheMC.hh"
#include "ViewSignal.hh"
#include "ViewField.hh"
#include "AvalancheMicroscopic.hh"
#include "ViewGeometry.hh"
#include "ViewCell.hh"
#include "DriftLineRKF.hh"

using namespace Garfield;

// Finds ion creation time based on distance from the electron cluster 
// to the wire (found by fitting the ion creation time of many avalanches)
double IonTiming(double dist) {

  double p0 =  1.49880e-13;   
  double p1 =  2.09250e+02;
  double p2 =  2.61998e+02;
  double p3 = -1.24766e+02;
  
  return p0 + p1 * dist + p2 * dist * dist + p3 * dist * dist * dist;

}

double transfer(double t) {

  const double tau =25;//ns
  return (t/tau)*std::exp(1-t/tau);

}

int main(int argc, char * argv[]) {

//  double dedx = 6.69e-3;//ion energy loss rate in gas MeV/mm (10torr C3F8)
  double dedx = 1.35e-2;//ion energy loss rate in gas MeV/mm (20torr C3F8)
//  double dedx = 2.03e-2;//ion energy loss rate in gas MeV/mm (30torr C3F8)
  double W = 25.;//gas work function eV
  const double Anode_Voltage = 800.;//V
  const double Frame_Space = 0.4;//cm
  const double tMax = 150;//ns

  TApplication app("app", &argc, argv);

  Garfield::MediumMagboltz* gas = new Garfield::MediumMagboltz();
//  gas->SetComposition("c3f8", 100.);
//  gas->SetTemperature(293.15);//K
//  gas->SetPressure(20.);//Torr
//  gas->SetFieldGrid(100,2000,2,false);
//  gas->GenerateGasTable(10,true);//# of collisions (x1e7), verbose
//  gas->SetMaxElectronEnergy(6000);
//  gas->WriteGasFile("C3F8_20torr.gas");
  gas->LoadGasFile("C3F8_20torr.gas");
  gas->LoadIonMobility("IonMobility_Ar+_Ar.txt");
//  gas->Initialise(false);

  Garfield::ComponentAnalyticField* cmp = new Garfield::ComponentAnalyticField();

  const double Electrode_Len = 10.0;//cm

  //Gas volume of ionisation chamber has the shape of a cube with 1 cm edge length in x,y,z  
  Garfield::SolidBox* solidBox = new Garfield::SolidBox(0,Frame_Space/2.,0,Electrode_Len/2.,Frame_Space/2.,Electrode_Len/2.);
  Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
  geo->AddSolid(solidBox,gas);

  cmp->SetGeometry(geo);

  cmp->AddPlaneY(Frame_Space,Anode_Voltage,"a0");//anode plane
  cmp->AddPlaneY(0,0,"c0");//cathode plane

/*
  ViewGeometry* view = new ViewGeometry();
  view->SetGeometry(geo);
  view->Plot();
*/
/*
  ViewCell* view = new ViewCell();
  view->SetComponent(cmp);
  view->Plot2d();
*/
  cmp->AddReadout("a0");
  
  //The sensor links the detector description and the particle transport
  //In the example here: 
  //Detector description: material (MediumMagboltz), geometry (GeometrySimple), fields (ComponentConstant) 

  Garfield::Sensor* sensor = new Garfield::Sensor();
  sensor->AddComponent(cmp);
  sensor->AddElectrode(cmp,"a0");
  sensor->AddElectrode(cmp,"c0");
  sensor->SetTransferFunction(transfer);

  const double tMin = 0;//ns
  const double tStep = (tMax-tMin)/1000.;
  const int nTimeBins = int((tMax-tMin)/tStep);
  sensor->SetTimeWindow(0.,tStep,nTimeBins);

//  std::vector<Avalanche*> Aval_Vect;
/*
  AvalancheMC* aval = new AvalancheMC();
  aval->SetSensor(sensor);
  aval->EnableSignalCalculation();
//  aval->SetTimeSteps(0.05);//us
  aval->SetDistanceSteps(0.1);//cm
//  aval->SetCollisionSteps(1000);
*/

  // This canvas will be used to display the drift lines and the field
  TCanvas * c = new TCanvas("c", "c", 10, 10, 1000, 700);
  c->Divide(2,1);
  // Construct object to visualise drift lines
  ViewDrift* viewdrift = new ViewDrift();
  viewdrift->SetArea(-Electrode_Len/2.,0,-Electrode_Len/2.,Electrode_Len/2.,Frame_Space,Electrode_Len/2.);
  viewdrift->SetClusterMarkerSize(0.1);
  viewdrift->SetCollisionMarkerSize(0.5);
  viewdrift->SetCanvas((TCanvas*)c->cd(1));

// construct avalanche class
  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  aval->EnableSignalCalculation(); 
  aval->SetTimeWindow(tMin,tMax ); 
  aval->SetCollisionSteps(2);
  //  aval->EnableAvalancheSizeLimit(aval_size);
  aval->EnablePlotting(viewdrift);
  aval->EnableDriftLines();
  aval->EnableMagneticField();

  //Electron Cluster info
  double xcls, ycls, zcls, tcls, e, extra;
  xcls = ycls = zcls = tcls = e = extra = -999.;
  // Electron info
  double xele, yele, zele, tele, eele, dxele, dyele, dzele;
  // Electron start and endpoints, momentum direction and status
  int n = 0; // number of electrons in cluster
  bool cls_ret;// return OK if cluster is OK

  // The initial impact position of the incoming ionising track
  double track_x = 0.0;// [cm]
  double track_y = 0.0;
  double track_z = 0.0;
  // Momentum direction of incoming track
  double track_dx = 0.0;
  double track_dy = 0.0; // Track coming downstream
  double track_dz = 1.0;

  double time;//ns
  double location;

  int ne, ni;
  int epoint;

  double dE = Frame_Space*10.*dedx*1e6;//ion energy loss in gas eV
  n = int(dE/W);//number of e- hole pairs formed

  double E = 5.5;//ion kinetic energy MeV
  double E_J = 1.602e-19*1e6*E;//energy in joules
  double m = 4.*1.66e-27;//ion mass kg

  double v = pow(2.*E_J/m,0.5);//m/s
  double t_tot = Frame_Space/(v*10.)*1e9;//time taken for ion to traverse detector gas (ns)

  std::cout <<"# of track n = "<< n << " energy transferred = "<< e <<std::endl;

//  srand (time(NULL));

  for (int j=1; j<=n; j++) {//loops over all evectrons in cascasde

    std::cout <<"# th of track n = "<< j << " is started"<<std::endl;
  // Find radial location of cluster

    xele=0;
    yele=0;
    zele=(rand()%10000)/10000.*Frame_Space;
    dxele=0;
    dyele=0;
    dzele=0;

    tele=zele/Frame_Space*t_tot;

  // Simulate an avalanche for the current electron
    aval->AvalancheElectron(xele, yele, zele, tele, eele, dxele, dyele, dzele);
    std::cout<<"Aval start point x = "<< xele<<" start y = "<<yele<<" start z = "<<zele<<" start time = "<< tele << std::endl;
    aval->GetAvalancheSize(ne, ni);
    epoint = aval->GetNumberOfElectronEndpoints();

    std::cout<<"GetAvalancheSize(int& ne, int& ni); e = "<< ne <<" i = "<< ni <<std::endl;
    std::cout<<"GetNumberOfElectronEndpoints() =  "<< epoint << std::endl; 

  }

  sensor->ConvoluteSignal();
  viewdrift->Plot();
  // Plot isopotential contours
/*
  ViewField* fView = new ViewField;
  fView->SetSensor(sensor);
  fView->SetComponent(cmp);
  fView->SetVoltageRange(0.,Anode_Voltage);
  fView->SetCanvas((TCanvas*)c->cd(2));
  fView->SetWeightingFieldRange(0.0, 10000.0);
  c->cd(2);
  fView->PlotContour();
*/
  ViewSignal* signalView = new ViewSignal();
  signalView->SetSensor(sensor);
  signalView->PlotSignal("c0");
  TH1D* hSignal = signalView->GetSignal();

  std::cout << hSignal->GetBinContent(hSignal->GetMinimumBin()) << std::endl;
  std::cout << hSignal->Integral() << std::endl;

  app.Run(kTRUE);

  return 0;

}
