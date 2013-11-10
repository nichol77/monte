#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include <math.h>
#include <string>
#include "TROOT.h"
#include <vector>
#include <iomanip>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "distances.h"
//#include "/Users/amyc/askSim_pared/NeutrinoEvent.h"
#include "TVector3.h"
//#include "signal.hh"

using namespace std;

#include "shared.hh"
#include "salt.hh"

// inputs 
int VERBOSE=1;
double EXPONENT; // energy exponent
int SEED; // seed
int HIST;  // whether or not to fill histograms and trees.
int SHAPE=0;    // Rectangular (0) or NFOLD-ogonal (1)
int DETECTOR=0; // 0=SalSA, 1=Salsita, 2=one hole
int TYPEOFANTENNAS; // antenna type  
//const double MAXXPOS=2000.;  // +- limits of possible interactions
//const double MAXYPOS=2000.;
//const double MAXZPOS=2000.;
//const double MAXPOS[3]={4000.,4000.,4000.};
double MAXPOS[3]={2000.,2000.,5000.}; // dimensions of the salt
//const double MAXPOS[3]={2000.,2000.,1500.};
int NRXX=10; // number of strings in the x dimension
int NRXY=10; // number of strings in the y dimension
int NHEX=3; // number of hexagons

int NNODESEACHSTRING=12; // number of nodes in the z direction per string
int NNODESTOTAL=NRXX*NRXY*NNODESEACHSTRING; // total number of nodes in the dettector
int NEACHNODE=3; // number of antennas in each module
int NSTRING=0; // how many "strings" on in the detector.  This is really how many modules in a x-y cross section
int NANTENNA=0; // total number of antennas
// nfold and nrad or for a bicycle wheel shaped detector, which is not working yet.
const int NFOLD=6;  // detector is NFOLD-fold symmetric.
// if NFOLD=6, detector is hexagonal.
const int NRAD=10;  // number of strings radially from center.

const double SPACINGWITHINMODULE=0.75; // separation of antennas between modules
int FOLLOWONEEVENT=1; // print out interesting quantities for one event, for debugging
int DEPTH_DEPENDENT_N=1; // whether or not to include a depth dependent index of refraction
int CONSTANTY=0; // whether or not to use a constant inelasticity y
int FORSECKEL=1; // Make array of strength of signal across frequencies for different viewing angles.
int SHOWERTYPE=0; // Type of shower for previous option
int SECONDARIES=0;  // whether or not to include secondaries
int ONLYFINAL=1; // only write to final histogram
int HIST_MAX_ENTRIES=10000; //maximum number of events to put in histograms
int TAUDECAY=1; // is tau decay counted as a secondary interaction
int SIDEWAYS=1; // whether or not to orient the dipoles/slots sideways
double MINE_DEPTH = 500.;    // depth of detector
//const double INSTRUMENTEDREGION_DEPTH=2500.; // depth of instrumented region.
double INSTRUMENTEDREGION_DEPTH=750.; // depth of instrumented region.
//  const double MAXARRAY=500.;   // +- limits of array
// 9 spaces*250m =2250m-> 2250/2=1125m (horizontally)
// 11 spaces*182m=2000m-> 2000/2=1000 (vertically)
double MAXARRAY[3]={1125.,1125.,1000.};// dimensions of the detector volume measured from the center of the array.   So if MAXARRAY[2]=1000., for example, then the detection volume is 2000 meters vertically
double LIVETIME=0.5*365.*24.*3600.; // livetime of the experiment for calculating the number of events expected.  Here it's half a year.

const bool USEBELOWHORIZ = true;      // include events below horizon?

// PG had 3 nodes and 4 ant
int NNODEHITSREQUIRED=5; // # of modules required to be hit
double DISTCUT=1.E3; // max distance between a cluster of hit nodes
int irxhitsrequired=5; // # of antennas required to be hit within a module
const double isigmarequired=2.8;
//const double isigmarequired=3.0;
//const double isigmarequired=2.5;

// for counting
int count_neutrinos=0; // count neutrinos that we loop over.
int count_getmine=0; // passes Getmine (neutrino's path through salt is more than 1 meter)
int count_getchord=0; // passes Getchord (they should all pass.)
int count_getchordweight=0; // passes Getchordweight (neutrino path length is less than 50 interaction lengths)
int count_viewangle_chanceinhell=0; // count antennas where the viewing angle is less than 20 cone widths from the cherenkov angle
int count_atten_factor_chanceinhell=0;  // count antennas that are less than 10 attenuation lengths away 
int count_passing_events=0;
int increment_viewangle_chanceinhell=0;
int increment_atten_factor_chanceinhell=0;
int count_anthit=0;
int count_modulehit[5]={0,0,0,0,0};
int count_passestrigger[5]={0,0,0,0,0};
int maxnnodehitsincluster;

// secondaries
int NPROB=100;
int secondary_e_noncons=0;
const int NPROB_MAX=200;
double pnu;
double plepton=pnu;  // energy of charged lepton after interaction(s)
double em_secondaries_max=0;  // em energy due to secondary interactions
double had_secondaries_max=0; // had energy due to secondary interactions
// first index=energy from 10^18 to 10^21 in half-decades
// second index=y (100 bins)
double dsdy_muon_brems[7][NPROB_MAX]; // probability distribution vs. y for brem
double dsdy_muon_epair[7][NPROB_MAX]; // prob. dist. vs. y for pair production
double dsdy_muon_pn[7][NPROB_MAX]; // photonuclear interactions

double y_muon_brems[7][NPROB_MAX]; // inelasticity corresponding to each bin,brem
double y_muon_epair[7][NPROB_MAX]; // same for pair production
double y_muon_pn[7][NPROB_MAX]; // same for photonuclear

vector<double> vy_muon_brems[7]; // vectors containing inelasticities distributed such
                                //that they follow the dsdy curves above.
vector<double> vy_muon_epair[7]; 
vector<double> vy_muon_pn[7]; 

double y_cumulative_muon_brems[7][NPROB_MAX];
double y_cumulative_muon_epair[7][NPROB_MAX];
double y_cumulative_muon_pn[7][NPROB_MAX];

double int_muon_brems[7];  // integral of prob. dist.=# of expected interactions.
                           // for each energy
double int_muon_epair[7]; // same for pair prod.
double int_muon_pn[7]; // same for photnuclear

double max_muon_brems=1000.; // maximum value of prob. dist., brems
double max_muon_epair=1000.;  // max value of prob. dist., pair prod.
double max_muon_pn=1000.; // same for photonucl.

double min_muon_brems=1000.; // minimum value of prob. dist., brems
double min_muon_epair=1000.;  // min value of prob. dist., pair prod.
double min_muon_pn=1000.; // same for photonucl.

double dsdy_tauon_brems[7][NPROB_MAX];  // same as above, but for taus
double dsdy_tauon_epair[7][NPROB_MAX];
double dsdy_tauon_pn[7][NPROB_MAX];
double dsdy_tauon_hadrdecay[7][NPROB_MAX]; // hadronic decay of taus
double dsdy_tauon_edecay[7][NPROB_MAX]; // tau decay to electrons
double dsdy_tauon_mudecay[7][NPROB_MAX]; // tau decay to muons.

double y_tauon_brems[7][NPROB_MAX];
double y_tauon_epair[7][NPROB_MAX];
double y_tauon_pn[7][NPROB_MAX];
double y_tauon_hadrdecay[7][NPROB_MAX];
double y_tauon_edecay[7][NPROB_MAX];
double y_tauon_mudecay[7][NPROB_MAX];

double y_cumulative_tauon_brems[7][NPROB_MAX];
double y_cumulative_tauon_epair[7][NPROB_MAX];
double y_cumulative_tauon_pn[7][NPROB_MAX];
double y_cumulative_tauon_hadrdecay[7][NPROB_MAX];
double y_cumulative_tauon_edecay[7][NPROB_MAX];
double y_cumulative_tauon_mudecay[7][NPROB_MAX];

vector<double> vy_tauon_brems[7]; // vectors containing inelasticities distributed such
                                //that they follow the dsdy curves above.
vector<double> vy_tauon_epair[7]; 
vector<double> vy_tauon_pn[7]; 
vector<double> vy_tauon_hadrdecay[7];
vector<double> vy_tauon_edecay[7];
vector<double> vy_tauon_mudecay[7];

double int_tauon_brems[7];
double int_tauon_epair[7];
double int_tauon_pn[7];
double int_tauon_hadrdecay[7];
double int_tauon_edecay[7];
double int_tauon_mudecay[7];

double max_tauon_brems=1000.;
double max_tauon_epair=1000.;
double max_tauon_pn=1000.;
double max_tauon_hadrdecay=1000.;
double max_tauon_edecay=1000.;
double max_tauon_mudecay=1000.;

double min_tauon_brems=1000.; // minimum value of prob. dist., brems
double min_tauon_epair=1000.;  // min value of prob. dist., pair prod.
double min_tauon_pn=1000.; // same for photonucl.
double min_tauon_hadrdecay=1000.; // minimum value of prob. dist., brems
double min_tauon_edecay=1000.;  // min value of prob. dist., pair prod.
double min_tauon_mudecay=1000.; // same for photonucl.





double vmmhz1m_perfreq[NFREQ_MAX];  // V/m/MHz per frequency
double freqerr[NFREQ_MAX]; // error for tgrapherrors
double yerr[NFREQ_MAX]; // error for tgrapherrors
double poyntingflux_1m; // power per unit area at 1 m
double forseckel[NVIEWANGLE][NFREQ_MAX];// Per Seckel's request, get strength of signal across frequencies for different viewing angles.
//  double viewangles[NVIEWANGLE]={acos(1/NSALT),
//  			       acos(1/NSALT)-2.*RADDEG,
//  			       acos(1/NSALT)-5.*RADDEG,
//  			       acos(1/NSALT)-10.*RADDEG,
//  			       acos(1/NSALT)-20.*RADDEG,
//  			       acos(1/NSALT)-30.*RADDEG,
//  			       90.*RADDEG};
double viewangles[NVIEWANGLE];


double vmmhz_temp;  // for making the forseckel array



fstream ftrigeff("trig.txt",ios::out);

double Square(double*,int);
int Min(int x,int y);
double Max(double x,double y);


string input="inputs.txt";
ifstream inputsfile(input.c_str());
ofstream fout("output.txt", ios::app);
ofstream fvmmhz1m("vmmhz1m.txt");
void GetNextNumberAsString(string& number);
void ReadInputs();
void ReadSecondaries(); // read in prob. dist. for secondary interactions

double GetMax(double *x,int n);
double GetMin(double *x,int n);

int TIR(double *n_surf,double *nrf2_iceside, double N_IN,double N_OUT);
int GetRayIceSide(double *n_exit2rx,
		  double *nsurf_rfexit,
		  double nexit,
		  double nenter,

		  double *nrf2_iceside);
double GetNSouthPole(double depth);
double GetNRonneIceShelf(double altitude);// Get index of refraction for a specific depth for the Ronne Ice Shelf
void EarthCurvature(double* array,double depth_temp); // adjusts coordinates within the mine to account for the curvature of the earth.
void GetDepth(double* posnu,double& depth); // input posnu array in coordinates where the origin is the center of the mine and it returns the depth relative to the top of the rock.
double depth_temp=0; // temporary variable used to find a depth in coordinates where the origin is the center of the mine and it returns the depth relative to the top of the rock.
double altitude_interaction=0; // altitude of the interaction relative to the surface (a negative number)
int istep_depth;
double costhray_depth;
double costhray_depth_beforerayiceside;
double depth_thisantenna=0.;

double maxflux;
double templength;
double temp[3];  // just used for taking cross products
// direction cross nnu cross direction

double signal_pol[3]; // polarization of the signal

//const int NSTRINGS_MAX=10001;           // Number of strings strings
//const int NMODULES_MAX=NSTRINGS_MAX*12; // 12 nodes per string
//const int NANTENNA_MAX=NMODULES_MAX*12; // 12 antenna per node


int nrays_detected=0; // number of rays (direct+reflected) that are detected.

double FREQ_LOW; // this should be lower than the lowest of freq_low_dipolesslots and freq_low_seaveys.  it is entered from the input file
double BW; // this should span the widest bandwidth 

double FREQ_LOW_DIPOLESSLOTS=200.E6;// low edge of bandwidth of dipoles and slots in air
double BW_DIPOLESSLOTS=400.E6; // antenna bandwidth of dipoles and slots in air


double FREQ_LOW_SEAVEYS=200.E6;// low edge of bandwidth of seaveys in air 
double BW_SEAVEYS=1000.E6; // antenna bandwidth of seaveys in air //KAR move this
double VNOISE_DIPOLESSLOTS;
double VNOISE_SEAVEYS;



int WHICHPARAMETERIZATION; // which parameterization, jaime's older papers or jaime's Dec. 2006 paper

int main(int argc, char **argv) {


  //  Signal *sig1=new Signal();

 cout << "I'm here.\n";

  //--------------------------------------------------------------
  //  MC a salt neutrino detector
  //
  // 10/23/03 Start w/ David's fortran code and convert it to c++
  // 
  //--------------------------------------------------------------



  double polarfactor; // effective of polarization on signal



  for (int i=0;i<NFREQ;i++) { // for making tgrapherrors
    freqerr[i]=0;
    yerr[i]=0;
  } 


  Zero(posnu_reltoearthcenter,3);
  Zero(posnu_reltominecenter,3);

  Zero(posnu_down,3);//the position of the mirror point of interaction.wufan

  MAX_LOGWEIGHT=0.;
  MIN_LOGWEIGHT=-1.;

  // local declarations

  const bool WEIGHTCHORDS=true;       // give weights to chords
  const bool ABSORBINEARTH=true;      // absorb neutrinos in earth 

  double rnd,rndlist[10];             // random numbers
  
  double posnux1,posnuy1,posnuz1; // for inserting into tree - position of interaction
  double posnux2,posnuy2,posnuz2; // for inserting into tree - position of mirror interaction


  //int MINMODULEZ,MAXMODULEZ;// for finding range of modules such that we
  // only loop over antennas within 10 radiation lengths.

  double nnu[3];               // direction of neutrino
  int inu;                            // for looping over neutrinos
 
  
     
  // note, lambda is "in salt" so divided by nsalt
  double thisfreq;// frequency

  double vmmhz1m[2];                      // V/m/MHz at 1m for each ray
  double emfrac,hadfrac;               // taking care of LPM
  double vmmhz[2];                        //  V/m/MHz
  double vmmhz0,vmmhz1;  // for inserting into tree
  double vmmhz0_max,vmmhz1_max; // for inserting into tree
  double volts[2];                        // antenna output
  double volts0,volts1; // for inserting into tree
  double volts0_max,volts1_max; // for inserting into tree
  double atten_factor[2];                 // factor to account for attenuation length for direct ray, reflected ray
  double atten_factor0,atten_factor1;
  double heff;                         // antenna effective height


  int ray_thishit;
  double viewangle[2];                    // angle between shower direction
  // and int->ant direction (for direct ray, reflected ray)
  double viewangle0,viewangle1;
  double hitangle_e,hitangle_h,e_component,h_component;


  

  double vnoise;                       // 1 sigma noise level
  int nrxhit_total=0; // total antennas hit in an event
  int nrxhit_dipole=0; // total dipoles hit in an event
  int nrxhit_x=0; // total +x slots hit in an event
  int nrxhit_y=0; // total +y slots hit in an event
  vector<int> nrxhit[3];        // number of antennas hit (module)
  int nmoduleshit[3];           // number of modules hit for one event
  int nmoduleshit_fortree; // put the number of modules hit (reflected + direct) into a tree
  // used to count by string but it never really meant anything
  int nhittotal[3]={0,0,0};     

  double eventsfound[3]={0.,0.,0.};   // how many events found
  double doublesfound=0;  //how many where both rays pass
  double eventsfound_below; // how many below the horizon.
  double eventsfound_onehit=0; // how many satisfy David Besson's one-hit trigger
  

  //double whichpol=0;
  
  double sigma;                       // for cross section
 
  string  nuflavor;                   // neutrino flavor
  string  current;                    //  CC or NC?
  string  taudecay;                   // tau decay type: e,m,h
  int nuflavorint = 0;                // For output purposes

    
  double y; 
  const bool DEBUG =true;                   // debugging option
     

  void Summarize(double,double* eventsfound,double doublesfound,double,double,double,double);
  void CloseTFile(TFile *hfile);
  int IsAbsorbed(double,double);
  void IsAbsorbed(double,double,double&);
  void AddNoise(int cantid,double *volts); // add the noise fluctuation to the signal

  double nabsorbed_earth = 0;
  double nabsorbed_noearth = 0;
  double nweighted;
  double mine_in[3];
  double mine_out[3];
  double mine[3];
     
  double earth_in[3];
  double earth_out[3];

  double firn_in[3];
  double firn_start[3];

  double maxdepth_firn;
  double stepsize_depth;
  double depth_thisstep;
  int NSTEPS_DEPTH=1000;
  double n_depth,n_depth_previous;

  int npass_noweight = 0;

  int Getmine(double*,double*, double*, double*);
  void GetEntryExit(double,double*,double*, double*, double*);
  int Getchord(double*,double*,double&, double&); // get chord lengths

  // weighting chords
  double d1,d2,chord_kgm2,chord;
  double chord_kgm2_noearth;
  double weight1=1;
  double weight2=1.;
  double weight3=1.;
     
  double remaining=0.;
  double len_int_kgm2;
  int GetChordWeight(double,double,double,double,double,double,double&,double&);
     
  int counting=0;

 

  // input parameters
  ReadInputs();
  
  // calculate 1 sigma noise
  VNOISE_DIPOLESSLOTS=sqrt(KBOLTZ*KELVINS*BW/N_RECEIVER*50.);
  VNOISE_SEAVEYS=sqrt(KBOLTZ*KELVINS*BW_SEAVEYS/N_RECEIVER*50.);




  if (FORSECKEL==1) {
    double viewangle_max=90.*RADDEG;
    double viewangle_min=30.*RADDEG;
    for (int i=0;i<NVIEWANGLE-2;i++) {
      viewangles[i]=viewangle_max-(viewangle_max-viewangle_min)/(double)(NVIEWANGLE-2)*(double)i;
    }
    viewangles[NVIEWANGLE-2]=acos(1/NMEDIUM);
    viewangles[NVIEWANGLE-1]=90.*RADDEG;
  }


  for (int i=0;i<NNODESTOTAL;i++) {
    for (int p=0;p<3;p++) {
      nrxhit[p].push_back(0); 
    }  
  }
  
  double freqlist[NFREQ_MAX];
  
  for (int ifreq=0;ifreq<NFREQ;ifreq++) {
    if (NFREQ==1) {
      freqlist[0]=200.E6;
      break;
    }
    if (NFREQ==2) {
      freqlist[0]=150.E6;
      freqlist[1]=250.E6;
      break;
    }
    if (NFREQ>2)
      freqlist[ifreq]=(FREQ_LOW+BW*ifreq/(NFREQ-1))/N_RECEIVER;
  }     

  GetBeamWidths(flare,gain,freqlist);
  
  // Antenna measured gain vs. frequency
  ReadGains();
  

   


  if (WHICHRAYS==1) { // just count direct rays
    MINRAY=0;
    MAXRAY=0;
  } //if (WHICHRAYS==1)
  if (WHICHRAYS==2) {// direct and reflected
    MINRAY=0;
    MAXRAY=1;
  } //if (WHICHRAYS==2) 
  if (WHICHRAYS==3) { // just reflected
    MINRAY=1;
    MAXRAY=1;
  } //if (WHICHRAYS==3
 

  cout << "Freq list ";
  for (int ifreq=0;ifreq<NFREQ;ifreq++)
    cout << freqlist[ifreq] << " ";
  cout << endl;


  void GetDetector(vector<double> rxpos[3], vector<int> rxtype[4]);
  //int WhichModule(double* posnu,int& MINMODULEZ, int& MAXMODULEZ);

  double costhetanu,sinthetanu,thetanu,phinu;      // neutrino direction

  double deltheta_em_max,deltheta_had_max; // for putting width of cone into tree
  double deltheta_em[NFREQ_MAX];
  double deltheta_had[NFREQ_MAX];

  vector<double> rxpos[3];                 // locations of antennas
  // 1,2,3=xpos,ypos,zpos
  vector<int> rxtype[4]; // info on each antenna
  // 1=orientation (0,1,2) 
  // 2=module#
  // 3=string# 
  // 4=polarization type (dipole slot etc)

  int irx;                            // for looping over antennas
  //vector<double> rxpos[NSTRINGS_MAX][NEACHNODE_MAX][3];
  //double rxpos[NSTRINGS_MAX][NEACHNODE][3][NNODES_MAX];                 // locations of antennas
  // 1,2,3=xpos,ypos,zpos
  double distance[2];                  // dist. from interaction to ant. along ray's path for direct, reflected rays
  double distance0,distance1;


  double rx_reltominecenter[3];                  // antenna position
  double rx_reltoearthcenter[3]; // antenna position rel to earth center
  double xhit,yhit,zhit;         //for filling tree
  double direction[2][3];        // direction from nu int. to antenna for direct, reflected rays. For the latter it is the upcoming ray after bouncing off the bottom.
  double backwards[2][3];        // minus the direction from nu int. to antenna for direct, reflected rays. For the latter it is the upcoming ray after bouncing off the bottom.
  double backwards_next[2][3];        // minus the direction of the ray at the top of the firn
  int tir[2]={0,0};
  int tir_thisantenna;
  double position_bent[2][3]; // keeps track of the position of the ray as it gets bent
  TVector3 vposition_bent[2];
  TVector3 vnsurf;
  TVector3 vbackwards[2];
  TVector3 vposnu_reltoearthcenter;
  TVector3 vfirn_start;
  TVector3 vdiff;
  TVector3 vtemp;
  TVector3 vperp;


  double xpos,zpos;

  double direction00,direction01,direction02; // direction vectors made into doubles for filling into trees
  double direction10,direction11,direction12;
 
  double nsurf[3];  // "surface" normal at interface between firn and ice

  double downwardray[3]; // direction of downward going ray, for events where the ray gets reflected
  double phiray[2],costhray[2]; // direction of direct ray as it hits the antenna for direct rays, reflected rays
  double costhray0,costhray1;
  double costhray0_atinteraction,costhray1_atinteraction;


  for (int i=0;i<3;i++) {
    rx_reltominecenter[i]=0.;
    rx_reltoearthcenter[i]=0.;
  }
  int nrx_fired[3]; // number of antennas fired for each ray for one event
  int ndipoles_fired[3]; // number of dipoles fired for each ray for one event
  int nslotsx_fired[3]; // number of slots (x) fired for each ray for one event
  int nslotsy_fired[3]; // number of slots (y) fired for each ray for one event
  int nrx_fired0,nrx_fired1; // for inserting into tree

  int nrx_fired_total[3]; // number of antennas fired for each ray for one event
  int ndipoles_fired_total[3]; // number of dipoles fired for each ray for one event
  int nslotsx_fired_total[3]; // number of slots (x) fired for each ray for one event
  int nslotsy_fired_total[3]; // number of slots (y) fired for each ray for one event
  int nhornsx_fired_total[3]; // number of horns (x polarization) fired for each ray for one event
  int nhornsy_fired_total[3]; // number of horns (y polarization) fired for each ray for one event

  int nmoduleshit_total[3]; // number of modules hit total


  double costh_pol[2]={0.,0.};  // cos(theta) where theta is the zenith angle of the polarization vector that triggers the antenna for direct, reflected rays.

  double costh_dir[2]={0.,0.};  // cos(theta) where theta is the zenith angle of the direction from the interaction to the antenna for direct, reflected rays.

  for (int p=0;p<3;p++) { // this counts the number of hits that fire for each ray in each event.
    nmoduleshit_total[p]=0;
    nrx_fired_total[p]=0;
    ndipoles_fired_total[p]=0;
    nslotsx_fired_total[p]=0;
    nslotsy_fired_total[p]=0;
    nhornsx_fired_total[p]=0;
    nhornsy_fired_total[p]=0;


  }

     
     
  if (EXPONENT>0)
    pnu=TMath::Power(10,EXPONENT);                          // neutrino momentum
  const double checkvalue=(int)((double)NNU*0.10);
  //-------------------
  // stuff for root
  //---------------
  double weight=0;
  TFile *hfile = new TFile("salt.root","RECREATE","salt");
  TGraphErrors *g1; // graph for looking at vmmhz1m at 1m


  TTree *tree_passing_events = new TTree("passing_events","passing_events");
  TTree *tree_beforetrigger = new TTree("beforetrigger","beforetrigger");
  TTree *tree_hitantennas = new TTree("hitantennas","hitantennas");
  TTree *tree_allantennas = new TTree("allantennas","allantennas");
  TTree *tree_steps=new TTree("steps","steps");
 
  tree_steps->Branch("inu",&inu,"inu/I",8*1024*1024);
  tree_steps->Branch("irx",&irx,"irx/I",8*1024*1024);
  tree_steps->Branch("xpos",&xpos,"xpos/D",8*1024*1024);
  tree_steps->Branch("zpos",&zpos,"zpos/D",8*1024*1024);
  tree_steps->Branch("istep_depth",&istep_depth,"istep_depth/I",8*1024*1024);
  tree_steps->Branch("costhray_depth",&costhray_depth,"costhray_depth/D",8*1024*1024);
  tree_steps->Branch("costhray_depth_beforerayiceside",&costhray_depth_beforerayiceside,"costhray_depth_beforerayiceside/D",8*1024*1024);

  tree_steps->Branch("costhray0_atinteraction",&costhray0_atinteraction,"costhray0_atinteraction/D",8*1024*1024);
  tree_steps->Branch("depth_thisstep",&depth_thisstep,"depth_thisstep/D",8*1024*1024);


  tree_passing_events->Branch("nmoduleshit",&nmoduleshit_fortree,"nmoduleshit_fortree/I",8*1024*1024);
  tree_passing_events->Branch("inu",&inu,"inu/I",8*1024*1024);
  tree_passing_events->Branch("whichray",&whichray,"whichray/I",8*1024*1024);
  tree_passing_events->Branch("phinu",&phinu,"phinu/D",8*1024*1024);
  tree_passing_events->Branch("phiray0",&phiray[0],"phiray0/D",8*1024*1024);
  tree_passing_events->Branch("costhray0",&costhray0,"costhray0/D",8*1024*1024);
  tree_passing_events->Branch("costhray1",&costhray1,"costhray1/D",8*1024*1024);
  tree_passing_events->Branch("thetanu",&thetanu,"thetanu/D",8*1024*1024);

  tree_passing_events->Branch("weight1",&weight1,"weight1/D",8*1024*1024);
  tree_passing_events->Branch("weight2",&weight2,"weight2/D",8*1024*1024);
  tree_passing_events->Branch("weight3",&weight3,"weight3/D",8*1024*1024);
  tree_passing_events->Branch("weight",&weight,"weight/D",8*1024*1024);
  tree_passing_events->Branch("posnux1",&posnux1,"posnux1/D",8*1024*1024);  // position of interaction relative to center of mine
  tree_passing_events->Branch("posnuy1",&posnuy1,"posnuy1/D",8*1024*1024);
  tree_passing_events->Branch("posnuz1",&posnuz1,"posnuz1/D",8*1024*1024);

  tree_passing_events->Branch("posnux2",&posnux2,"posnux2/D",8*1024*1024);  // position of interaction relative to center of mine
  tree_passing_events->Branch("posnuy2",&posnuy2,"posnuy2/D",8*1024*1024);
  tree_passing_events->Branch("posnuz2",&posnuz2,"posnuz2/D",8*1024*1024);
  tree_passing_events->Branch("nuflavor",&nuflavorint,"nuflavorint/I");//1=electron, 2=muon, 3=tau
  tree_passing_events->Branch("chord",&chord,"chord/D",8*1024*1024);
  tree_passing_events->Branch("d1",&d1,"d1/D",8*1024*1024);
  tree_passing_events->Branch("d2",&d2,"d2/D",8*1024*1024);
  tree_passing_events->Branch("emfrac",&emfrac,"emfrac/D",8*1024*1024);
  tree_passing_events->Branch("hadfrac",&hadfrac,"hadfrac/D",8*1024*1024);
  tree_passing_events->Branch("pnu",&pnu,"pnu/D",8*1024*1024);
  tree_passing_events->Branch("doublesfound",&doublesfound,"doublesfound/D",8*1024*1024);
  tree_passing_events->Branch("altitude_interaction",&altitude_interaction,"altitude_interaction/D",8*1024*1024);
  

  tree_beforetrigger->Branch("inu",&inu,"inu/I",8*1024*1024);
  tree_beforetrigger->Branch("phinu",&phinu,"phinu/D",8*1024*1024);
  tree_beforetrigger->Branch("vmmhz0_max",&vmmhz0_max,"vmmhz0_max/D",8*1024*1024);
  tree_beforetrigger->Branch("vmmhz1_max",&vmmhz1_max,"vmmhz1_max/D",8*1024*1024);
  tree_beforetrigger->Branch("volts0_max",&volts0_max,"volts0_max/D",8*1024*1024);
  tree_beforetrigger->Branch("volts1_max",&volts1_max,"volts1_max/D",8*1024*1024);
  tree_beforetrigger->Branch("thetanu",&thetanu,"thetanu/D",8*1024*1024);
  tree_beforetrigger->Branch("nrx_fired0",&nrx_fired0,"nrx_fired0/I",8*1024*1024);
   tree_beforetrigger->Branch("nrx_fired1",&nrx_fired1,"nrx_fired1/I",8*1024*1024);
   tree_beforetrigger->Branch("changle",&changle,"changle/D",8*1024*1024);
   tree_beforetrigger->Branch("deltheta_em_max",&deltheta_em_max,"deltheta_em_max/D",8*1024*1024);
   tree_beforetrigger->Branch("deltheta_had_max",&deltheta_had_max,"deltheta_had_max/D",8*1024*1024);
   tree_beforetrigger->Branch("weight",&weight,"weight/D",8*1024*1024);
   tree_beforetrigger->Branch("maxnnodehitsincluster",&maxnnodehitsincluster,"maxnnodehitsincluster/I",8*1024*1024);
  tree_beforetrigger->Branch("posnux1",&posnux1,"posnux1/D",8*1024*1024);  // position of interaction relative to center of mine
  tree_beforetrigger->Branch("posnuy1",&posnuy1,"posnuy1/D",8*1024*1024);
  tree_beforetrigger->Branch("posnuz1",&posnuz1,"posnuz1/D",8*1024*1024);   

  tree_allantennas->Branch("NANTENNA",&NANTENNA,"NANTENNA/I",8*1024*1024);
  tree_allantennas->Branch("xhit",&xhit,"xhit/D",8*1024*1024);
  tree_allantennas->Branch("yhit",&yhit,"yhit/D",8*1024*1024);
  tree_allantennas->Branch("zhit",&zhit,"zhit/D",8*1024*1024);

  tree_allantennas->Branch("irx",&irx,"irx/I",8*1024*1024);
  tree_allantennas->Branch("inu",&inu,"inu/I",8*1024*1024);
  tree_allantennas->Branch("phiray0",&phiray[0],"phiray0/D",8*1024*1024);
  tree_allantennas->Branch("vmmhz0",&vmmhz0,"vmmhz0/D",8*1024*1024);
  tree_allantennas->Branch("vmmhz1",&vmmhz1,"vmmhz1/D",8*1024*1024);
  tree_allantennas->Branch("volts0",&volts0,"volts0/D",8*1024*1024);
  tree_allantennas->Branch("volts1",&volts1,"volts1/D",8*1024*1024);


  tree_allantennas->Branch("viewangle0",&viewangle0,"viewangle0/D",8*1024*1024);
  tree_allantennas->Branch("viewangle1",&viewangle1,"viewangle1/D",8*1024*1024);
  tree_allantennas->Branch("atten_factor0",&atten_factor0,"atten_factor0/D",8*1024*1024);
  tree_allantennas->Branch("atten_factor1",&atten_factor1,"atten_factor1/D",8*1024*1024);
  tree_allantennas->Branch("tir_thisantenna",&tir_thisantenna,"tir_thisantenna/I",8*1024*1024);
  tree_allantennas->Branch("costhray0",&costhray0,"costhray0/D",8*1024*1024);
  tree_allantennas->Branch("costhray1",&costhray1,"costhray1/D",8*1024*1024);
  tree_allantennas->Branch("costhray0_atinteraction",&costhray0_atinteraction,"costhray0_atinteraction/D",8*1024*1024);
  tree_allantennas->Branch("costhray1_atinteraction",&costhray1_atinteraction,"costhray1_atinteraction/D",8*1024*1024);
  

  Int_t cantid,anthit,cmodule,cstring,cantpolar;
  TTree *reil1 = new TTree("mcinfo","mcinfo");
  TTree *reil2 = new TTree("rawdata","rawdata");
  Int_t reilinu, reilray, reilantnum,reilhit;
  Double_t reilpnu,reilnux,reilnuy,reilnuz,reilem,reilhad;
  Double_t reilnnu0,reilnnu1,reilnnu2;
  Double_t reilvolts,reilvnoise,reiltime,reiltnoise;
  Double_t fluctime;

  int anyanthit=0; // for each event, keeps track of whether any antenna was hit
  int anymodulehit[5]={0,0,0,0,0}; // for each event, keeps track of whether any module was hit
  Double_t dtime;
  if (HIST==1) {
    tree_hitantennas->Branch("inu",&inu,"inu/I",8*1024*1024);
    tree_hitantennas->Branch("direction00",&direction00,"direction00/D",8*1024*1024);
    tree_hitantennas->Branch("direction01",&direction01,"direction01/D",8*1024*1024);
    tree_hitantennas->Branch("direction02",&direction02,"direction02/D",8*1024*1024);
    tree_hitantennas->Branch("direction10",&direction10,"direction10/D",8*1024*1024);
    tree_hitantennas->Branch("direction11",&direction11,"direction11/D",8*1024*1024);
    tree_hitantennas->Branch("direction12",&direction12,"direction12/D",8*1024*1024);
    tree_hitantennas->Branch("cantpolar",&cantpolar,"cantpolar/I",8*1024*1024);
    
    
    tree_hitantennas->Branch("xhit",&xhit,"xhit/D",8*1024*1024);
    tree_hitantennas->Branch("yhit",&yhit,"yhit/D",8*1024*1024);
    tree_hitantennas->Branch("zhit",&zhit,"zhit/D",8*1024*1024);
    tree_hitantennas->Branch("volts0",&volts0,"volts0/D",8*1024*1024);
    tree_hitantennas->Branch("vnoise",&VNOISE_DIPOLESSLOTS,"vnoise/D",8*1024*1024);
    tree_hitantennas->Branch("dtime",&dtime,"dtime/D",8*1024*1024);
    tree_hitantennas->Branch("antid",&cantid,"antid/I",8*1024*1024);
    tree_hitantennas->Branch("antpolar",&cantpolar,"antpolar/I",8*1024*1024);
    tree_hitantennas->Branch("anthit",&anthit,"anthit/I",8*1024*1024);
    tree_hitantennas->Branch("module",&cmodule,"module/I",8*1024*1024);
    tree_hitantennas->Branch("string",&cstring,"string/I",8*1024*1024);
    tree_hitantennas->Branch("anthit",&anthit,"anthit/I",8*1024*1024);	
    tree_hitantennas->Branch("costh_pol0",&costh_pol[0],"costh_pol0/D",8*1024*1024);
    tree_hitantennas->Branch("costh_pol1",&costh_pol[1],"costh_pol1/D",8*1024*1024);
    
    tree_hitantennas->Branch("costh_dir0",&costh_dir[0],"costh_dir0/D",8*1024*1024);
    tree_hitantennas->Branch("costh_dir1",&costh_dir[1],"costh_dir1/D",8*1024*1024);
    tree_hitantennas->Branch("weight",&weight,"weight/D",8*1024*1024);
    tree_hitantennas->Branch("phiray0",&phiray[0],"phiray0/D",8*1024*1024);
    tree_hitantennas->Branch("costhray0",&costhray0,"costhray0/D",8*1024*1024);
    tree_hitantennas->Branch("costhray1",&costhray1,"costhray1/D",8*1024*1024);
    tree_hitantennas->Branch("costhray0_atinteraction",&costhray0_atinteraction,"costhray0_atinteraction/D",8*1024*1024);
    tree_hitantennas->Branch("costhray1_atinteraction",&costhray1_atinteraction,"costhray1_atinteraction/D",8*1024*1024);
    tree_hitantennas->Branch("altitude_interaction",&altitude_interaction,"altitude_interaction/D",8*1024*1024);

    


    tree_hitantennas->Branch("atten_factor0",&atten_factor0,"atten_factor0/D",8*1024*1024);
    tree_hitantennas->Branch("atten_factor1",&atten_factor1,"atten_factor1/D",8*1024*1024);
    tree_hitantennas->Branch("viewangle0",&viewangle0,"viewangle0/D",8*1024*1024);
    tree_hitantennas->Branch("viewangle1",&viewangle1,"viewangle1/D",8*1024*1024);
    tree_hitantennas->Branch("distance0",&distance0,"distance0/D",8*1024*1024);
    tree_hitantennas->Branch("distance1",&distance1,"distance1/D",8*1024*1024);
    
    
    //    tree_hitantennas->Branch("whichpol",&whichpol,"whichpol/D",8*1024*1024);
    tree_hitantennas->Branch("ray_thishit",&ray_thishit,"ray_thishit/I",8*1024*1024);
    
    tree_hitantennas->Branch("posnux1",&posnux1,"posnux1/D",8*1024*1024);  // position of interaction relative to center of mine
    tree_hitantennas->Branch("posnuy1",&posnuy1,"posnuy1/D",8*1024*1024);
    tree_hitantennas->Branch("posnuz1",&posnuz1,"posnuz1/D",8*1024*1024);
    
    tree_hitantennas->Branch("posnux2",&posnux2,"posnux2/D",8*1024*1024);  // position of interaction relative to center of mine
    tree_hitantennas->Branch("posnuy2",&posnuy2,"posnuy2/D",8*1024*1024);
    tree_hitantennas->Branch("posnuz2",&posnuz2,"posnuz2/D",8*1024*1024);
    
       
    // neutrino number
    reil1->Branch("inu",&reilinu,"inu/I",8*1024*1024);
    // position of interaction relative to center of mine
    reil1->Branch("posnux",&reilnux,"posnux/D",8*1024*1024);  
    reil1->Branch("posnuy",&reilnuy,"posnuy/D",8*1024*1024);
    reil1->Branch("posnuz",&reilnuz,"posnuz/D",8*1024*1024);
    // neutrino energy and branching
    reil1->Branch("pnu",&reilpnu,"pnu/D");
    reil1->Branch("emfrac",&reilem,"emfrac/D");
    reil1->Branch("hadfrac",&reilhad,"hadfrac/D");
    // neutrino direction
    reil1->Branch("nnu0",&reilnnu0,"nnu0/D");
    reil1->Branch("nnu1",&reilnnu1,"nnu1/D");
    reil1->Branch("nnu2",&reilnnu2,"nnu2/D");
    // the "ray number" of the hit
    reil1->Branch("pray",&reilray,"pray/I",8*1024*1024);
    
    // neutrino number
    reil2->Branch("inu",&reilinu,"inu/I",8*1024*1024);
    // volts
    reil2->Branch("volts",&reilvolts,"volts0/D",8*1024*1024);
    // flucuation used
    reil2->Branch("vnoise",&reilvnoise,"vnoise/D",8*1024*1024);
    // a "relative" time
    reil2->Branch("dtime",&reiltime,"dtime/D",8*1024*1024);
    // flucatuation on time
    reil2->Branch("tnoise",&reiltnoise,"tnoise/D",8*1024*1024);
    // Flag for if it was hit or not
    reil2->Branch("anthit",&reilhit,"anthit/I",8*1024*1024);
    // The "ray number" of the hit
    reil2->Branch("pray",&reilray,"pray/I",8*1024*1024);
    reil2->Branch("xhit",&xhit,"xhit/D",8*1024*1024);
    reil2->Branch("yhit",&yhit,"yhit/D",8*1024*1024);
    reil2->Branch("zhit",&zhit,"zhit/D",8*1024*1024);
    
    // the antenna number (which will convey all other antenna info)
    reil2->Branch("antnum",&reilantnum,"antnum/I",8*1024*1024);
    reil2->Branch("antid",&cantid,"antid/I",8*1024*1024);
    reil2->Branch("antpolar",&cantpolar,"antpolar/I",8*1024*1024);
    reil2->Branch("anthit",&anthit,"anthit/I",8*1024*1024);
    reil2->Branch("module",&cmodule,"module/I",8*1024*1024);
    reil2->Branch("string",&cstring,"string/I",8*1024*1024);



  }

  TH1F *hpnu=new TH1F("hpnu","hpnu",5,16.,21.);


  //TH1F *hwhichpols=new TH1F("hwhichpols","hwhichpols",5,0.,5.);

  TProfile *h1=new TProfile("thetanu_vs_posnuz","thetanu_vs_posnuz",50,3.,5.,0.,3.2);

     
  if (DEBUG) cout << "changle(deg)=" << changle*DEGRAD << "\n";

  //Now define the detector center with respect to the geocentric
  //coordinate system.

  mine[0] = 0.;
  mine[1] = 0.;
  mine[2] = R_EARTH - MINE_DEPTH-(2*MAXPOS[2])/2.;


  // loop over possible trigger requirements.
  //	 ftrigeff << "Hits:  " << ihitsrequired << "\n";
  // ftrigeff << "Nsigma:  " << isigmarequired << "\n";
	 
  GetDetector(rxpos,rxtype); // get antenna positions with respect to center of detector.
  // Initialize nrxhit[3]
  for (int imodule=0; imodule<NNODESTOTAL; imodule++) {
    for (int p=0;p<3;p++) {
      nrxhit[p].push_back(0);
    }
  }
	
 
  // clear some things
  for (int p=0;p<3;p++) {
    eventsfound[p]=0;
  }
  eventsfound_below=0;
  nabsorbed_earth=0;
  nabsorbed_noearth=0;
  nweighted=0;
	 
	 
  // Need a neutrino spectrum instead of a unique energy
  // pick a neutrino  momentum
  // pnu=1e19
  //cout << "Enter neutrino momentum (eV).  Example:  10^19 = 1E19.\n";
  //cin >> pnu;
  //cout << "neutrino momentum is " << pnu << "\n";
	 
	 
  // if we are picking from a GZK spectrum, get pnu for each event.


  if (EXPONENT>1) {
    ierr=GetSigma(pnu,sigma,len_int_kgm2);	// for a specific energy, just get the cross section and interaction lengths once.   
  }
	 

  // loop over events
	 
  for (int irun=0;irun<NRUNS;irun++) { // loop over runs
    for (inu=0;inu<NNU;inu++) { // loop over neutrinos



    Distances dist(DISTCUT,NNODEHITSREQUIRED); // create an instance of the Distances vector, for keeping track of the minimum distances between hit antennas.  First input is the cut on the distance between hits nodes.  Second input is the number of nodes required to be hit within this distance.

    //    cout << "starting new event. inu is " << inu << "\n";

    count_neutrinos++; // count neutrinos that we loop over 


    for (int p=0;p<3;p++) { // this counts the number of hits that fire for each ray in each event.
      nrx_fired[p]=0;
      ndipoles_fired[p]=0;
      nslotsx_fired[p]=0;
      nslotsy_fired[p]=0;
    }

    if (EXPONENT==-1) {
      //pnu=TMath::Power(10,16.0+5.0*Rand3.Rndm()*Rand3.Rndm()); 
      pnu=TMath::Power(10,16.0+5.0*Rand3.Rndm()*Rand3.Rndm()); 
      //pnu=TMath::Power(10,16.0+5.0*Rand3.Rndm()*Rand3.Rndm()*Rand3.Rndm()); 
   
    }
    if (EXPONENT==0) {
      pnu=GetGZK();	 
    	 
    }
    
    if (EXPONENT==1) {
      pnu=GetGZKIron();	 
      
    }
    if (EXPONENT==2) {
      pnu=GetE2();	 
      
    }
    if (EXPONENT==3) {
      pnu=GetE2();	 
      
    }

    if (EXPONENT<=3) // if it is a standard energy, then we already calculated these quantities above.
      ierr=GetSigma(pnu,sigma,len_int_kgm2);	// given neutrino momentum, standard deviation and interaction length of neutrino.
    // ierr=0 if the energy is too low for the parameterization
    // ierr=1 otherwise
    //cout << "pnu, len_int are " << pnu << " " << len_int_kgm2/densities[2] << "\n";
    



    if (inu%(int)checkvalue == 0) cout << "inu= " << inu << ", being " << (int)(100*inu/NNU) << "%\n";
    
    if (ierr!=0 || FORSECKEL) {
    
      // This is here to force events to start the same each time
      static int localcount=0;
      if (VERBOSE>=3 || localcount < 3) {
	cout << "FIX SEED FIX SEED" << endl;
	cout << inu << " " << SEED+inu << endl;
	++ localcount;
      }
      Rand3.SetSeed(SEED+inu);
      // get random numbers
    Rand3.RndmArray(10,rndlist);
    if (VERBOSE>=5) {      
      cout << "Random List "; Print(rndlist,10);
    }
    GetNuFlavor(nuflavor);
     

    if (FORSECKEL==1) // For making array of signal vs. freq, viewangle, just use electron neutrinos
      nuflavor="nue";
    if (nuflavor=="nue")  //For outputting to file
      nuflavorint=1;
    else if (nuflavor=="numu")
      nuflavorint=2;
    else if (nuflavor=="nutau")
      nuflavorint=3;



    GetCurrent(current);
      
    // if nu_tau choose tau decay type (never used)
    if (nuflavor=="nutau" && current=="cc") {
      rnd = Rand3.Rndm(1);
      if (rnd<=0.18) 
   	taudecay="m";
      else if(rnd<=0.36) 
   	taudecay="e";
      else if(rnd<=1.0) 
   	taudecay="h";
      else 
   	cout << "error in taudecay";
    }
    
    // ***** muon showers?
    // ***** tau showers?
    // keep track of type of interaction
    y=Gety();

    if (FORSECKEL==1) {
      if (SHOWERTYPE==0) // all had shower
	y=1.;
      if (SHOWERTYPE==1) // all em shower
	y=0.;
    } //if (FORSECKEL)


    if (CONSTANTY==1) { // if we ask to make y a constant=0.2
      y=0.2;
      nuflavor="nue";
      current="cc";
    }


    if (nuflavor=="nue" && current=="cc") {
      emfrac=1.-y;
      hadfrac=y;
    
    }
    else if(nuflavor=="numu" || current=="nc") {
      emfrac=1.E-10;
      hadfrac=y;
    }
    else if(nuflavor=="nutau" && current=="cc") {
      // ***** assume only see first bang?
      // behaves like a muon
      emfrac=1.E-10;
      hadfrac=y;
    }




    // ***** Need to modify?
    // pick a point in the volume that it interacts
    for (int i=0;i<3;i++) {
      posnu_reltominecenter[i]=MAXPOS[i]*(2.*rndlist[i+2]-1.); // relative to center of mine
      if (i==0 || i==1)
	posnu_down[i]=posnu_reltominecenter[i];
      if (i==2)
	posnu_down[i]=-MAXPOS[i]-(posnu_reltominecenter[i]-(-MAXPOS[i])); // posnu_down is the mirror reflection of the interaction point, reflected about the bottom "plane" of the salt dome.
    }

    if (FOLLOWONEEVENT==1 && inu==303) {
      cout << "y is " << y << "\n";
      cout << "nuflavor is " << nuflavor << "\n";
      cout << "emfrac, hadfrac are " << emfrac << " " << hadfrac << "\n";
      cout << "posnu rel. to center of min is ";Print(posnu_reltominecenter,3);
    }

    //    cout << "MAXPOS[2], posnu, posnudown is " << MAXPOS[2] << " " << posnu[2] << " " << posnu_down[2] << "\n";
  

    // for plotting
    posnux2=posnu_down[0];
    posnuy2=posnu_down[1];
    posnuz2=posnu_down[2];
    

    // for plotting
    posnux1=posnu_reltominecenter[0];
    posnuy1=posnu_reltominecenter[1];
    posnuz1=posnu_reltominecenter[2];
    
    //    cout << "I'm here 1.\n";
    

    // pick a neutrino polar angle relative to surface
    // **** looks like it's measured relative to vertical
    //  flat in cos theta
           
    if (USEBELOWHORIZ) 
      costhetanu=2*rndlist[5]-1;
    else
      costhetanu=-rndlist[5];

    // pick a neutrino azimuthal angle
    phinu=TWOPI*rndlist[6];

    // check that these give the right result
    thetanu=acos(costhetanu);

    if (FOLLOWONEEVENT==1 && inu==303) {
      cout << "phinu is " << phinu*DEGRAD << "\n";
      cout << "thetanu is " << thetanu*DEGRAD << "\n";
    }
    sinthetanu=sin(thetanu);

    // find direction vector of neutrino
    // **** are cosine and sine flipped?
    nnu[0]=sinthetanu*cos(phinu);
    nnu[1]=sinthetanu*sin(phinu);
    nnu[2]=costhetanu;
	
    if (VERBOSE>=5) {
	cout << nuflavor << " " << current << " ";
	cout << phinu << " theta " << thetanu << " ";
	Print(nnu,3);
	Print(posnu_reltominecenter,3);
	Print(posnu_down,3);
	cout << "em had y " << emfrac << " " << hadfrac << " " << y << endl;
      }
    // ****  put in energy dependence
    // for now assume we always have an interaction...later could put in
    //    the energy dependence

    // Here we find the six possible lengths of the vector from 
    //the neutrino position to the wall of the detector. Note that
    //there are 2 such possible vectors.


    //    if (dGetTheta(nnu)>1.7 && posnu[2]<-200.)
    //cout << "inu, theta, posnu are " << inu << " " << dGetTheta(nnu) << " " << posnu[2] << "\n";


    if (Getmine(nnu,posnu_reltominecenter, mine_in, mine_out)) { // get the mine entrance and exit points given neutrino direction and interaction point.  Returns zero if the neutrino's path through the mine is less than 1 meter.

      // NeutrinoEvent::NeutrinoEvent(irun,inu, // why long64's
// 				   mine_out, // this needs to be a TVector3
// 		    pnu, 
// 				   AskCons::ParticleType_t flavour,
// 				   AskCons::MaterialType_t matType,
// 				   posnu, // this needs to be a tvector3 
// 				   nnu, // this needs to be a tvector3
// 				   AskCons::InteractionType_t intType,
// 				   y); // this should be a Double_t



      count_getmine++; // count events that pass Getmine

      // get the depth of the interaction
      GetDepth(posnu_reltominecenter,depth_temp);

      //Redefine the position of the neutrino in the geocentric frame.
      //Redefine mine entrance/exit in Earth frame

      for (int i=0;i<3;i++) {
	posnu_reltoearthcenter[i] = mine[i]+posnu_reltominecenter[i]; // relative to center of earth
      }
      
      // convert posnu for flat earth surface to curved earth surface
      EarthCurvature(posnu_reltoearthcenter,depth_temp);
      
      altitude_interaction=-1.*depth_temp; // altitude of interaction (negative number)

      // get the depth of the place where the neutrino enters the mine
      GetDepth(mine_in,depth_temp);

      for (int i=0;i<3;i++) {
	mine_in[i] += mine[i];
      }

      // convert mine_in for flat earth surface to curved earth surface
      EarthCurvature(mine_in,depth_temp);

       // get the depth of the place where the neutrino would exit the mine
      GetDepth(mine_out,depth_temp);

      

      for (int i=0;i<3;i++) {
	mine_out[i] +=mine[i];
      }
      
      // convert mine_out for flat earth surface to curved earth surface
      EarthCurvature(mine_out,depth_temp);

      //Now determine the two vector lengths from the neutrino position to 
      //the Earth's surface.
      // (place where neutrino enters and would exit the earth)
      GetEntryExit(R_EARTH,nnu,posnu_reltoearthcenter,earth_in,earth_out);



      //Now get the chord length and the chord_kgm3 for use in choosing the
      //neutrino interaction point in the next section.

	  
      if (Getchord(earth_in,posnu_reltoearthcenter,chord,chord_kgm2)) { 
	
       
	count_getchord++; // count events that pass Getchord.  They should all pass.

	// chord[3] is the 3-d vector from the earth entrance point to the interaction
	
	// 11/04/03
	// add weighting for chords.  
	// inputs are:  d1,d2,chord,len_int_kgm2 
	// (these variables are set in the code that Dawn wrote)
	// outputs:  weight2,remaining
	    
	    
	d1=sqrt(TMath::Power(mine_in[0]-earth_in[0],2)+TMath::Power(mine_in[1]-earth_in[1],2)+TMath::Power(mine_in[2]-earth_in[2],2)); // from earth entrance point to salt dome entrance point
	d2=sqrt(TMath::Power(posnu_reltoearthcenter[0]-mine_in[0],2)+TMath::Power(posnu_reltoearthcenter[1]-mine_in[1],2)+TMath::Power(posnu_reltoearthcenter[2]-mine_in[2],2)); // from salt dome entrance point to interaction point
	if (d2>chord) {
	  cout<<"error chord"<<" "<<d2<<" "<<chord<<"\n";
	  cout << "inu is " << inu << "\n";
	}  
	    
	chord_kgm2_noearth=d2*RHOMEDIUM; // chord traverse just within dome
	    
	len_int_kgm2=(M_NUCL/sigma); // interaction length is kg/m^2
	IsAbsorbed(chord_kgm2_noearth,len_int_kgm2,weight3); // chance of survival in chord just through salt
	    
	nabsorbed_noearth+=weight3; // number that would survive trip through salt if there were no earth
	    
	if (ABSORBINEARTH) {
	  IsAbsorbed(chord_kgm2,len_int_kgm2,weight1);  //Is it absorbed between the entrance point and the interaction 
	  
	}
	
       

	nabsorbed_earth+=weight1; // number absorbed on its trip to the interaction point
      
	if (d2>chord) {
		
	  for (int i=0;i<3;i++) {
	    cout<<mine_in[i] <<" "<<posnu_reltoearthcenter[i] << " " << mine_out[i]<<" "<<earth_in[i]<< " " << nnu[i] << "\n";
	  }
	  cout << "lengths are " << sqrt(dSquare(mine_in)) << " " << sqrt(dSquare(posnu_reltoearthcenter)) << " " << sqrt(dSquare(mine_out)) << " " << sqrt(dSquare(earth_in)) << " " << sqrt(dSquare(nnu)) << "\n"; 


	}


	if (VERBOSE>=5) {
	  cout << pnu << " a " << d1 << " b " << d2 << " c ";
	  cout << chord_kgm2 << " d ";
	  cout << chord << " e " << len_int_kgm2 << " f ";
	  cout << weight2 << " " << remaining << endl;
	}

	if (!WEIGHTCHORDS || GetChordWeight(pnu,d1,d2,chord_kgm2,chord,len_int_kgm2,weight2,remaining) || FORSECKEL) { // weight2 is probability interaction occurs in salt, of all places along its path.  Returns 0 when the neutrino's path length is greater than 50 interaction lengths.
	  
	  count_getchordweight++; // when WEIGHTCHORDS==1 (default), this is the number of neutrinos that pass GetChordWeight (their chord length is less than 50 interaction lengths)

	  
	  
	  npass_noweight++;
	  //nweighted+= weight2*weight1/weight3;
	  //weight=weight2*weight1/weight3;
	  //  if (inu<50)
//  	    cout << "Messed with weights!  inu is " << inu << "\n";
//  	  if (thetanu>1.7)
//  	    weight1=0.;

	  nweighted+=(1-weight1);
	  weight=(1-weight1);

		
	  // now loop over antennas to find distance to shower and response
		
	  //Redefine the position of the neutrino in the geocentric frame.
		
	  

	  
	  //WhichModule(posnu,MINMODULEZ,MAXMODULEZ);

	  vmmhz0_max=0;
	  vmmhz1_max=0;
	  volts0_max=0;
	  volts1_max=0;
	       
	  increment_atten_factor_chanceinhell=0; // this keeps track of whether or not you have already incremented count_atten_factor_chanceinhell for this event 
	  increment_viewangle_chanceinhell=0;// this keeps track of whether or not you have already incremented count_viewangle_chanceinhell for this event 

	  
	    // We will need the index of refraction at the interaction depth
	  if (DEPTH_DEPENDENT_N) {
	    N_DEPTH=GetNRonneIceShelf(altitude_interaction);
	    changle=acos(1/N_DEPTH);
	  }
	  else 
	    N_DEPTH=NMEDIUM;

	  if (FORSECKEL==1)
	    changle=acos(1/NSALT);



	  //clear out the nhit array
	  for (int p=0;p<3;p++) {
	    nmoduleshit[p]=0;//[ifreq]=0;
	    nrx_fired[p]=0;
	    ndipoles_fired[p]=0;
	    nslotsx_fired[p]=0;
	    nslotsy_fired[p]=0;
	  }
	
	  for (int imodule=0;imodule<NNODESTOTAL;imodule++) {
	    for (int p=0;p<3;p++) {
	      nrxhit[p][imodule]=0;
	    }
	  }

	  //	  cout << "I'm here 2.\n";

	  nrxhit_total=0; // total antennas hit in a event
	  nrxhit_dipole=0; // total dipoles hit in an event
	  nrxhit_x=0;  // total x slots hit in an event
	  nrxhit_y=0; // total y slots hit in an event

	  anyanthit=0; // this checks if any antenna in the event is hit
	  for (int k=0;k<5;k++) {
	    anymodulehit[k]=0;// this checks if any module in the event is hit
	  }
	  for (irx=0;irx<NANTENNA;irx++) {		
	    for (int i=0;i<3;i++) {
	      rx_reltominecenter[i]=rxpos[i][irx]; 
	      //if (i==2)
	      //cout << "rx is " << rx[2] << "\n";
	    }
	    
	    for (int p=0;p<2;p++) {
	      tir[p]=0; // keeps track of whether a ray is totally internally reflected
	      tir_thisantenna=0;
	    }



	    if (VERBOSE>=5) {
	      cout << rxtype[0][irx] << " ";
	      cout << rxtype[1][irx] << " ";
	      cout << rxtype[2][irx] << " ";
	      cout << rxtype[3][irx] << " ";
	      Print(rx_reltominecenter,3);
	    }

	    cantid=rxtype[0][irx]; // antenna orientation
	    //cout << "at the top, cmodule is " << cmodule << "\n";
	    cmodule=rxtype[1][irx]; // antenna string
	    cstring=rxtype[2][irx]; // antenna module
	    cantpolar=rxtype[3][irx];
	
	    xhit=rx_reltominecenter[0];
	    yhit=rx_reltominecenter[1];
	    zhit=rx_reltominecenter[2];

	    GetDepth(rx_reltominecenter,depth_temp);

//  	    if (inu==9) {
//  	      cout << "before:  rx, depth are ";Print(rx,3);
//  	      cout << "depth_temp " << depth_temp << "\n";
//  	    }


	    // we are inside a loop over antennas (irx)

	    for (int i=0;i<3;i++) {
	      rx_reltoearthcenter[i] = rx_reltominecenter[i]+mine[i]; // relative to center of earth
	    }

	    if (inu==0 && irx==896) { 

	      cout << "rx is ";Print(rx_reltoearthcenter,3);
	    }
	
	    EarthCurvature(rx_reltoearthcenter,depth_temp);


	    


//  	    if (inu==9) {
//  	      cout << "after:  rx, depth are ";Print(rx,3);
//  	      cout << "depth_temp " << depth_temp << "\n";
//  	      cout << "posnu is ";Print(posnu,3);
//  	    }

	    
	    //		  if (rx[0]>1000)
	    //cout << "rx is " << rx[0] << "\n";
	    
	    
	    // first get distance
	    

	    if (inu==0 && irx==896) { 

	      cout << "rx is ";Print(rx_reltoearthcenter,3);
	    }
	       	    
	    distance[0]=GetDistance(rx_reltominecenter,posnu_reltominecenter);
	    distance[1]=GetDistance(rx_reltominecenter,posnu_down);
	    

	    //cout << "distances are " << distance[0] << " " << distance[1] << "\n";

	    //if (posnu[2]>4000.)
	    //cout << "rx[2], posnu[2] are " << rx[2] << " " << posnu[2] << "\n";
	    // we are inside a loop over antennas (irx)
	    
	    distance0=distance[0];
	    distance1=distance[1];
	    
	    //		if (inu==18 && distance[0]<2000) {
	    //cout << "posnu is ";Print(posnu,3);
	    //cout << "rx is ";Print(rx,3);
	    //cout << "distance[0] is " << distance[0] << "\n";
	    
	    //cout << "distance[1] is " << distance[1] << "\n";
	    //}
	    // now get direction vector of line from interaction point to antenna
	    if (VERBOSE>=5) {
	      cout << "dist " << distance[0] << " " << distance[1] << endl;
	    }
	    if (distance[0]!=0 && distance[1]!=0) {
	      

	      // the nsurf is upside down because we're using GetRayIceSide backwards
	      for (int i=0;i<3;i++) {
		nsurf[0]=0.;
		nsurf[1]=0.;
		nsurf[2]=-1.;
	      }


	      if (distance[0]==0 || distance[1]==0) {

		cout << "Warning!  About to divide by zero.  distances are " << distance[0] << " " << distance[1] << "\n";
	      }
	      else{
		for (int i=0;i<3;i++) { // for direct rays.  This points from the interaction point to the antenna
		  direction[0][i]=(rx_reltominecenter[i]-posnu_reltominecenter[i])/distance[0];
		}
		for (int i=0;i<3;i++) { // for reflected rays.  This points from the mirror interaction point to the antenna.  It is the direction the ray is going as it arrives at the antenna.
		  direction[1][i]=(rx_reltominecenter[i]-posnu_down[i])/distance[1];
		}
	      }
	      

	      // we are inside a loop over antennas (irx) inside a
	      // if (distance[0]!=0 && distance[1]!=0)

	      costhray0_atinteraction=cos(dGetTheta(direction[0]));
	      costhray1_atinteraction=cos(dGetTheta(direction[1]));

	      vposnu_reltoearthcenter.SetXYZ(posnu_reltoearthcenter[0],posnu_reltoearthcenter[1],posnu_reltoearthcenter[2]);
	      vnsurf.SetXYZ(nsurf[0],nsurf[1],nsurf[2]);

	      

	      // if we're doing depth dependent index of refraction for ARIANNA
	      if (DEPTH_DEPENDENT_N) {
		
		depth_thisantenna=sqrt(dSquare(earth_out))-sqrt(dSquare(rx_reltoearthcenter)); // we're inside a loop over antennas
		
		// loop through rays we're interested in
		// the loop over direct and reflected rays was designed for arianna.
		for (int p=MINRAY;p<=MAXRAY;p++) {
		  
		  // direction is the direct line from the interaction to the antenna
		  for (int i=0;i<3;i++) { // backwards is opposite the direction of the ray in ice
		    backwards[p][i]=-1.*direction[p][i];
		  }
		  
		  // if the interaction has occured inside the firn
		  if (sqrt(dSquare(earth_out))-sqrt(dSquare(posnu_reltoearthcenter))<FIRNDEPTH) {
		    for (int m=0;m<3;m++) {
		      firn_start[m]=posnu_reltoearthcenter[m];// set the firn entrance point to be the interaction point
		    }

		  }
		  else // if the interaction occured at a greater depth then find the point where it enters the firn
		    // getentryexit takes a radius rel. to the earth center, a neutrino direction and an interaction position and gets the place where the neutrino entered and exit at that radius
		    GetEntryExit(R_EARTH-FIRNDEPTH,nnu,posnu_reltoearthcenter,firn_in,firn_start); // firn_start is supposed to be the place where the signal enters the firn.  It looks like I did it wrong and found the place where the neutrino enters the firn.

		  // this is the same as firn_start, but a vector
		  vfirn_start.SetXYZ(firn_start[0],firn_start[1],firn_start[2]);

		  for (int i=0;i<3;i++) {
		    position_bent[p][i]=firn_start[i]; // this keeps track of the position of the bent ray at each step.
		  }

		  vposition_bent[p].SetXYZ(position_bent[p][0],position_bent[p][1],position_bent[p][2]); // this is the same thing as the above but a vector.

		  maxdepth_firn=sqrt(dSquare(earth_out))-sqrt(dSquare(firn_start));// maximum depth of the firn

		  stepsize_depth=maxdepth_firn/NSTEPS_DEPTH;// size of our steps in depth

		  //		  if (!TIR(nsurf,backwards[p],NICE,NFIRN_RONNE)) { // if this ray isn't TIR'ed
		  // the reason why we need to use the backwards ray is that GetRayIceSide calculates the direction in the medium with the higher index of refraction given the direction in the lower index of refraction.
		    // find the ray at the top of the firn
		  n_depth_previous=NICE; // initialize the index of refraction of the previous step to that of the ice

		  for (istep_depth=0;istep_depth<NSTEPS_DEPTH;istep_depth++) { // step in depth and calculate the bending

		    depth_thisstep=maxdepth_firn-stepsize_depth*(double)istep_depth; // find the depth at this step
		    
		    n_depth=GetNSouthPole(depth_thisstep); // get the index of refraction at this depth
		    // find the index of refraction at each step

// 		    cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
// 		    cout << "firn_start is " << firn_start[0] << " " << firn_start[1] << " " << firn_start[2] << "\n";
// 		    cout << "nsurf is " << nsurf[0] << " " << nsurf[1] << " " << nsurf[2] << "\n";
		    if (!TIR(nsurf,backwards[p],n_depth_previous,n_depth) && depth_thisstep>depth_thisantenna) { // if this ray isn't TIR'ed and it's not above the surface
		      
		      
		      costhray_depth_beforerayiceside=cos(dGetTheta(backwards[p])); // this is for plotting

		      // getrayiceside does snell's law
		      // the output is a 3d array backwards_next that 
		      // is in direction of the ray at the next step
		      GetRayIceSide(backwards[p],
				    nsurf,
				    n_depth_previous,
				    n_depth,
				    
				    backwards_next[p]);
		      
		      costhray_depth=cos(dGetTheta(backwards_next[p])); // for plotting

// 		      for (int i=0;i<3;i++) {
// 			// there is a minus sign because backwards is backwards
// 			position_bent[p][i]-=stepsize_depth/dDot(backwards_next[p],nsurf,3)*backwards_next[p][i];
// 		      }
		      
		      vbackwards[p].SetXYZ(backwards_next[p][0],backwards_next[p][1],backwards_next[p][2]);
		      vposition_bent[p]=vposition_bent[p]-stepsize_depth/dDot(backwards_next[p],nsurf,3)*vbackwards[p];



		      
		      // the following is for findin xpos and zpos for my bending plots
		      // need to find what this is in the plane of nsurf and backwards
		      vdiff=vposition_bent[p]-vfirn_start;


		      vperp=(vnsurf.Cross(vbackwards[p])).Unit(); // unit vector perp to both vbackwards and vnsurf
		      vtemp=vperp.Cross(vdiff); // this should be in the plane of nsurf and vbackwards and have the correct length
		      vtemp=vtemp.Cross(vperp);
		      zpos=-1.*vtemp.Dot(vnsurf);
		      xpos=sqrt(vtemp.Mag()*vtemp.Mag()-zpos*zpos);

// 		    cout << "angle between backwards and backwards_next is " << dDot(backwards[p],backwards_next[p],3) << "\n";
// 		    cout << "angle between backwards_next and nsurf is " << dDot(backwards_next[p],nsurf,3) << "\n";
		      if (tree_steps->GetEntries()<HIST_MAX_ENTRIES*100 && HIST==1)
			tree_steps->Fill();
		    }
		    else if (TIR(nsurf,backwards[p],n_depth_previous,n_depth)) {
		      istep_depth=NSTEPS_DEPTH; // if it would be tir'ed, stop looping.
		      tir[p]=1;
		      tir_thisantenna=1;

		    }
		    else if (depth_thisantenna>depth_thisstep) {
		      istep_depth=NSTEPS_DEPTH;
		    }

		    for (int i=0;i<3;i++) {
		      backwards[p][i]=backwards_next[p][i];
		    }

		    n_depth_previous=n_depth;

		  } // end stepping over depths from the bottom of the firn to the top

		  for (int i=0;i<3;i++) {
		    // flip it back around so it is in the direction the ray travels
		    direction[p][i]=-1.*backwards_next[p][i]; 
		  }
		  


// 		  } // end if not tir

	       		  
		} // end loop over rays
		
		
	      } // if we're using depth dependent index of refraction
	      

	      // we are inside a loop over antennas (irx) inside a
	      // if (distance[0]!=0 && distance[1]!=0)


	      direction00=direction[0][0];
	      direction01=direction[0][1];
	      direction02=direction[0][2];
	      direction10=direction[1][0];
	      direction11=direction[1][1];
	      direction12=direction[1][2];
	      
	      
	      // find the direction of the ray before it reflects from the bottom
	      downwardray[0]=rx_reltominecenter[0]-posnu_reltominecenter[0];
	      downwardray[1]=rx_reltominecenter[1]-posnu_reltominecenter[1];
	      downwardray[2]=posnu_down[2]-rx_reltominecenter[2];
	      
	      // normalize it.
  	      templength=sqrt(dSquare(downwardray));
  	      for (int i=0;i<3;i++) {
  		downwardray[i]=downwardray[i]/templength;
  	      }

	      // we are inside a loop over antennas (irx) inside a
	      // if (distance[0]!=0 && distance[1]!=0)

	      phiray[0]=dGetPhi(direction[0]);
	      costhray[0]=cos(dGetTheta(direction[0]));
	      
      	      phiray[1]=dGetPhi(direction[1]);
	      costhray[1]=cos(dGetTheta(direction[1]));
	      costhray0=costhray[0]; // this is the direction as it hits the antenna
	      costhray1=costhray[1]; 
	    } // end if distance != 0
	    else {
	      cout << "danger distance=0!  crapping out!\n";
	      return 0;
	    }
	    
	    
	    atten_factor[0]=GetFac2(rx_reltominecenter,posnu_reltominecenter,attenlength);       
	    atten_factor[1]=GetFac2(rx_reltominecenter,posnu_down,attenlength);
	    atten_factor0=atten_factor[0];  // for inserting into a tree
	    atten_factor1=atten_factor[1]; // for inserting into a tree.

	    double fac1[2];
	    fac1[0]=GetFac1(rx_reltominecenter,posnu_reltominecenter);
	    fac1[1]=GetFac1(rx_reltominecenter,posnu_down);
	  

	    // we are inside a loop over antennas (irx)

	    viewangle[0]=GetViewAngle(direction[0],nnu);
	    viewangle0=viewangle[0];  // for inserting into tree
	    viewangle[1]=GetViewAngle(downwardray,nnu);
	    viewangle1=viewangle[1];  // for inserting into tree
	    
	    double fac3[2];
	    //cout << viewangle0 << " " << changle << " "<< pnu << " ";
	    //cout << emfrac << " " << hadfrac << " " << FREQ_LOW << endl;
	    fac3[0]=GetFac3(viewangle0, changle, pnu, emfrac, hadfrac, 
			    freqlist[0],deltheta_em_max,deltheta_had_max);
	    fac3[1]=GetFac3(viewangle1, changle, pnu, emfrac, hadfrac, 
			    freqlist[0],deltheta_em_max,deltheta_had_max);

	    


	    
	    for (int p=0;p<2;p++) {
	      volts[p]=0;
	    }
	    fluctime=0;
	    
	    if (VERBOSE>=5) {
	      cout << pnu << " " << emfrac << " " << hadfrac << " ";
	      cout << changle << " " << viewangle0 << " " << viewangle1 << " ";
	      cout << "Fac " << fac3[0] << " " << fac3[1] << endl;
	    }

	    //	    cout << "I'm here 3.\n";

	    // we are inside a loop over antennas (irx)

	    if (Max(fac3[0],fac3[1])>0) { // fac3's are greater than 0 when the viewing angle is less than 20 cone widths away from the cherenkov angle
	      if (increment_viewangle_chanceinhell==0) {// only want to increment this once per antenna
		count_viewangle_chanceinhell++; // counting how many antennas pass this cut
		increment_viewangle_chanceinhell++; // only want to do this once per event so this switch tells us that we have already incremented
	      }
	      if (atten_factor[0]!=0
		  || atten_factor[1]!=0) { // these are greater than 0 when you are less than 10 attenuation lengths away 
		
		if (increment_atten_factor_chanceinhell==0) { // only want to increment this once per antenna
		  count_atten_factor_chanceinhell++; // counting how many antennas pass this cut
		  increment_atten_factor_chanceinhell++; // only want to do th is once per antenna so this switch tells us that we have already incremented
		}	



	    // we are inside a loop over antennas (irx)
		// and  if (Max(fac3[0],fac3[1])>0) { 
		// and  if (atten_factor[0]!=0
		//		  || atten_factor[1]!=0) {


		for (int ifreq=0;ifreq<NFREQ;ifreq++) { // loop through frequencies
		  thisfreq=freqlist[ifreq];
		  
		  // calculate e-field/MHz at this wavelength at 1m.
		  
		  //		  cout << "pnu is " << pnu << "\n";
		  
		  // loop through rays we're interested in
		  for (whichray=MINRAY;whichray<=MAXRAY;whichray++) {

		    //for (int p=0;p<1;p++) {

		    // if this ray isn't tir'ed
		    if (!tir[whichray]) {
		    vmmhz1m[whichray]=GetVmMHz1m(pnu,thisfreq,X0MEDIUM,ECMEDIUM,N_DEPTH,AEXMEDIUM,WHICHPARAMETERIZATION);
		    vmmhz1m_perfreq[ifreq]=vmmhz1m[whichray];


	    // we are inside a loop over antennas (irx)
		// and  if (Max(fac3[0],fac3[1])>0) { 
		// and  if (atten_factor[0]!=0
		    //		  || atten_factor[1]!=0) {
		    // and for (int ifreq=0;ifreq<NFREQ;ifreq++) {
		    // and for (whichray=MINRAY;whichray<=MAXRAY;whichray++) 
		    // and if (!tir[whichray]) {
		    
		    
		    
		    if (VERBOSE>=5) {
		      cout << " 1m " << pnu << " " << whichray << " " << vmmhz1m[whichray] << endl;
		    }
		    // Note: Comments about GetVmMHz are in shared.cc at EOF
		    if (FOLLOWONEEVENT==1 && inu==303 && ifreq==0 && (irx==3888 || irx==3889 || irx==3891)) {
		      cout << "thisfreq, irx, vmmhz1m[whichray] are " << thisfreq << " " << irx << " " << vmmhz1m[whichray] << "\n";
		    }
		    

		    fac3[whichray]=GetFac3(viewangle0, changle, pnu, emfrac, hadfrac, 
				    thisfreq,deltheta_em[ifreq],deltheta_had[ifreq]);


	    // we are inside a loop over antennas (irx)
		// and  if (Max(fac3[0],fac3[1])>0) { 
		// and  if (atten_factor[0]!=0
		    //		  || atten_factor[1]!=0) {
		    // and for (int ifreq=0;ifreq<NFREQ;ifreq++) {
		    // and for (whichray=MINRAY;whichray<=MAXRAY;whichray++) 
		    // and if (!tir[whichray]) {



		    if (FORSECKEL==1) {
		      for (int iviewangle=0;iviewangle<NVIEWANGLE;iviewangle++) {
			
			vmmhz_temp=vmmhz1m_perfreq[ifreq];
			
			
			TaperVmMHz(viewangles[iviewangle],deltheta_em[ifreq],deltheta_had[ifreq],emfrac,hadfrac,WHICHPARAMETERIZATION,
				   vmmhz_temp);
			
			forseckel[iviewangle][ifreq]=vmmhz_temp;
			//if (k==0 && iviewangle==0)
			//cout << "viewangle, deltheta_em[k], deltheta_had[k], emfrac, hadfrac, vmmhz_temp are " << viewangle << " " << deltheta_em[k] << " " << deltheta_had[k] << " " << emfrac << " " << hadfrac << " " << vmmhz_temp << "\n";
		      } //for (loop over viewing angles)
		    } //if (FORSECKEL==1)


		    // we are inside a loop over antennas (irx)
		    // and  if (Max(fac3[0],fac3[1])>0) { 
		    // and  if (atten_factor[0]!=0
		    //		  || atten_factor[1]!=0) {
		    // and for (int ifreq=0;ifreq<NFREQ;ifreq++) {
		    // and for (whichray=MINRAY;whichray<=MAXRAY;whichray++) 
		    // and if (!tir[whichray]) {



		    // sum
		    vmmhz1m[whichray]*=fac1[whichray]*fac3[whichray];
		    vmmhz[whichray]=vmmhz1m[whichray];
		  
		    // account for attenuation in salt
		    vmmhz[whichray]=vmmhz[whichray]*atten_factor[whichray];
		    
		    
		    //if (TYPEOFANTENNAS==0 && (cmodule%2==1))
		    //cantpolar*=-1; // for the dipoles and slots, for every other x slot, for example, turn it to be polarized in the -x direction.  same for y.




		    polarfactor=GetFac4(direction[whichray],
					hitangle_e,hitangle_h,
					e_component,h_component,
					nnu,cantid,cantpolar,
					thisfreq, ifreq,
					costh_pol[whichray], costh_dir[whichray], heff);




	    // we are inside a loop over antennas (irx)
		// and  if (Max(fac3[0],fac3[1])>0) { 
		// and  if (atten_factor[0]!=0
		    //		  || atten_factor[1]!=0) {
		    // and for (int ifreq=0;ifreq<NFREQ;ifreq++) {
		    // and for (whichray=MINRAY;whichray<=MAXRAY;whichray++) 
		    // and if (!tir[whichray]) {



		    if (VERBOSE>=5) {
		      cout << setprecision(3) << "1234 " << fac1[whichray] << " " << atten_factor[whichray] << " " << fac3[whichray] << " " << polarfactor << endl;
		      cout << "heff BW NFREQ" << heff << " " << BW << " " << polarfactor << endl;
		    }
		    volts[whichray]+=vmmhz[whichray]*0.5*heff*((BW/1E6)/NFREQ)*polarfactor;

		    //KAR cout << irx << " volts " << volts[whichray] << endl;
		    if (VERBOSE>=5) {
		      cout << volts[whichray] << " " << whichray << endl;
		    }
		    
		    } // end if not tir
		   //else
		      // cout << "tir is " << tir[whichray] << "\n";
		  } // end loop over rays	
		} // end loop over freq...

		//		cout << "I'm here 4.\n";

	    // we are inside a loop over antennas (irx)
		// and  if (Max(fac3[0],fac3[1])>0) { 
		// and  if (atten_factor[0]!=0
		    //		  || atten_factor[1]!=0) {

		if (FOLLOWONEEVENT==1 && inu==303 && (irx==3888 || irx==3889 || irx==3891)) {
		  cout << "volts is " << volts[0] << "\n";
		}
	      } // end if attenuation factor>0
	    } // end if it's too far off the cerenkov cone
		




	    if (FOLLOWONEEVENT==1 && inu==303 && (irx==3888 || irx==3889 || irx==3891)) {
	      cout << "vnoise_dipolesslots is " << VNOISE_DIPOLESSLOTS << "\n";
	    }
	
			   



	    // we are inside a loop over antennas (irx)

	    AddNoise(cantid,volts);


	    //if (inu==0)
	    //cout << "warning! took out fluctuations.\n";  
	    
	    //		cout << "inu, volts are " << inu << " " << volts << "\n";
	    
	    //  		    if (inu<=1000) {
	    
	    //  		      tree1->Fill();
	    //  		    }
	    
	    // ******* probably also want to record where these hits occur
	    // ******* get position resolution
	    // ******* get angular resolution
	    // find out how many antennas have s/n >6

	    // we are inside a loop over antennas (irx)

	    
	    counting++;
	    anthit=0;  // set to zero the number of antennas hit
	  	    
	    
	    for (int p=MINRAY;p<=MAXRAY;p++) { //p=0 is direct ray
	      //p=1 is reflected ray

	      if (TYPEOFANTENNAS==1 && (cantpolar==0 || cantpolar==1))
		vnoise=VNOISE_SEAVEYS;
	      else
		vnoise=VNOISE_DIPOLESSLOTS;

	      //	      cout << "cmodule is " << cmodule << "\n";
	      if (fabs(volts[p])/vnoise>=isigmarequired) {


		if (VERBOSE>=4) {
		  cout << "Trigger irx= " << irx << " " << inu << endl;
		}

	    // we are inside a loop over antennas (irx)
		// and for (int p=MINRAY;p<=MAXRAY;p++) {
		// and if (volts[p]/vnoise>=isigmarequired)


		if (FOLLOWONEEVENT==1 && inu==303) {
		  if (1==1) {
		    for (int ti=0; ti<3; ++ti)
		      cout << rx_reltominecenter[ti] << " ";
		    cout << endl;
		    for (int ti=0; ti<3; ++ti)
		      cout << posnu_reltominecenter[ti] << " ";
		    cout << endl;
		    for (int ti=0; ti<3; ++ti)
		      cout << posnu_reltominecenter[ti] - rx_reltominecenter[ti] << " ";
		    cout << endl;
		    
		    cout << viewangle0 << " " << changle  << " " << pnu  << " " << emfrac  << " " << hadfrac  << " " << 200.0E6 << endl;
		    cout << vmmhz1m[0] << " ";
		    cout << "(" << direction[0][0] << "," << direction[0][1] << "," << direction[0][2] << ") ";
		    cout << "(" << nnu[0] << "," << nnu[1] << "," << nnu[2] << ") ";
		  
		    cout << fac1[0] << " " << atten_factor[0] << " ";
		    cout << fac3[0] << " ";
		    cout << polarfactor << " " << endl;
		    cout << emfrac << " " << hadfrac << " heff" << heff << endl;
		    cout << "p " << p << " " << volts[p] << endl;
		  }
		  cout << "hit irx is " << irx << "\n";
		} // end FOLLOWONEEVENT loop


	    // we are inside a loop over antennas (irx)
		// and for (int p=MINRAY;p<=MAXRAY;p++) {
		// and if (volts[p]/vnoise>=isigmarequired)

	      
		anthit++; 

		if (anyanthit==0) {
		  count_anthit++; // count events where at least one antenna is hit
		  anyanthit=1;
		}

		nrx_fired[p]++; // number of antennas fired for one event
		nrx_fired_total[p]++; // total number of antennas fired

		if (TYPEOFANTENNAS==0) {
		  if (cantpolar==0) {
		    ndipoles_fired[p]++; // number of dipoles fired for one event
		    ndipoles_fired_total[p]++; // total
		  }
		  if (abs(cantpolar)==1) {
		    nslotsx_fired[p]++;// # slots in x direction fired for one event
		    nslotsx_fired_total[p]++; // total
		  }
		  if (abs(cantpolar)==2) {
		    nslotsy_fired[p]++;// # slots in y direction fired for one event
		    nslotsy_fired_total[p]++; // total
		    
		  }
		}
		else if (TYPEOFANTENNAS==1) {

		  if (cantpolar==0) {
		  
		    nhornsx_fired_total[p]++; // total
		  }
		  if (abs(cantpolar)==1) {
	       
		    nhornsy_fired_total[p]++; // total
		  }
		  if (abs(cantpolar)==2) {
		
		    ndipoles_fired_total[p]++; // total
		    
		  }

		}


	    // we are inside a loop over antennas (irx)
		// and for (int p=MINRAY;p<=MAXRAY;p++) {
		// and if (volts[p]/vnoise>=isigmarequired)

		
		//		cout << "I'm here 4.5\n";
		nrxhit[p][cmodule]++; // records how many rx hit on each module by each ray
		//cout << "I'm here 4.51\n";
		for (int k=0;k<5;k++) {
		  if (nrxhit[p][cmodule]==k+1 && anymodulehit[k]==0 && p==0) {
		    count_modulehit[k]++; // count events where at least one module is hit
		    anymodulehit[k]=1;
		  }
		}

		//		cout << "I'm here 4.7\n";

		if (anthit<2) { // if the antenna hasn't already been hit by another ray 
		  nrxhit[2][cmodule]++; // record that an antenna on this module has been hit
		  nrx_fired[2]++; // the 2nd element of all of these array is the union of reflected and direct hits
		  ndipoles_fired[2]++;
		  nslotsx_fired[2]++;
		  nslotsy_fired[2]++;

		  

		  nrx_fired_total[2]++;
		  ndipoles_fired_total[2]++;
		  nslotsx_fired_total[2]++;
		  nslotsy_fired_total[2]++;


		  nhornsx_fired_total[2]++;
		  nhornsy_fired_total[2]++;

		}
		

		//		cout << "I'm here 4.8\n";
	    // we are inside a loop over antennas (irx)
		// and for (int p=MINRAY;p<=MAXRAY;p++) {
		// and if (volts[p]/vnoise>=isigmarequired)



		ray_thishit=p; // which ray causes this antenna to be hit
		volts0=volts[0];
		volts1=volts[1];
		vmmhz0=vmmhz[0];
		vmmhz1=vmmhz[1];
		dtime=distance[p]/(CLIGHT/NMEDIUM);
		Double_t tnoise=100E-9;          
		fluctime=Rand3.Gaus(0,tnoise);    
		dtime+=fluctime;                  

		if (HIST) {
		  reilinu=inu;
		  reilray=p;
		  reilantnum=irx;
		  reilhit=1;
		  reilvolts=volts[p];
		  reilvnoise=VNOISE_DIPOLESSLOTS; // better to have actual fluc
		  reiltime=dtime;
		  reiltnoise=fluctime;
		  //		  cout << "I'm here 4.88\n";
		  if (reil2->GetEntries()<HIST_MAX_ENTRIES && HIST==1)
		    reil2->Fill();
		}
		//		cout << "I'm here 4.89\n";
		if (tree_hitantennas->GetEntries()<HIST_MAX_ENTRIES && HIST==1) {
		  tree_hitantennas->Fill();

		}
	      } // end if this antenna passes the trigger
	    } // end loop over rays
		
	    //	    cout << "I'm here 4.9\n";

	    // we are inside a loop over antennas (irx)


	    vmmhz0=vmmhz[0];
	    vmmhz1=vmmhz[1];
	    
	    volts0=volts[0];
	    volts1=volts[1];
	    
	    if (vmmhz0>vmmhz0_max)
	      vmmhz0_max=vmmhz0;
	    if (vmmhz1>vmmhz1_max)
	      vmmhz1_max=vmmhz1;
	    
	    if (volts0>volts0_max)
	      volts0_max=volts0;
	    if (volts1>volts1_max)
	      volts1_max=volts1;		
	    
	    if (tree_allantennas->GetEntries()<HIST_MAX_ENTRIES && HIST==1)
	      tree_allantennas->Fill();
	    
	  } // antennas

	  //	  cout << "I'm here 5.\n";

	  //	  cout << "event " << inu << "\n";
	  //cout << "NNODES, MINRAY, MAXRAY are " << NNODES << " " << MINRAY << " " << MAXRAY << "\n";
	  int thismodulehit=0; // keeps track of whether this module is hit
	  for (Int_t imodule=0;imodule<NNODESTOTAL;imodule++) {	 
	    thismodulehit=0; // initialize
	    for (int p=MINRAY;p<=MAXRAY;p++) {
	      //for (Int_t ifreq=0;ifreq<NFREQ;ifreq++) {		    
	      //  cout << "nrxhit, irxhitsrequired are " << nrxhit[p][imodule] << " " << irxhitsrequired << "\n";
	      if (nrxhit[p][imodule]>0)		
		if (nrxhit[p][imodule]>=irxhitsrequired) {
		  if (VERBOSE>=2) {
		    cout << "Module " << imodule << " had " << nrxhit[p][imodule];
		    cout << " hits." << endl;		  
		  }
		  nmoduleshit[p]++; // records number of modules hit for one event
		  nmoduleshit_total[p]++; // records number of modules hit for one event
	       
		
		  for (int icomp=0;icomp<3;icomp++) {
		    rx_reltominecenter[icomp]=rxpos[icomp][imodule*NEACHNODE];
		  }

		  dist.AddHit(rx_reltominecenter,p); // add this node position to the vector containing all of the nodes that are hit
		  
		
		  thismodulehit++; // flag this module as hit
		} // the module passes
	    }// end loop over rays


	    if (nrxhit[2][imodule]>=irxhitsrequired) {
	      dist.AddHit(rx_reltominecenter,2);
	      nmoduleshit[2]++; //
	      nmoduleshit_total[2]++; // records number of modules hit for one event. The 2nd element is the union of direct and reflected hits.
	    }
	  
	    // we are inside for (Int_t imodule=0;imodule<NMODULE;imodule++)
	  
	    for (int p=0;p<3;p++) {
	      nhittotal[p]=0;
	      nhittotal[p]+=nmoduleshit[p]; // number of modules hit, total
	    }
	    
	  } // end looping through modules
	  nmoduleshit_fortree=nmoduleshit[2];


	  nrx_fired0=nrx_fired[0];
	  nrx_fired1=nrx_fired[1];
	  
		  
	  Bool_t hitfound=kFALSE;
	  nrays_detected=0;
	  for (int p=0;p<3;p++) { // loop over rays.  p=2 keeps track of antennas/modules that were hit by either the direct or reflected ray.

	    for (int k=0;k<5;k++) {
	      if (nhittotal[p]>=k+1 && p==2) {
		count_passestrigger[k]++; // count events where at least k+1 modules are hit
	      }
	    }

	    // we are inside for (int p=0;p<3;p++)


	    if (dist.DoesItPass(p)) { // see if the hits from this ray pass the trigger


	      //	    if (nhittotal[p]>=NNODEHITSREQUIRED) {
	      if (VERBOSE>=3) {
		cout << "Trigger inu= " << inu << endl;
	      }
	      
	      hitfound=kTRUE;
	      if (p==2)
		count_passing_events++;  // increment this once per event
	      whichray=p;
	      nrays_detected++;
	      eventsfound[p]+=(1-weight1);
	      
	      
	      if (p==1 && nrays_detected==2) {
		doublesfound+=(1-weight1);
		whichray=2;
		
	      }
	      if (p==2) {
		//	    logweight2=log10(weight2*weight1/weight3);
		logweight2=log10(1-weight1);
		
		// for calculating errors on sensitivity
		// need to find how many events as a function of weight
		// here, we find how to index weight2
		if (logweight2<MIN_LOGWEIGHT)  // underflows, set to 0th bin
		  index_weights=0;
		else if (logweight2>MAX_LOGWEIGHT) // overflows, set to last bin
		  index_weights=NBINS-1;
		else
		  index_weights=(int)(((logweight2-MIN_LOGWEIGHT)/(MAX_LOGWEIGHT-MIN_LOGWEIGHT))*(double)NBINS);  // which index weight2 corresponds to.
		


	    // we are inside for (int p=0;p<3;p++)
		// and if (nhittotal[p]>=NNODEHITSREQUIRED)
		// and if (p==2)

		// count number of events that pass, binned in weight
		if (index_weights<NBINS)
		  eventsfound_binned[index_weights]++;
		
		// increment for each flavor
		if (nuflavor=="nue") {
		  //	      sum[0]+=weight2*weight1/weight3;
		  sum[0]+=(1-weight1);
		  eventsfound_binned_e[index_weights]++;
		}
		if (nuflavor=="numu") {
		  //	      sum[1]+=weight2*weight1/weight3;
		  sum[1]+=(1-weight1);
		  eventsfound_binned_mu[index_weights]++;
		}
		if (nuflavor=="nutau") {
		  //	      sum[2]+=weight2*weight1/weight3;
		  sum[2]+=(1-weight1);
		  eventsfound_binned_tau[index_weights]++;
		}
	      

	    // we are inside for (int p=0;p<3;p++)
		// and if (nhittotal[p]>=NNODEHITSREQUIRED)
		// and if (p==2)

		
		
		if (thetanu>90.*RADDEG)
		  //eventsfound_below+=weight2*weight1/weight3;
		  eventsfound_below+=(1-weight1);		
		
		if (HIST==1 && tree_passing_events->GetEntries()<HIST_MAX_ENTRIES) {
		  tree_passing_events->Fill();
		  
		  		  
		  h1->Fill(posnu_reltominecenter[2]/1000.,thetanu);
		  hpnu->Fill(log10(pnu));
		  
		}
	      } // if it passes by either ray, or combination of both
	    }  // passes the trigger

	    maxnnodehitsincluster=dist.vnnodes_close[0]; // grab the maximum number of hits in a cluster for this event for filling the tree
	    // want to fill this whether it passes the trigger or not
	    if (tree_beforetrigger->GetEntries()<HIST_MAX_ENTRIES && HIST==1)
	      tree_beforetrigger->Fill();



	    // we are inside for (int p=0;p<3;p++)

	    if (FOLLOWONEEVENT==1 && inu==303) { 
	      cout << "Number of nodes hit is " << nmoduleshit[p] << "\n";
	      cout << "Number of antennas hit is " << nrx_fired[p] << "\n";
	      cout << "Number of dipoles hit is " << ndipoles_fired[p] << "\n";
	      cout << "Number of x slots hit is " << nslotsx_fired[p] << "\n";
	      cout << "Number of y slots hit is " << nslotsy_fired[p] << "\n";
	    }
  
	    if (hitfound) {
	      reilinu=inu;
	      reilray=p;
	      reilpnu=pnu;
	      if (p==0 || p==2) {
		reilnux=posnu_reltominecenter[0];
		reilnuy=posnu_reltominecenter[1];
		reilnuz=posnu_reltominecenter[2];
	      } else {
		reilnux=posnu_down[0];
		reilnuy=posnu_down[1];
		reilnuz=posnu_down[2];
	      }
	      reilem=emfrac;
	      reilhad=hadfrac;
	      reilnnu0=nnu[0];
	      reilnnu1=nnu[1];
	      reilnnu2=nnu[2];
	      if (reil1->GetEntries()<HIST_MAX_ENTRIES)
		reil1->Fill();	    
	    }	    
	  } // end loop over rays
	}  // end GetChordWeight 
		
	//} // if the neutrino is absorbed
      } //if good chord, chord> 1m
	
    } // end if Getmine
        
    //    if (anyanthit>0)
      //      eventsfound_onehit+=weight2*weight1/weight3;
    //eventsfound+=(1-weight1);

    }// end if neutrino energy within the region of parameterization of the cross section    



  }    // neutrinos
  } // end loop over runs
    


  poyntingflux_1m=0; // sum the max poynting flux at 1 m from the interaction.

  for (int iviewangle=0;iviewangle<NVIEWANGLE;iviewangle++) {
    for (int k=0;k<NFREQ;k++) {
      foutseckel << viewangles[iviewangle] << " " << freqlist[k] << " " << forseckel[iviewangle][k] << "\n";
      if (iviewangle==0)
	poyntingflux_1m+=forseckel[iviewangle][k]*forseckel[iviewangle][k];
    }
  }
  poyntingflux_1m*=BW/(double)NFREQ/1.E6;
  
  g1=new TGraphErrors(NFREQ,freqlist,vmmhz1m_perfreq,freqerr,yerr); // graph showing V/m/MHz

  for (int ifreq=0;ifreq<NFREQ;ifreq++) {
    fvmmhz1m << freqlist[ifreq] << "\t" << vmmhz1m_perfreq[ifreq] << "\n"; 
  }
TCanvas *c1 = new TCanvas("offaxis","offaxis",500,500);

 TH2F *htemp=new TH2F("htemp","htemp",100,0.,1500000000.,100,0.,0.002);
 htemp->Draw();
 
 g1->SetMarkerStyle(1);
 g1->SetMarkerSize(0.5);
 g1->Draw("same");

c1->Print("vmmhz1m.eps");

  

  



  Summarize(pnu,eventsfound,doublesfound,nabsorbed_noearth,nabsorbed_earth,nweighted,sigma);

  // loop through rays
//    for (int p=0;p<3;p++) {
//      cout << "Number of nodes hit is " << nmoduleshit_total[p] << "\n";
//      cout << "Number of antennas hit is " << nrx_fired_total[p] << "\n";
//      cout << "Number of dipoles hit is " << ndipoles_fired_total[p] << "\n";
//      cout << "Number of x slots hit is " << nslotsx_fired_total[p] << "\n";
//      cout << "Number of y slots hit is " << nslotsy_fired_total[p] << "\n";
//      cout << "Number of x horns hit is " << nhornsx_fired_total[p] << "\n";
//      cout << "Number of y horns hit is " << nhornsy_fired_total[p] << "\n\n";

//    }


  //cout << "Doubles found:  " << doublesfound << "\n";
       
      
  fout.close();
  if (HIST==1)
    CloseTFile(hfile);
  cout << __LINE__ << " " << __FILE__ << endl;
  return 0;

}

void AddNoise(int cantid,double *volts) {

  double fluc;
	
  if (SIGNAL_FLUCT==1) { // only do this if it's requested in the input file.
    for (int p=MINRAY;p<=MAXRAY;p++) { // loop through rays
      if (cantid==0) // if we're simulating dipoles and slots
	fluc=Rand3.Gaus(0,VNOISE_DIPOLESSLOTS);
      else if (cantid==1 || cantid==3) // if we're simulating seaveys or log periodics
	fluc=Rand3.Gaus(0,VNOISE_SEAVEYS); // add fluctuation appropriate for their bandwidth
      
      volts[p]+=fluc; // add noise to the signal
    }
  } 
  

}

//int WhichModule(double* posnu,int& MINMODULEZ,int& MAXMODULEZ) {



//double spacing=(2*MAXARRAY[2]-0.1)/(double)NNODES;
//double posnuz_reltobottomofdetector=posnu[2]-(-1*MAXARRAY[2]+MAXPOS[2]-INSTRUMENTEDREGION_DEPTH+MINE_DEPTH);
  
//int irx;
//if (posnuz_reltobottomofdetector<0)
//irx=0;
//else if (posnuz_reltobottomofdetector>2*MAXARRAY[2])
//irx=NNODES;
//else
//irx=(int)(posnuz_reltobottomofdetector/spacing);
  

//cout << "spacing is " << spacing << "\n";
// if it's below the bottom of the detector or it's within 10 attenuation lengths of the bottom
//if (posnuz_reltobottomofdetector<0 || (double)irx*spacing/attenlength<10)
//MINMODULEZ=0;
// otherwise calculate the minimum - 10 attenuation lengths below the position
//else
//MINMODULEZ=irx-(int)(attenlength*10./spacing);

  // if the interaction position is within 10 attenuation lengths of the top of the detector or anywhere above the top of the detector
//if ((posnuz_reltobottomofdetector>0 && (double)(irx-(NNODES-1))*spacing/attenlength<10) || posnuz_reltobottomofdetector>2*MAXARRAY[2])
//MAXMODULEZ=NNODES-1;
//else if (posnuz_reltobottomofdetector<0 && posnuz_reltobottomofdetector/attenlength<-10.)
//MAXMODULEZ=0;
//else {
//MAXMODULEZ=Min(NNODES-1,(int)(posnuz_reltobottomofdetector+attenlength*10.));
//}

//cout << "posnuz_reltobottomofdetector, irx, MINMODULEZ, MAXMODULEZ are " << posnuz_reltobottomofdetector << " " << irx << " " <<  MINMODULEZ << " " << MAXMODULEZ << "\n";
//return irx;
//}

int Min(int x,int y) {
  if (x<y) return x;
  else return y;
}

double Max(double x,double y) {
  if (x>y) return x;
  else return y;
}

void GetDetector(vector<double> rxpos[3], vector<int> rxtype[4]) {
  int writeout=1;
  if (DETECTOR==1) {
    NRXX=2;
    NRXY=2;
  }
  if (DETECTOR==2) {
    NRXX=1;
    NRXY=1;
  }
  

  int NNODESTOTAL=0;

   ofstream fout2("inputs.geom");
    if (writeout!=1) fout2.close();
  if (SHAPE==0) { // rectangular shape
 
    // set up the antennas

    for (int i=0;i<NRXX;i++) { // loop over antennas in the x dimension
   
      for (int j=0;j<NRXY;j++) { // loop over antennas in the y dimension
  

	for (int k=0;k<NNODESEACHSTRING;k++) { // loop over nodes along z
	   
	  for (int n=0;n<NEACHNODE;n++) { // number of antennas in each node
	    // rxpos will be converted to mine-centric coordinates later
	    // these coordinates are detector-centric

	    // position the antennas in x,y
	    //rxpos[0].push_back(-MAXARRAY[0]+(2.*MAXARRAY[0]/(NRXX-1))*((double)i+0.5));
	    //rxpos[1].push_back(-MAXARRAY[1]+(2.*MAXARRAY[1]/(NRXY-1))*((double)j+0.5));
	    rxpos[0].push_back(-MAXARRAY[0]+(2.*MAXARRAY[0]/NRXX)*((double)i+0.5));
	    rxpos[1].push_back(-MAXARRAY[1]+(2.*MAXARRAY[1]/NRXY)*((double)j+0.5));

	    // next we position the nodes in z
	    // for z>1, the position of the highest node will be at +MAXARRAY[2]
	    // for z=1, we position the node in the center of the mine.
	    if (NNODESEACHSTRING>1)
	      rxpos[2].push_back(-MAXARRAY[2] 
				 // begin at bottom of array
				 +(2.*MAXARRAY[2]/(NNODESEACHSTRING-1))*((double)k) 
				 // to z position of module
				 -(SPACINGWITHINMODULE*(double)NEACHNODE)/2. 
				 // to bottom of module
				 +SPACINGWITHINMODULE*((double)n+0.5));
	    else
	      rxpos[2].push_back(-(SPACINGWITHINMODULE*(double)NEACHNODE)/2. 
				 // to bottom of module
				 +SPACINGWITHINMODULE*((double)n+0.5));


	    rxtype[0].push_back(TYPEOFANTENNAS); //cantid
	    rxtype[1].push_back(NNODESTOTAL);        //cmodule
	
	    rxtype[2].push_back(NSTRING);        //cstring

	    int cantpolar=-50;
	    if (TYPEOFANTENNAS==0) {
//  	      if (n%2==0) { 
//  		cantpolar=0;  //0,2,4,6,8,10
//  	      } 
//  	      else {
//  		if (n>=7)
//  		  cantpolar=1; // 7,9,11
//  		else
//  		  cantpolar=2; // 1,3,5
//  	      }
//  	      if (n==3 || n==9) cantpolar*=-1; // 3 and 9

	      cantpolar=n%3; // 3 different types of antennas in each module, dipoles and x and y polarized slots

	    }
	    if (TYPEOFANTENNAS==1) {
	      cantpolar=n; // one seavey antenna and one dipole per station.  3 polarizations total (2 from the seavey and one from the dipole)
	    }
	    if (TYPEOFANTENNAS==2) {
	      cantpolar=n%2; // just seavey antennas, so two polarizations
	    }
	    if (TYPEOFANTENNAS==3) {
	      cantpolar=n%8; // seavey antennas in 8 different orientations
	    }


	    rxtype[3].push_back(cantpolar); //cantpolar

	    //cout << rxpos[0][NANTENNA] << " ";
	    //cout << rxpos[1][NANTENNA] << " ";
	    //cout << endl;
	    if (writeout==1) {
	      fout2 << rxpos[0][NANTENNA] << " ";
	      fout2 << rxpos[1][NANTENNA] << " ";
	      fout2 << rxpos[2][NANTENNA] << " ";
	      fout2 << rxtype[0][NANTENNA] << " ";
	      fout2 << rxtype[1][NANTENNA] << " ";
	      fout2 << rxtype[2][NANTENNA] << " ";
	      fout2 << rxtype[3][NANTENNA] << " ";
	      fout2 << endl;
	    }
	    NANTENNA++;
	  }
	  NNODESTOTAL++; // counting nodes
	}
	NSTRING++;
	 
      }
    }
    if (writeout==1) fout2.close();
  }
  else if (SHAPE==1) {// attempt at making a bicycle wheel-type detector shape.  Haven't gotten this to work yet.
    double phi;
    double r;
   
    for (int i=0;i<NFOLD;i++) {
      phi=(double)i/(double)NFOLD*TWOPI;
      for (int j=0;j<NRAD;j++) {
	 
	r=MAXARRAY[0]*((double)(j+1)/(double)NRAD);
	// each of these (in the following loop) are on the same string
	for (int k=0;k<NNODESEACHSTRING;k++) {
	  for (int n=0;n<NEACHNODE;n++) {
	    rxpos[0].push_back(r*cos(phi));
	    rxpos[1].push_back(r*sin(phi));
	    rxpos[2].push_back(-MAXARRAY[2] 
			       // begin at bottom of array
			       +(2.*MAXARRAY[2]/NNODESEACHSTRING)*((double)k+0.5) 
			       // to z position of module
			       -(SPACINGWITHINMODULE*(double)NEACHNODE)/2. 
			       // to bottom of module
			       +SPACINGWITHINMODULE*((double)n+0.5)/(double)NEACHNODE); 
	    // to position of antenna within module

	    rxtype[0].push_back(TYPEOFANTENNAS); //cantid
	    rxtype[1].push_back(NNODESTOTAL);        //cmodule
	    rxtype[2].push_back(NSTRING);        //cstring


	    int cantpolar=-50;
	    if (TYPEOFANTENNAS==0) {
	      if (n%2==0) { 
		cantpolar=0;  //0,2,4,6,8,10
	      } else {
		if (n>6)
		  cantpolar=1; // 7,9,11
		else
		  cantpolar=2; // 1,3,5
	      }
	      if (n==3 || n==9) cantpolar*=-1; // 3 and 9
	    }
	    if (TYPEOFANTENNAS==1) {
	      cantpolar=n;
	    }
	    if (TYPEOFANTENNAS==2) {
	      cantpolar=n%2;
	    }

	    rxtype[3].push_back(cantpolar); //cantpolar

	     
	    NANTENNA++;
	  }
	  NNODESTOTAL++;
	}
	NSTRING++;
      }
    }
  } else if (SHAPE==2) { // take shape from the inputs.geom file
    ifstream fin2("inputs.geom");
    TString wordin[7];
    while (!fin2.eof()) {
      for (int i=0; i<7; ++i)
	fin2 >> wordin[i];
      if (fin2.eof()) break;
      rxpos[0].push_back((double)atof(wordin[0]));
      rxpos[1].push_back((double)atof(wordin[1]));
      rxpos[2].push_back((double)atof(wordin[2]));
      rxtype[0].push_back(atoi(wordin[3]));
      rxtype[1].push_back(atoi(wordin[4]));
      rxtype[2].push_back(atoi(wordin[5]));
      rxtype[3].push_back(atoi(wordin[6]));

      if (1==2) {
	for (int i=0; i<7; ++i)
	  cout << wordin[i] << " ";
	cout << endl;

	cout << NANTENNA << " ";
	cout << rxpos[0][NANTENNA] << " ";
	cout << rxpos[1][NANTENNA] << " ";
	cout << rxpos[2][NANTENNA] << " ";
	cout << rxtype[0][NANTENNA] << " ";
	cout << rxtype[1][NANTENNA] << " ";
	cout << rxtype[2][NANTENNA] << " ";
	cout << rxtype[3][NANTENNA] << " ";
	cout << endl;
      }
      ++NANTENNA;
    }
    fin2.close();
    
    NSTRING=1;
    NNODESTOTAL=1;
    for (int ii=0;ii<NANTENNA-1;ii++) {
      if (rxtype[2][ii]!=rxtype[2][ii+1]) ++NSTRING; // new string number, so increment
      if (rxtype[1][ii]!=rxtype[1][ii+1]) ++NNODESTOTAL; // new module number, so increment
    }
  }
  else if (SHAPE==3) {

    double spacing=1333.;
    vector<double> rhex[3];
    vector<double> rhex_layer[3];

    rhex[0].push_back(-0.5*spacing);
    rhex[1].push_back(sqrt(3.)/2.*spacing);
    rhex[2].push_back(0.);

    rhex[0].push_back(0.5*spacing);
    rhex[1].push_back(sqrt(3.)/2.*spacing);
    rhex[2].push_back(0.);

    rhex[0].push_back(1.*spacing);
    rhex[1].push_back(0.);
    rhex[2].push_back(0.);

    rhex[0].push_back(0.5*spacing);
    rhex[1].push_back(-1.*sqrt(3.)/2.*spacing);
    rhex[2].push_back(0.);

    rhex[0].push_back(-0.5*spacing);
    rhex[1].push_back(-1.*sqrt(3.)/2.*spacing);
    rhex[2].push_back(0.);

    rhex[0].push_back(-1.*spacing);
    rhex[1].push_back(0.);
    rhex[2].push_back(0.);

    double rhex_next;
    double rhex_this;
   for (int i=0;i<NHEX;i++) { // loop over hexagons - for now we'll do 3
   
     for (int idim=0;idim<3;idim++) {
       rhex_layer[idim].clear();
     }
     int count[3]={0,0,0};
     for (int j=0;j<6;j++) { 

       for (int idim=0;idim<3;idim++) { // loop over 3 dimensions

	 rhex_layer[idim].push_back(rhex[idim][j]*(double)(i+1));
	 if (i==1)
	   cout << "idim, j, rhex_layer is " << idim << " " << j << " " << rhex_layer[idim][rhex_layer[idim].size()-1] << "\n";

	 // now have to fill in the stations in between
	 // on the zeroth layer there are 0 to fill in
	 // 1st layer there is 1
	 // 2nd layer there are 2
	 rhex_next=rhex[idim][(j+1)%6]*((double)(i+1));
	 rhex_this=rhex[idim][j]*((double)(i+1));


	 for (int ifill=0;ifill<i;ifill++) {


	   rhex_layer[idim].push_back(rhex_this+(rhex_next-rhex_this)*(double)(ifill+1)/(double)(i+1));
	   cout << "idim, ifill, i, rhex_this, rhex_next, rhex_layer is " << idim << " " << ifill << " " << i << " " << rhex_this << " " << rhex_next << " " <<  rhex_this+(rhex_next-rhex_this)*(double)(ifill+1)/(double)(i+1) << " " << rhex_layer[idim][rhex_layer[idim].size()-1] << "\n";
	 } // number of stations to fill in
       } // loop over 3 dimensions	   
     } // loop over 6 corners of the hexagon

     cout << "size, nodeseachstring, neachnode are " << rhex_layer[0].size() << " " << NNODESEACHSTRING << " " << NEACHNODE << "\n";
     for (int istation=0;istation<rhex_layer[0].size();istation++) {
       for (int k=0;k<NNODESEACHSTRING;k++) { // loop over nodes along z
	 
	 for (int n=0;n<NEACHNODE;n++) { // number of antennas in each node
	   // rxpos will be converted to mine-centric coordinates later
	   // these coordinates are detector-centric
	   
	   // position the antennas in x,y
	   //rxpos[0].push_back(-MAXARRAY[0]+(2.*MAXARRAY[0]/(NRXX-1))*((double)i+0.5));
	   //rxpos[1].push_back(-MAXARRAY[1]+(2.*MAXARRAY[1]/(NRXY-1))*((double)j+0.5));
	       //rxpos[0].push_back(-MAXARRAY[0]+(2.*MAXARRAY[0]/NRXX)*((double)i+0.5));
	       //rxpos[1].push_back(-MAXARRAY[1]+(2.*MAXARRAY[1]/NRXY)*((double)j+0.5));
	   rxpos[0].push_back(rhex_layer[0][istation]);
	   rxpos[1].push_back(rhex_layer[1][istation]);
	   if (i==1) 
	     cout << "rxpos's are " << i << " " << istation << " " << rxpos[0][istation] << " " << rxpos[1][istation] << "\n";
	   
	   // next we position the nodes in z
	   // for z>1, the position of the highest node will be at +MAXARRAY[2]
	   // for z=1, we position the node in the center of the mine.
	   if (NNODESEACHSTRING>1)
	     rxpos[2].push_back(-MAXARRAY[2] 
				// begin at bottom of array
				+(2.*MAXARRAY[2]/(NNODESEACHSTRING-1))*((double)k) 
				// to z position of module
				-(SPACINGWITHINMODULE*(double)NEACHNODE)/2. 
				// to bottom of module
				+SPACINGWITHINMODULE*((double)n+0.5));
	   else
	     rxpos[2].push_back(-(SPACINGWITHINMODULE*(double)NEACHNODE)/2. 
				// to bottom of module
				+SPACINGWITHINMODULE*((double)n+0.5));
	   
	   
	   rxtype[0].push_back(TYPEOFANTENNAS); //cantid
	   rxtype[1].push_back(NNODESTOTAL);        //cmodule
	   
	   rxtype[2].push_back(NSTRING);        //cstring
	       
	   int cantpolar=-50;
	   if (TYPEOFANTENNAS==0) {
	     //  	      if (n%2==0) { 
	     //  		cantpolar=0;  //0,2,4,6,8,10
	     //  	      } 
	     //  	      else {
	     //  		if (n>=7)
	     //  		  cantpolar=1; // 7,9,11
	     //  		else
	     //  		  cantpolar=2; // 1,3,5
	     //  	      }
	     //  	      if (n==3 || n==9) cantpolar*=-1; // 3 and 9
	     
	     cantpolar=n%3; // 3 different types of antennas in each module, dipoles and x and y polarized slots
	     
	   }
	   if (TYPEOFANTENNAS==1) {
	     cantpolar=n; // one seavey antenna and one dipole per station.  3 polarizations total (2 from the seavey and one from the dipole)
	   }
	   if (TYPEOFANTENNAS==2) {
	     cantpolar=n%2; // just seavey antennas, so two polarizations
	   }
	   if (TYPEOFANTENNAS==3) {
	     cantpolar=n%8; // seavey antennas in 8 different orientations
	   }
	   
	   
	   rxtype[3].push_back(cantpolar); //cantpolar
	   
	   //cout << rxpos[0][NANTENNA] << " ";
	   //cout << rxpos[1][NANTENNA] << " ";
	   //cout << endl;
	   if (writeout==1) {
	     fout2 << rxpos[0][NANTENNA] << " ";
	     fout2 << rxpos[1][NANTENNA] << " ";
	     fout2 << rxpos[2][NANTENNA] << " ";
	     fout2 << rxtype[0][NANTENNA] << " ";
	     fout2 << rxtype[1][NANTENNA] << " ";
	     fout2 << rxtype[2][NANTENNA] << " ";
	     fout2 << rxtype[3][NANTENNA] << " ";
	     fout2 << endl;
	   }
	   cout << "NANTENNA is " << NANTENNA << "\n";
	   NANTENNA++;
	     }
	 NNODESTOTAL++; // counting nodes
       }
       NSTRING++;
       
     }
   }
  }


  // Offset in z
  for (int ii=0;ii<NANTENNA;ii++) {
    if (1==2) {
      cout << ii << " ";
      cout << rxpos[2][ii] << " ";
      cout << -1*MAXARRAY[2] << " ";
      cout << +MAXPOS[2] << " ";
      cout << -INSTRUMENTEDREGION_DEPTH << " ";
      cout << +MINE_DEPTH;
      cout << endl;
    }
    rxpos[2][ii]+=
      (-1*MAXARRAY[2]+MAXPOS[2]-INSTRUMENTEDREGION_DEPTH+MINE_DEPTH); // put position of the antennas in the mine-centric coordinate system.
    //-MAXARRAY[2]-INSTRUMENTEDREGION_DEPTH is the depth of the center of instrumented region
    //MAXPOS[2]+MINE_DEPTH is the depth of the center of the mine
    
  }

  cout << NANTENNA << " antennas" << endl;
  cout << NSTRING << " strings" << endl;
  cout << NNODESTOTAL << " nodes" << endl;
  cout << NNODESTOTAL/NSTRING << " nodes/string" << endl;
  cout << NANTENNA/NNODESEACHSTRING << " antennas/node/string" << endl;


  Double_t maxz=-1E9;
  for (Int_t i=0; i<NANTENNA; ++i) {
    if (rxpos[2][i]>maxz) maxz=rxpos[2][i];
  }

  Double_t v1=MAXPOS[2];
  Double_t v2=MINE_DEPTH;
  Double_t v3=INSTRUMENTEDREGION_DEPTH;
  Double_t v4=MAXPOS[2] + MINE_DEPTH - INSTRUMENTEDREGION_DEPTH;


  cout << "Height of the top of the mine relative to the center of the mine is " << MAXPOS[2] << "\n";
  cout << "Depth of the mine relative to the surface is " << MINE_DEPTH << "\n";
  cout << "Depth of instrumented region relative to the surface is " << INSTRUMENTEDREGION_DEPTH << "\n";
  cout << "maximum value of rxpos[2][i] relative to the center of the mine is " << maxz << "\n";

  cout << "The uppermost z of any antenna is " << maxz << endl;
  cout << "Should be MAXPOS + MINE_DEPTH - INSTRUMENTEDREGION_DEPTH = " <<endl;
  cout << "   " << v1 << " + " << v2 << " - " << v3 << " = " << v4 <<endl;

  if (fabs(v4-maxz)>5) { // within 5 m and I will not complain
    cout << "ERROR GEOMETRY IS BAD " << endl;
    cout << "ERROR GEOMETRY IS BAD " << endl;
    cout << "ERROR GEOMETRY IS BAD " << endl;
    cout << "ERROR GEOMETRY IS BAD " << endl;
  }

}


void CloseTFile(TFile *hfile) {
  hfile->Write();
  hfile->Close();
  
}


int GetChordWeight(double pnu,double d1,double d2,double chord_kgm2,double chord,double len_int_kgm2,double &weight2,double &remaining) {

  // only bother continuing if interaction is within shell (meters) of 
  //  the neutrino leaving the surface
  //  so first find out where it interacts.  This will be in an
  //  exponential distribution.  (Is there are Ticho correction?) 
  //        
  // do this part as a weighted MC.  Because this is an inefficient
  //  process and would only keep a small fraction of events up to
  //  this point and we have already invested a lot of CPU.  So
  //  rather than quitting event when not in last shell meters, just try again.
  //  the weight takes care of the implicit bias.  For large chords,
  //   the weight is smaller than for small chords

  // for calculating the weight---see my notes

  // get exponential from 0 to shell with characteristic length len_int.   
  //  "remaining" is how far before nu would emerge other side of moon 


  double rtemp=0;
  double len_int;
  double shell_kgm2;
  remaining = 1.;

  len_int = len_int_kgm2/2200.;  

  shell_kgm2 = d2*2200.;
 
  //  cout << "len_int is " << len_int*2200/2900 << "\n";

  //rnd1=Rand3.Rndm(1);        //first or last "shell"
  //rnd1=0.2;
   
  //r1 = exp(-d1/len_int);
  //r2 = exp(-(d2+d1)/len_int);

  //rtemp=fabs(d2+d1)*1000.;     //force rtemp>d1 to start the while loop. 
  //while (rtemp >= d2+d1 || rtemp<=d1) {
  // rnd=Rand3.Rndm(1);
  //rnd=0.3;
  //rtemp=-len_int*log(r2+((r1-r2)*rnd));
   
  //}
  //remaining=chord-rtemp;
  
  
  //if (SANITY) {
  //if (remaining<chord-(d1+d2)) {
  //  cout << "error: remaining=" << remaining << "\n";
  //}
  //}
  
  // This is the weight fcn
  rtemp=chord_kgm2/len_int_kgm2; // number of interaction lengths the neutrino travels
  if (rtemp<50) {
    //  weight2=(exp(-(d2+d1)/len_int_kgm2)-exp(-d1/len_int_kgm2))/(exp(-chord_kgm2/len_int_kgm2)-1.);
    weight2 = (exp(-chord_kgm2/len_int_kgm2)-exp(-(chord_kgm2-shell_kgm2)/len_int_kgm2))/(exp(-chord_kgm2/len_int_kgm2)-1);
  }
  else { // if number of interaction lengths is more than 50, then set weight to 0 and return 0.
    weight2=0.;
    return 0;
  }
  
  // Make sure weight2 is between 0 and 1
  if (weight2<0) cout << "error weight="<<weight2<<"\n";
  if (weight2>1) cout << "error weight="<<weight2<<" "<<d1<<" "<<d2<<" "<<chord<<" "<<chord_kgm2<<" "<<shell_kgm2<<"\n";
  
    
  return 1;
}
void Summarize(double pnu,double* eventsfound,double doublesfound,double nabsorbed_noearth,double nabsorbed_earth,double nweighted,double sigma) {

  double volume=0;                       // active volume
  double ses;                          // single-event sensitivity
  double km2sr;                        // aperture km**2-sr
  double ster_temp;                           // 2pi or 4pi, depending on above
  double km3sr = 0;

  fout << "Livetime is " << LIVETIME << " seconds or " << LIVETIME/24/3600/365 << " years.\n";

  // write out summary
  fout << "Analyzed " << NNU << " events with pnu= " << pnu << "\n";
  fout << "freq,     num found  frac found";
  fout << "------------------------\n";
  //for (int ifreq=0;ifreq<NFREQ;ifreq++) {
  //if (double(NNU)!=0) 
  //fout << freqlist[ifreq]/1.E6 << " " << eventsfound[ifreq] << " " << (double)eventsfound[ifreq]/(double)NNU << "\n";
  //fout << freqlist[ifreq]/1.E6;
  //else
  //fout << "error dividing by nnnu.\n";	 
  //}

  // single event sensitivity for 150 MHz array per year
  //  everything is normalized to the active volume
  //  assume dF/dE propto E**-2
  //  assume 2pi sr or 4pi sr depending on whether below horizon is used
  //  assume 3.16e7 seconds/year
  //------------------------------------------------
  volume=1;
  for (int i=0;i<3;i++) {
    volume*=2.*MAXPOS[i];  //m**3
  }
  // ses=1/(A_eff * sr * t)
  //ses=1.0/(sigma*volume*RHOMEDIUM*(1./M_NUCL)*3.16E7);  // 1 per year
  
  
  double nweights=0; // number of events before the trigger is applied
  int maxntheta=0;
  int minntheta=0;
  if (USEBELOWHORIZ) {
    minntheta=0;
    maxntheta=100;
  }
  else {
    minntheta=50;
    maxntheta=100;
  }
  
 
 
  nweights+=nweighted;  // number of events before the trigger is applied
  
  error_plus=0;
  error_e_plus=0;
  error_mu_plus=0;
  error_tau_plus=0;
  error_minus=0;
  error_e_minus=0;
  error_mu_minus=0;
  error_tau_minus=0;

//    cout << "eventsfound_binned is ";Print(eventsfound_binned,NBINS);
//    cout << "eventsfound_binned_e is ";Print(eventsfound_binned_e,NBINS);
//    cout << "eventsfound_binned_mu is ";Print(eventsfound_binned_mu,NBINS);
//    cout << "eventsfound_binned_tau is ";Print(eventsfound_binned_tau,NBINS);
  


  for (int i=0;i<NBINS;i++) {

    if (eventsfound_binned[i]<=20) {
      error_plus+=pow(poissonerror_plus[(int)eventsfound_binned[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_e_plus+=pow(poissonerror_plus[(int)eventsfound_binned_e[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_mu_plus+=pow(poissonerror_plus[(int)eventsfound_binned_mu[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_tau_plus+=pow(poissonerror_plus[(int)eventsfound_binned_tau[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_minus+=pow(poissonerror_minus[(int)eventsfound_binned[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_e_minus+=pow(poissonerror_minus[(int)eventsfound_binned_e[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_mu_minus+=pow(poissonerror_minus[(int)eventsfound_binned_mu[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_tau_minus+=pow(poissonerror_minus[(int)eventsfound_binned_tau[i]]*pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
    }
    else {
      error_plus+=eventsfound_binned[i]*pow(pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_e_plus+=eventsfound_binned_e[i]*pow(pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_mu_plus+=eventsfound_binned_mu[i]*pow(pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);
      error_tau_plus+=eventsfound_binned_tau[i]*pow(pow(10.,((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT),2);

      error_minus=error_plus;
      error_e_minus=error_e_plus;
      error_mu_minus=error_mu_plus;
      error_tau_minus=error_tau_plus;

    }
  }
    error_plus=sqrt(error_plus);
    error_e_plus=sqrt(error_e_plus);
    error_mu_plus=sqrt(error_mu_plus);
    error_tau_plus=sqrt(error_tau_plus);
    
    error_minus=sqrt(error_minus);
    error_e_minus=sqrt(error_e_minus);
    error_mu_minus=sqrt(error_mu_minus);
    error_tau_minus=sqrt(error_tau_minus);


  
  if (USEBELOWHORIZ)
    ster_temp=4*PI;
  else
    ster_temp=2*PI;

  
  // account for efficiency
  //if (NNU != 0 && nevents!=0 && EXPONENT!=0) {
	 
  // alternate calculation
    
  fout << "eventsfound (with weights) direct, reflected, triggered by each independently, total triggered " << eventsfound[0] << " " << eventsfound[1] << " " << doublesfound << " " << eventsfound[2] << "\n";
  fout << "eventsfound (with weights) direct, reflected, triggered by a mix of the two, triggered by each independently, total triggered " << eventsfound[0] << " " << eventsfound[1] << " " << eventsfound[2]-(eventsfound[0]+eventsfound[1]-doublesfound) << " " << doublesfound << " " << eventsfound[2] << "\n";
  fout << "fractions - direct, reflected, triggered by a mix of the two, triggered by each independently, total triggered\t" << (eventsfound[0]-doublesfound)/eventsfound[2] << " " << (eventsfound[1]-doublesfound)/eventsfound[2] << " " << (eventsfound[2]-eventsfound[0]-eventsfound[1]+doublesfound)/eventsfound[2] << " " << doublesfound/eventsfound[2] << " " << eventsfound[2]/eventsfound[2] << "\n";


  km3sr=volume*TMath::Power(1.E-3,3)*RHOMEDIUM/RHOH20*(eventsfound[2])/(double)NNU*ster_temp;
  error_plus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_plus/(double)NNU;
  error_minus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_minus/(double)NNU;
  
  km3sr_e = (pow(1.e-3,3))*volume*ster_temp*(sum[0]/(double)nnu_e)*RHOMEDIUM/RHOH20;
  


  error_e_plus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_e_plus/(double)nnu_e;
  error_e_minus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_e_minus/(double)nnu_e;



     km3sr_mu = (pow(1.e-3,3))*volume*ster_temp*(sum[1]/(double)nnu_mu)*RHOMEDIUM/RHOH20;
     error_mu_plus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_mu_plus/(double)nnu_mu;
     error_mu_minus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_mu_minus/(double)nnu_mu;


     km3sr_tau = (pow(1.e-3,3))*volume*ster_temp*(sum[2]/(double)nnu_tau)*RHOMEDIUM/RHOH20;

     error_tau_plus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_tau_plus/(double)nnu_tau;

     error_tau_minus=volume*pow(1.E-3,3)*RHOMEDIUM/RHOH20*ster_temp*error_tau_minus/(double)nnu_tau;

 
     ftrigeff << "trigger efficiency is " << (eventsfound[2])/nweights << "\n";
     
     fout << "Total volume * solid angle is\t" << km3sr << " + " << error_plus << " - " << error_minus << " km^3 str\n";
     
     fout << "Total volume * solid angle for electron neutrinos is \t" << km3sr_e << " + " << error_e_plus << " - " << error_e_minus << " km^3 str\n";
     fout << "Total volume * solid angle for muon neutrinos is \t" << km3sr_mu << " + " << error_mu_plus << " - " << error_mu_minus << " km^3 str\n";
     fout << "Total volume * solid angle for tau neutrinos is \t" << km3sr_tau << " + " << error_tau_plus << " - " << error_tau_minus << " km^3 str\n";



     double len_int_salt=1.0/(sigma*RHOMEDIUM*(1./M_NUCL)*1000);
     km2sr=km3sr/len_int_salt;
     
     fout << "interaction length (salt) is " << len_int_salt << "\n";
     fout << "km2sr is " << km2sr << "\n";
     
     ses=(pnu/1.E9)/(km2sr*3.16E7);
     fout << "ses is " << ses*TMath::Power(1.E-2,2) << " GeV**2/cm**2/s/sr/GeV\n";
     //}
     //else 

     double sum_events=0;     
     if (EXPONENT==1) {
       

     double thisenergy=0;
     double thislen_int_kgm2=0;
     for (int i=0;i<8;i++) {
       thisenergy=TMath::Power(10.,12.+(double)i);
       ierr=GetSigma(thisenergy,sigma,thislen_int_kgm2);
       
       if (thislen_int_kgm2>0)
       sum_events+=EdNdEdAdt[i]/(thislen_int_kgm2/RHOH20);
     }
     sum_events*=volume*LIVETIME*RHOMEDIUM/RHOH20*(eventsfound[2])/NNU*ster_temp*TMath::Power(10.,4);

     }
     if (EXPONENT==0) {
       

     double thisenergy=0;
     double thislen_int_kgm2=0;
     for (int i=0;i<12;i++) {
       thisenergy=TMath::Power(10.,16.+((double)i)*0.5);
       ierr=GetSigma(thisenergy,sigma,thislen_int_kgm2);

       sum_events+=EdNdEdAdt[i]/(thislen_int_kgm2/RHOH20);
     }
     sum_events*=volume*LIVETIME*RHOMEDIUM/RHOH20*(eventsfound[2])/NNU*ster_temp*TMath::Power(10.,4);

     }

     fout << "Number of events detected is " << sum_events << " in " << LIVETIME/3600./24./365. << " years.\n";

     fout << "Loop over \t" << setprecision(4) << (double)count_neutrinos << " events.\n";
     fout << "Pass GetMine: \t" << setprecision(4) << (double)count_getmine << "\t" << (double)count_getmine/(double)count_neutrinos << " events.\n";
     fout << "Pass GetChord: \t" << setprecision(4) << (double)count_getchord << "\t" << (double)count_getchord/(double)count_getmine << " events.\n";
     fout << "Pass GetChordWeight: \t" << setprecision(4) << (double)count_getchordweight << "\t" << (double)count_getchordweight/(double)count_getchord << " events.\n";
     fout << "Pass Viewangle ChanceinHell cut: \t" << setprecision(4) << (double)count_viewangle_chanceinhell << "\t" << (double)count_viewangle_chanceinhell/(double)count_getchordweight << " events.\n";
     fout << "Pass Atten Factor ChanceinHell cut: \t" << setprecision(4) << (double)count_atten_factor_chanceinhell << "\t" << (double)count_atten_factor_chanceinhell/(double)count_viewangle_chanceinhell << " events.\n";
     fout << "Pass Trigger: \t" << setprecision(4) << (double)count_passing_events << "\t" << (double)count_passing_events/(double)count_atten_factor_chanceinhell << " events.\n\n";
     
     fout << "fraction of events where an antenna is hit: \t" << setprecision(4) << (double)count_anthit << "\t" << (double)count_anthit/(double)count_getchordweight << " events.\n";     
     fout << "fraction of events where a module has an antenna hit: \t" << setprecision(4);
     for (int k=0;k<5;k++) {
       fout << (double)count_modulehit[k] << "\t" << (double)count_modulehit[k]/(double)count_getchordweight << "\t";
     }
     fout << " events.\n\n";
     fout << "fraction of events where the event passes: \t" << setprecision(4);
     for (int k=0;k<5;k++) {
       fout << (double)count_passestrigger[k] << "\t" << (double)count_passestrigger[k]/(double)count_getchordweight << "\t";
     }
     fout << " events.\n\n";



}


int IsAbsorbed(double chord_kgm2,double len_int_kgm2) {
  // see if neutrino is absorbed 
  //  weighting works, but not to much purpose since nu's always
  //   interact at these energies.
  double rtemp;
  rtemp=chord_kgm2/len_int_kgm2;
  if (rtemp<=20.) {  // prevent FP warnings
    if (Rand3.Rndm()>(1.-exp(-rtemp))) {
     
      return 0;
    }
  }
 
  return 1;
}
void IsAbsorbed(double chord_kgm2,double len_int_kgm2,double &weight1) {
  // see if neutrino is absorbed 
  //  weighting works, but not to much purpose since nu's always
  //   interact at these energies.
  double rtemp;
  rtemp=chord_kgm2/len_int_kgm2;
  if (rtemp<=20)
    weight1=1.-exp(-rtemp);
  else weight1=1;

}
int Getmine(double *nnu,double *posnu, double *mine_in, double *mine_out) {

  double k1[6],kvec,k2,k3,k_in,k_out,d;

  for (int i=0;i<3;i++) {
    k1[2*i] = (MAXPOS[i]-posnu[i])/nnu[i];
    k1[2*i+1] = (-MAXPOS[i]-posnu[i])/nnu[i];
  }
		

  k2=0;
  k3=0;
  for (int i=0;i<6;i++) {
    kvec = k1[i];
    if ((posnu[0]+kvec*nnu[0])<1.0001*MAXPOS[0] && 
	(posnu[0]+kvec*nnu[0])>-1.0001*MAXPOS[0] && 
	(posnu[1]+kvec*nnu[1])<1.0001*MAXPOS[1] && 
	(posnu[1]+kvec*nnu[1])>-1.0001*MAXPOS[1] && 
	(posnu[2]+kvec*nnu[2])<1.0001*MAXPOS[2] && 
	(posnu[2]+kvec*nnu[2])>-1.0001*MAXPOS[2]) {
      if(k2==0) 
	k2 = kvec;
      else 
	k3 = kvec; 
    }
  }
 
  //Determine which is incoming, and which is outgoing.
  
  //	if ((x+k2*nx<x+k3*nx)&&nx>0) {
  if (k3>0) {
    k_in = k2;
    k_out = k3;}
  else {
    k_out = k2;
    k_in = k3;}
  
  //Now we know the entrance point at the mine wall, and the exit point
  //at the detector wall.
  
  for (int i=0;i<3;i++) {
    mine_in[i] = posnu[i]+k_in*nnu[i];
    mine_out[i] = posnu[i]+k_out*nnu[i];
  }

  d = sqrt(TMath::Power((mine_in[0]-mine_out[0]),2)+TMath::Power((mine_in[1]-mine_out[1]),2)+TMath::Power((mine_in[2]-mine_out[2]),2)); // distance between the entrance and exit point
  if (d<=1) return 0; // return 0 is the neutrino is only in the salt volume for less than 1 meter distance
  return 1;  
  
}

void GetEntryExit(double r,double *nnu,double *posnu, double *earth_in, double *earth_out){

  double k2 = 0;
  double k3 = 0;
  double a = 0;
  double b = 0;
  double c = 0;
  double k_in = 0;
  double k_out = 0;

  c = Square(posnu,3)-TMath::Power(r,2);

  if (c==0)
    cout << "WARNING!  c=0 inside Getearth.\n";

  for (int i=0;i<3;i++) {
    b+= nnu[i]*posnu[i];
  }
  b = 2.*b;
  a = 1.;
  
  k2 = (-b+sqrt(TMath::Power(b,2)-4.*a*c))/(2.*a);
  k3 = (-b-sqrt(TMath::Power(b,2)-4.*a*c))/(2.*a);
  
  
  //Determine which of these two vectors is the entrance and 
  //which is the exit and then determine the entrance and exit 
  //points on the Earth's surface.
  
  if (k3>0) {
    k_in = k2;
    k_out = k3;}
  else {
    k_out = k2;
    k_in = k3;}
  
  for (int i=0;i<3;i++) {
    earth_in[i] = posnu[i]+k_in*nnu[i];
    earth_out[i] = posnu[i]+k_out*nnu[i];
  }
  
}


int Getchord(double *earth_in,double *posnu,double& c, double& ck) {

  double chord[3];
  double step=5000;
  double sum,rtemp;
  int nsteps,istep;
  double step3[3];
  double posn[3];
  double radius2;
  

  double normalization = 1;
  for (int i=0;i<3;i++) {
    chord[i]=posnu[i]-earth_in[i];
  }
  c=sqrt(Square(chord,3)); // distance between neutrino's earth entrance point and the interaction point in meters.
  
  //if (c<=1) return 0;
  if (c>2.*R_EARTH) {
    cout << "bad chord" << " " << c << "\n";
  }

 

  if (c>step*1.1) {

    sum=0.;
    nsteps=int((c/step)+0.5);

    rtemp=1./double(nsteps);
    for (int i=0;i<3;i++) {
      step3[i]=chord[i]*rtemp;
    }
    for (int i=0;i<3;i++) {
      posn[i]=earth_in[i]-step3[i]/2.;
    }

    for (istep=1;istep<=nsteps;istep++) { 
      for (int i=0;i<3;i++) {
	posn[i]+=step3[i];
      }
      
      //check radius, but without the sqrt. using CPU time
      radius2=Square(posn,3);
      if (radius2<radii[0]) 
	sum+=densities[0];   
      else {
	if (radius2<radii[1])   
	  sum+=densities[1]; 
	else {
	  if (radius2<radii[2])   
	    sum+=densities[2]; 
	  else {
	    cout << "bad distance " << radius2<< " " << sqrt(radius2)<< " " << chord[0]<< " " <<chord[1]<< " " <<chord[2]<< "\n";
	  }}}}
    
    ck=((sum*normalization)/double(nsteps))*c;}
  
  else {
    ck=c*densities[2];}
  return 1; // shouldn't ever get here.
  
}
double Square(double *p,int i) {
  double sum = 0;
  for (int j=0;j<i;j++) {
    sum += p[j]*p[j];
  }
  return sum;

}

void ReadInputs() {

  // read inputs to the code.  
  // See comments in input file


  string number;
  string junk;

  getline(inputsfile,junk);

  fout << "\n\n";


  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );



  fout << "Current date and time are: " << asctime (timeinfo) << "\n";


  GetNextNumberAsString(inputsfile,foutput,number);

  NNU=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  EXPONENT=(double)atof(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  LIVETIME=(double)atof(number.c_str())*365.*24.*3600.;

  GetNextNumberAsString(inputsfile,foutput,number);

  HIST=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  SEED=(int)atoi(number.c_str());

  Rand3.SetSeed(SEED);

  GetNextNumberAsString(inputsfile,foutput,number);

  SHAPE=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  NFREQ=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  DETECTOR=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  SIGNAL_FLUCT=(int)atoi(number.c_str());
  
  if (SIGNAL_FLUCT!=1)
    cout << "Non-default setting:  Not adding noise to the signal!\n";

  GetNextNumberAsString(inputsfile,foutput,number);

  MAXPOS[0]=(double)atof(number.c_str());

  if (MAXPOS[0]!=2000.)
    cout << "Non-default setting:  Extent in x= +/-" << MAXPOS[0] << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);

  MAXPOS[1]=(double)atof(number.c_str());

  if (MAXPOS[1]!=2000.)
    cout << "Non-default setting:  Extent in y= +/-" << MAXPOS[1] << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  MAXPOS[2]=(double)atof(number.c_str());
  
  if (MAXPOS[2]!=5000.)
    cout << "Non-default setting:  Extent in z= +/-" << MAXPOS[2] << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  MINE_DEPTH=(double)atof(number.c_str());

  if (MINE_DEPTH!=500.)
    cout << "Non-default setting:  Mine depth is " << MINE_DEPTH << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);

  MAXARRAY[0]=(double)atof(number.c_str());

  if (MAXARRAY[0]!=1125.)
    cout << "Non-default setting:  Extent of detector in x= +/-" << MAXARRAY[0] << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);

  MAXARRAY[1]=(double)atof(number.c_str());

  if (MAXARRAY[1]!=1125.)
    cout << "Non-default setting:  Extent of detector in y= +/-" << MAXARRAY[1] << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  MAXARRAY[2]=(double)atof(number.c_str());
  
  if (MAXARRAY[2]!=1000.)
    cout << "Non-default setting:  Extent of detector in z= +/-" << MAXARRAY[2] << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  INSTRUMENTEDREGION_DEPTH=(double)atof(number.c_str());

  if (INSTRUMENTEDREGION_DEPTH!=750.)
    cout << "Non-default setting:  Instrumented region depth is " << INSTRUMENTEDREGION_DEPTH << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  NRXX=atoi(number.c_str());

  if (NRXX!=10)
    cout << "Non-default setting:  Number of strings in x is " << NRXX << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  NRXY=atoi(number.c_str());

  if (NRXY!=10)
    cout << "Non-default setting:  Number of strings in y is " << NRXY << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  NNODESEACHSTRING=atoi(number.c_str());

  if (NNODESEACHSTRING!=12)
    cout << "Non-default setting:  Number of strings in z is " << NNODESEACHSTRING << "\n";




  GetNextNumberAsString(inputsfile,foutput,number);
  
  NNODEHITSREQUIRED=atoi(number.c_str());

  if (NNODEHITSREQUIRED!=5)
    cout << "Non-default setting:  Number of modules required to be hit " << NNODEHITSREQUIRED << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  DISTCUT=atoi(number.c_str());

  if (DISTCUT!=1.E3)
    cout << "Non-default setting:  Max distance between a cluster of hit modules " << DISTCUT << "\n";

  


  GetNextNumberAsString(inputsfile,foutput,number);
  
  irxhitsrequired=atoi(number.c_str());

  if (irxhitsrequired!=5)
    cout << "Non-default setting:  Number of antennas required to be hit on one module" << irxhitsrequired << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  attenlength=(double)atof(number.c_str());

  if (attenlength!=250.)
    cout << "Non-default setting:  Attenuation length is " << attenlength << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  SALT=atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);
  
  N_NONDEFAULT=(double)atof(number.c_str());

  if (SALT!=1)
    cout << "Non-default setting:  Ice instead of salt!\n";

  if (SALT==1) {
    KELVINS=KELVINS_SALT;
    changle=changle_salt;
    RHOMEDIUM=RHOSALT;

    if (N_NONDEFAULT==0)
      NMEDIUM=NSALT;
    else
      NMEDIUM=N_NONDEFAULT;

    X0MEDIUM=X0SALT;
    ECMEDIUM=ECSALT;
    AEXMEDIUM=AEX_SALT;
    ALPHAMEDIUM=ALPHASALT;
    RM_MEDIUM=RM_SALT;
    KE_MEDIUM=KE_SALT; // constant in jaime's parameterization, in V/cm/MHz
    KL_MEDIUM=KL_SALT; //constant in jaime's parameterization
    KDELTA_MEDIUM=KDELTA_SALT; // constant in jaime's parameterization
    KR_MEDIUM=KR_SALT; // constant in jaime's parameterization
    BETAMEDIUM=BETASALT; // exponent, in jaime's parameterization
    N_RECEIVER=NMEDIUM; // since this code is for embedded detectors, assume n_receiver is the same as nmedium
  }
  else if (SALT==0) {
    KELVINS=KELVINS_ICE;
    changle=changle_ice;
    RHOMEDIUM=RHOICE;

    if (N_NONDEFAULT==0)
      NMEDIUM=NICE;
    else
      NMEDIUM=N_NONDEFAULT;

    X0MEDIUM=X0ICE;
    ECMEDIUM=ECICE;
    AEXMEDIUM=AEX_ICE;
    ALPHAMEDIUM=ALPHAICE;
    N_RECEIVER=NMEDIUM; // since this code is for embedded detectors, assume n_receiver is the same as nmedium

  }

  GetNextNumberAsString(inputsfile,foutput,number);
  
  TYPEOFANTENNAS=atoi(number.c_str());

  
  if (TYPEOFANTENNAS==0)
    NEACHNODE=3; // default, 3 polarizations per node, 1 dipole and 2 slots
  if (TYPEOFANTENNAS==1)
    NEACHNODE=3; // 3 polarizations per node, 2 from seavey and 1 dipole
  else if (TYPEOFANTENNAS==2)
    NEACHNODE=4;  // 4 polarizations per node, 2 from each of 2 seaveys
  else if (TYPEOFANTENNAS==3) 
    NEACHNODE=8; // 8  polarizations per node.
    
  GetNextNumberAsString(inputsfile,foutput,number);
  
  if (NEACHNODE!=0)
    NEACHNODE=atoi(number.c_str());



  GetNextNumberAsString(inputsfile,foutput,number);
  
  SIGMA_FACTOR=(double)atof(number.c_str());
  SIGMA_FACTOR=(double)atof(number.c_str());
  
  if (SIGMA_FACTOR!=1.)
    cout << "Warning!  Using non-default value for SIGMA_FACTOR = " << SIGMA_FACTOR << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);
  
  FOLLOWONEEVENT=atoi(number.c_str());


  GetNextNumberAsString(inputsfile,foutput,number);

  
  DEPTH_DEPENDENT_N=atoi(number.c_str());

  if ((DEPTH_DEPENDENT_N==1 && SALT==1) ||
      (DEPTH_DEPENDENT_N==0 && SALT==0))
    cout << "Warning!  Medium is " << SALT << " and " << "depth_dependent_n is " << DEPTH_DEPENDENT_N << "\n";

 GetNextNumberAsString(inputsfile,foutput,number);
  VERBOSE=atoi(number.c_str());
  cout << "VERBOSE level is = " << VERBOSE << "\n";

  GetNextNumberAsString(inputsfile,foutput,number);

  WHICHRAYS=(int)atoi(number.c_str()); // which rays to include, direct, reflected 

  GetNextNumberAsString(inputsfile,foutput,number);

  CONSTANTY=(int)atoi(number.c_str()); // whether to use contant of 0.2 for y (1) yes or (0) no

  GetNextNumberAsString(inputsfile,foutput,number);

  FORSECKEL=(int)atoi(number.c_str()); 

  GetNextNumberAsString(inputsfile,foutput,number);
  
  SHOWERTYPE=(int)atoi(number.c_str()); 

  GetNextNumberAsString(inputsfile,foutput,number);

  FREQ_LOW=(double)atof(number.c_str()); 

  GetNextNumberAsString(inputsfile,foutput,number);
  
  BW=(double)atof(number.c_str()); 


  if (FREQ_LOW>FREQ_LOW_SEAVEYS || FREQ_LOW>FREQ_LOW_DIPOLESSLOTS)
    cout << "WARNING!! FREQ_LOW is not lower than the lowest frequency that the antennas can reach!\n";

  if (FREQ_LOW+BW<FREQ_LOW_SEAVEYS+BW_SEAVEYS || FREQ_LOW+BW<FREQ_LOW_DIPOLESSLOTS+BW_DIPOLESSLOTS)
    cout << "WARNING!! BW does not span the entire bandwidth of all of the antennas!\n";




  GetNextNumberAsString(inputsfile,foutput,number);
  
  WHICHPARAMETERIZATION=atoi(number.c_str()); 

  GetNextNumberAsString(inputsfile,foutput,number);
  
  SECONDARIES=atoi(number.c_str()); 

  GetNextNumberAsString(inputsfile,foutput,number);

  ONLYFINAL=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  HIST_MAX_ENTRIES=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  TAUDECAY=(int)atoi(number.c_str());

  GetNextNumberAsString(inputsfile,foutput,number);

  SIDEWAYS=(int)atoi(number.c_str());
 

}
void GetNextNumberAsString(string& number) {

string temp;
 getline(inputsfile,temp); // get next line of input file
 
 fout << temp << "\n"; // output to output.txt

  int place=0;
  place=temp.find_first_of(" "); // find first space

  number=temp.substr(0,place); // everything up until the space is the value

}



//  double gety() {
  

//  // THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM 
//  //  Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364
//  //  (the curves are not in their later article.)
//  //  There is also a slow energy dependence.

//    double rnd;
//    double x = 0;
//    const double R1=0.36787944;  // 1/e
//    const double R2=0.63212056;  // 1-r1

//  // generate according to Ghandi fig. 6 
//  // adjust exponent until looks like the curve
//  //  and has right mean.
//  //  (Note this is not the fcn, but the inverse of the integral...)

//    rnd = Rand3.Rndm(1); // (0,1)
//    x=TMath::Power(-log(R1+rnd*R2),2.5);  
//    return x;
    
//  }
void GetDepth(double* array,double& depth) {

  depth=MAXPOS[2]-array[2]+MINE_DEPTH;

}

void EarthCurvature(double* array,double depth_temp) {

  // adjust array coordinates so that it fits to a curved earth surface at a specific depth 
  double length=R_EARTH-depth_temp; // length=distance from center of earth
  double rxposx=array[0]; // x coordinate of antenna position
  double rxposy=array[1]; // y coordinate of antenna position
  double rxdr=sqrt(rxposx*rxposx+rxposy*rxposy); // distance in horizontal plane from the center of the detector to the antenna
  if (dSquare(array)==0) cout << "Attempt of divide by zero in Earth curvature!!\n";
  double rxdtheta=asin(rxdr/sqrt(dSquare(array)));
  double rxdphi=atan2(rxposy,rxposx);

  array[0]=length*sin(rxdtheta)*cos(rxdphi);// have the array sit on a sphere of radius "length"
  array[1]=length*sin(rxdtheta)*sin(rxdphi);
  array[2]=length*cos(rxdtheta);

}
double GetNSouthPole(double depth) {

//   double n_0=1.79;
//   double a=0.465;
//   double b=0.0144;

//   return n_0-a*exp(-1.*b*depth);

  double n1=1.325;
  double n2=1.79;
  double d1=0.;
  double d2=150.;

  return n1+(n2-n1)/(d2-d1)*depth;


}
double GetNRonneIceShelf(double altitude) {
  // these are Peter's fit parameters
  double a1=0.37146;
  double b1=0.0140157;
  double n=0;

  if (-1.*altitude < FIRNDEPTH_RONNE) 
    n=NICE;
  else if (-1.*altitude >= FIRNDEPTH && altitude <=0) 
    //    N_DEPTH=NFIRN-(4.6198+13.62*(altitude_int/1000.))*
    //(altitude_int/1000.);   // Besson's equation for n(z)
    n=NFIRN_RONNE+a1*(1.0-exp(b1*altitude));   // Peter's equation for n(z)
  else if (altitude > 0)
    cout<<"Error!  N requested for position in air!\n";

  return n;
} //GetN(altitude)
int GetRayIceSide(double* n_exit2rx,
			 double* nsurf_rfexit,
			 double nexit,
			 double nenter,
			 
			 double* nrf2_iceside) {

  double costh=0;
  
  double NRATIO=nexit/nenter;

  costh=dDot(n_exit2rx,nsurf_rfexit,3)/(sqrt(Square(n_exit2rx,3)) * sqrt(Square(nsurf_rfexit,3)));
  
  if (costh<0)
    return 0;
  
  double sinth=sqrt(1 - costh*costh);
  

  double factor=NRATIO*costh-sqrt(1-(NRATIO*sinth*NRATIO*sinth));
 

  for (int i=0;i<3;i++) {
    nrf2_iceside[i] = -factor*nsurf_rfexit[i] + NRATIO*n_exit2rx[i];
  }

  for (int i=0;i<3;i++) {
    nrf2_iceside[i] = nrf2_iceside[i]/sqrt(Square(nrf2_iceside,3));
  }

  return 1;
}//GetRayIceSide
int TIR(double *n_surf,double *nrf_in,double N_IN,double N_OUT) {
   
  double rtemp=acos(dDot(nrf_in,n_surf,3));

  double test=sin(rtemp)*N_IN/N_OUT;

  if(test>=1) 
    return 1;
  else 
    return 0;
  
  return 0; 
}//TIR
