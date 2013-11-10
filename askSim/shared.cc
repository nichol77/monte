// This file contains functions that are shared by the salt and ice simulations.
inline double dMinNotZero(const double *x,int n) {
  double min=dMax(x,n);
  if (min==0)
    cout << "max is 0.\n";
  for (int k=1;k<n;k++) {
    if (x[k]<min && x[k]!=0)
      min=x[k];
  }
  return min;
} //dMinNotZero(double*, int)
inline double dMax(const double *x,int n) {
  double max=x[0];
  for (int k=1;k<n;k++) {
    if (x[k]>max)
      max=x[k];
  }
  return max;
} //dMax(double*, int)

inline double dMax(double a,double b) {
  if (a>b)
    return a;
  else if (a<b)
    return b;
  else if (a==b)
    return a;
  return 0;
} //dMax(double,double)

inline void Picky(double *y_cumulative,int NPROB,double rnd,double& y) {
  for (int i=0;i<NPROB;i++) {
    if (y_cumulative[i]<=rnd && y_cumulative[i+1]>rnd) {
      y=(double)i/(double)NPROB;
      continue; // once you found the right bin, stop looping.
    } //if
  } //for
} //Picky
inline int GetEMFrac(string nuflavor,
		     string current,
		     string taudecay,	      
		     double y,
		     TH1F *hy,

		     double& emfrac,
		     double& hadfrac,
		     int& n_interactions) {

  if (current=="cc")
    plepton=(1.-y)*pnu;
  else
    plepton=0.;
  
  if (nuflavor=="nue" && current=="cc") {
    emfrac=1.-y;
    hadfrac=y;
  }
  else if(nuflavor=="numu" && current=="cc") {
    emfrac=1.E-10;
    hadfrac=y;
  }
  else if(nuflavor=="nutau" && current=="cc") {
    // behaves like a muon
    emfrac=1.E-10;
    hadfrac=y;
  }
  else if (current=="nc") {
    emfrac=1.E-10;
    hadfrac=y;
  }

  em_secondaries_max =emfrac; // initialize search for maximum signal among primary, secondary interactions.
  had_secondaries_max=hadfrac;


  if (SECONDARIES==1 && current=="cc" && FORSECKEL!=1) {
    while (1) {
      GetSecondaries(nuflavor,plepton,em_secondaries_max,had_secondaries_max,n_interactions,hy); // find how much em and hadronic energies comes from secondary interactions.  keep picking until you get a bunch of secondary interactions that conserve energy
      if (em_secondaries_max+had_secondaries_max<=plepton) // if conserves energy, break.
	break;
      else {
	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
	em_secondaries_max=emfrac;
	had_secondaries_max=hadfrac;
      } //else
    } //while(1)

    if ((em_secondaries_max+had_secondaries_max)>(emfrac+hadfrac)*pnu) { // if maximum signal from secondaries is larger than
                                                                         // signal from primary interaction
      emfrac=em_secondaries_max/pnu; // then use that one.
      hadfrac=had_secondaries_max/pnu;
      if (emfrac <= 1.E-10)
	emfrac=1.E-10;
      if (hadfrac<= 1.E-10)
	hadfrac=1.E-10;
    } //if
  } //if (charged current, secondaries on)

  if (nuflavor=="numu" && current=="cc" && n_interactions==0)
    cout << "Look at this one.  inu is " << inu << "\n";
  
  if ((y<0 || y>1) && y != -999.) 
    cout <<  "illegal y=" << y << "\n";
          
  if (emfrac+hadfrac>1.00001) {
    cout << "error emfrac,hadfrac=" << emfrac << " " << hadfrac << " " << emfrac+hadfrac << "\n";
    cout << "nuflavor,taudecay=" << nuflavor << " " << taudecay << "\n";
  } //if
  
  return 1;

} //GetEMFrac
inline void GetSecondaries(string nuflavor,double plepton,double &em_secondaries_max,double &had_secondaries_max,int &n_interactions,TH1F *hy) {


  em_secondaries_max=0.;
  had_secondaries_max=0.;

  int i=(int)((log10(plepton)-18.)*2.);
  if (i>6)
    i=6;
  if (i<0)
    i=0;

  int n_brems,n_epair,n_pn; // number of interactions of each type.
  int index_y; // index along the horizontal axis of ped's plots
  double rnd1=1000.;
  double rnd2=1000.;  // random numbers for throwing at dart board
  double y = 0; // inelasticity
 
  string whichtype; // which type of interaction corresponds to that index
  
  if (nuflavor=="numu") {   
    n_brems=Rand3.Poisson(int_muon_brems[i]); // pick number of brem interactions
    n_epair=Rand3.Poisson(int_muon_epair[i]); // # of pair production
    n_pn=Rand3.Poisson(int_muon_pn[i]); // # photonuclear interactions   
    
    n_interactions+=(n_brems+n_epair+n_pn);					  

    for (int j=0;j<n_brems+n_epair+n_pn;j++) {
      rnd1=Rand3.Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";

      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {	
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_muon_brems[i],NPROB,rnd1,y);
      }
      else if (whichtype=="epair") {	
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_muon_epair[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="pn") {
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_muon_pn[i],NPROB,rnd1,y);
      }
     
      if (y*plepton>(em_secondaries_max+had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (whichtype=="brems" || whichtype=="epair") {  // save it
	  em_secondaries_max=y*plepton;

	}
	if (whichtype=="pn")
	  had_secondaries_max=y*plepton;	
      }
    } // loop over secondary interactions
  } // end if it was a muon neutrino
  if (nuflavor=="nutau") {
    n_brems=Rand3.Poisson(int_tauon_brems[i]);
    n_epair=Rand3.Poisson(int_tauon_epair[i]);
    n_pn=Rand3.Poisson(int_tauon_pn[i]);

    n_interactions+=(n_brems+n_epair+n_pn); // increment number of secondary interactions.

    for (int j=0;j<n_brems+n_epair+n_pn;j++) { // loop over secondary interactions. 
      
      rnd1=Rand3.Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";
  
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {  // bremstrahlung interaction
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_tauon_brems[i],NPROB,rnd1,y);
      }
      if (whichtype=="epair") { // pair production
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_tauon_epair[i],NPROB,rnd1,y);
      }
      if (whichtype=="pn") {
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_tauon_pn[i],NPROB,rnd1,y);
      }

      if (HIST==1 && !ONLYFINAL && hy->GetEntries()<HIST_MAX_ENTRIES)
	hy->Fill(y);
      if (y*plepton>(em_secondaries_max+had_secondaries_max)) { // if this is the biggest secondary signal yet,
	if (whichtype=="brems" || whichtype=="epair") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="pn")
	  had_secondaries_max=y*plepton;
      }
    }
   

    if (TAUDECAY) {
      n_interactions++; // increment number of interactions, for plotting.

      rnd1=Rand3.Rndm();
      if (rnd1<0.65011)  // pick which type of decay it is.
	whichtype="hadrdecay";
      if (rnd1>=0.65011 && rnd1<0.8219)
	whichtype="mudecay";
      if (rnd1>=0.8219)
	whichtype="edecay";
           
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;     
      
      if (whichtype=="hadrdecay") { // hadronic decay
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_tauon_hadrdecay[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="edecay") { // e decay	
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_tauon_edecay[i],NPROB,rnd1,y);
      }
      else if (whichtype=="mudecay") { // mu decay
	rnd1=Rand3.Rndm();
	Picky(y_cumulative_tauon_mudecay[i],NPROB,rnd1,y);
      }
      
     
      if (y*plepton>(em_secondaries_max+had_secondaries_max)) {  // if this is the biggest interaction yet,    
	if (whichtype=="edecay") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="hadrdecay")
	  had_secondaries_max=y*plepton;
      } //if     
    } //if (TAUDECAY)
  } //if (nutau)

} //GetSecondaries

inline void ReadSecondaries() {
  // reading in data for secondary interactions
  string istring;
  char buffer[50];
  int n;  // counter

  ifstream ifile;
  int index;
  cout<<"Reading in data on secondary interactions.\n";
  for (int j=0;j<2;j++) {
    // for each energy
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/muons/dsdy_brems_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/muons/dsdy_brems_1e%d.5.vec", i);

      istring=buffer;

      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_brems[index][NPROB] >> dsdy_muon_brems[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
  }

  for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/muons/dsdy_epair_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/muons/dsdy_epair_1e%d.5.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_epair[index][NPROB] >> dsdy_muon_epair[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
  }

  for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {      
      if (j==0)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/muons/dsdy_pn_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/muons/dsdy_pn_1e%d.5.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_pn[index][NPROB] >> dsdy_muon_pn[index][NPROB];
	NPROB++;
      }
      ifile.close();   
       }
    }
  }

   for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_brems_1e%d_tau.vec", i);
      if (j==1)
	n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_brems_1e%d.5_tau.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_tauon_brems[index][NPROB] >> dsdy_tauon_brems[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
   }

   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {       

       if (j==0)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_epair_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_epair_1e%d.5_tau.vec", i);
       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_epair[index][NPROB] >> dsdy_tauon_epair[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {  
       if (j==0)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_pn_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_pn_1e%d.5_tau.vec", i);
       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_pn[index][NPROB] >> dsdy_tauon_pn[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {

       if (j==0)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_hadrdecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_hadrdecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_hadrdecay[index][NPROB] >> dsdy_tauon_hadrdecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       
       if (!(i==21 && j==1)) {
       if (j==0)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_edecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_edecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_edecay[index][NPROB] >> dsdy_tauon_edecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {
       if (j==0)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_mudecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "/home/connolly/icemc/secondary/tauon/dsdy_mudecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_mudecay[index][NPROB] >> dsdy_tauon_mudecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
 
   // for filling vectors with y values distributed so that they follow
   // dsdy distributions.

   for (int j=0;j<7;j++) {

     int_muon_brems[j]=0;
     int_muon_epair[j]=0;
     int_muon_pn[j]=0;   
     int_tauon_brems[j]=0;
     int_tauon_epair[j]=0;
     int_tauon_pn[j]=0;
     int_tauon_hadrdecay[j]=0;
     int_tauon_edecay[j]=0;
     int_tauon_mudecay[j]=0;
     
     // integrating prob. distributions.
     for (int i=0;i<NPROB;i++) {
       int_muon_brems[j]+=dsdy_muon_brems[j][i];
       int_muon_epair[j]+=dsdy_muon_epair[j][i];
       int_muon_pn[j]+=dsdy_muon_pn[j][i];
       int_tauon_brems[j]+=dsdy_tauon_brems[j][i];
       int_tauon_epair[j]+=dsdy_tauon_epair[j][i];
       int_tauon_pn[j]+=dsdy_tauon_pn[j][i];
       int_tauon_hadrdecay[j]+=dsdy_tauon_hadrdecay[j][i];
       int_tauon_edecay[j]+=dsdy_tauon_edecay[j][i];
       int_tauon_mudecay[j]+=dsdy_tauon_mudecay[j][i];
     }

     // maximum value of prob. dist. 
     max_muon_brems=dMax(dsdy_muon_brems[j],NPROB);
     max_muon_epair=dMax(dsdy_muon_epair[j],NPROB);
     max_muon_pn=dMax(dsdy_muon_pn[j],NPROB);   
     max_tauon_brems=dMax(dsdy_tauon_brems[j],NPROB);
     max_tauon_epair=dMax(dsdy_tauon_epair[j],NPROB);
     max_tauon_pn=dMax(dsdy_tauon_pn[j],NPROB);
     max_tauon_hadrdecay=dMax(dsdy_tauon_hadrdecay[j],NPROB);
     max_tauon_edecay=dMax(dsdy_tauon_edecay[j],NPROB);
     max_tauon_mudecay=dMax(dsdy_tauon_mudecay[j],NPROB);
     
     // minimum value of prob. dist.
     min_muon_brems=dMinNotZero(dsdy_muon_brems[j],NPROB);
     min_muon_epair=dMinNotZero(dsdy_muon_epair[j],NPROB);
     min_muon_pn=dMinNotZero(dsdy_muon_pn[j],NPROB);   
     min_tauon_brems=dMinNotZero(dsdy_tauon_brems[j],NPROB);
     min_tauon_epair=dMinNotZero(dsdy_tauon_epair[j],NPROB);
     min_tauon_pn=dMinNotZero(dsdy_tauon_pn[j],NPROB);
     min_tauon_hadrdecay=dMinNotZero(dsdy_tauon_hadrdecay[j],NPROB);
     min_tauon_edecay=dMinNotZero(dsdy_tauon_edecay[j],NPROB);
     min_tauon_mudecay=dMinNotZero(dsdy_tauon_mudecay[j],NPROB);
     
     if (min_muon_brems<=0)
       cout << "Minimum probability is <=0!\n";

     
     // for each y bin in dsdy curve, fill vector y_muon_brem with as
     // many of y's as you need to get the right distribution.
     for (int i=0;i<NPROB;i++) {
       y_cumulative_muon_brems[j][i]      = dSum(dsdy_muon_brems[j],i+1);
       y_cumulative_muon_epair[j][i]      = dSum(dsdy_muon_epair[j],i+1);
       y_cumulative_muon_pn[j][i]         = dSum(dsdy_muon_pn[j],i+1);
       y_cumulative_tauon_brems[j][i]     = dSum(dsdy_tauon_brems[j],i+1);
       y_cumulative_tauon_epair[j][i]     = dSum(dsdy_tauon_epair[j],i+1);
       y_cumulative_tauon_pn[j][i]        = dSum(dsdy_tauon_pn[j],i+1);
       y_cumulative_tauon_hadrdecay[j][i] = dSum(dsdy_tauon_hadrdecay[j],i+1);
       y_cumulative_tauon_mudecay[j][i]   = dSum(dsdy_tauon_mudecay[j],i+1);
       y_cumulative_tauon_edecay[j][i]    = dSum(dsdy_tauon_edecay[j],i+1);
     } //for

     // normalize the distributions
     for (int i=0;i<NPROB;i++) {
       y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][NPROB-1];
       y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][NPROB-1];
       y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][NPROB-1];
       y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][NPROB-1];
       y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][NPROB-1];
       y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][NPROB-1];
       y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][NPROB-1];
       y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][NPROB-1];
       y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][NPROB-1];
     } //for
     
   }

   cout<<"Finished reading secondary interaction data.\n";
} //end method ReadSecondaries


double Step(double x) {
  if (x>=0.)
    return 1.;
  if (x<0.)
    return 0.;

  return 0;

}
void TaperVmMHz(double viewangle,
		double deltheta_em,
		double deltheta_had,
		double emfrac,
		double hadfrac,
		int WHICHPARAMETERIZATION,

		double& vmmhz1m) {

  //--EM
  
    vmmhz1m_em=0;
    vmmhz1m_had=0;
    double rtemp=0.5*(viewangle-changle)*(viewangle-changle)/(deltheta_em*deltheta_em);
    
    if (emfrac!=0) {
      if (rtemp<=20) {
	vmmhz1m_em=vmmhz1m*exp(-rtemp);
	//if (fabs(viewangle-CHANGLE)/deltheta_em>2.7)
	//vmmhz1m_em=vmmhz1m/1000;
	//if (fabs(viewangle-CHANGLE)>0.44)
	//vmmhz1m_em=0;
      }
      else
	vmmhz1m_em=0.;
    }
    else
      vmmhz1m_em=0;
    
    //--HAD
    rtemp=ALOG2*(viewangle-changle)*(viewangle-changle)/(deltheta_had*deltheta_had);
    
    if (hadfrac!=0) {
      if (rtemp<20) {
	vmmhz1m_had=vmmhz1m*exp(-rtemp);
	
	//if (fabs(viewangle-CHANGLE)/deltheta_had>2.7)
      //vmmhz1m_had=vmmhz1m/1000;
	//if (fabs(viewangle-CHANGLE)>0.44)
	//vmmhz1m_had=0;
	
      }
      else
	vmmhz1m_had=0.;
    }
    else 
      vmmhz1m_had=0.;
    
    logscalefactor_taper=log10((emfrac*vmmhz1m_em+hadfrac*vmmhz1m_had)/vmmhz1m);
    
    //if (logscalefactor_taper<-7)
    //cout << "this taper factor is small.  event is " << inu << "\n";
    
    vmmhz1m=sin(viewangle)*(emfrac*vmmhz1m_em+hadfrac*vmmhz1m_had);

  
} //TaperVmMHz

void ReadGains(){

  // gains from university of hawaii measurements.
  string sfrequency;
  string sgainv;
  string sgainh;
  string junk;

  getline(gainsfile,junk);

  NPOINTS_GAIN=0;
  while (!gainsfile.eof()) {
    gainsfile >> sfrequency >> sgainh >> sgainv;


    gainv_measured[NPOINTS_GAIN]=(double)atof(sgainv.c_str());
    gainh_measured[NPOINTS_GAIN]=(double)atof(sgainh.c_str());
    frequency_forgain_measured[NPOINTS_GAIN]=(double)atof(sfrequency.c_str())*1.E9;
    
    NPOINTS_GAIN++;
    getline(gainsfile,junk);
  } //while (not at end of antenna gain data file)
} //ReadGains

double GetGainV(double freq) {

  int whichbin=(int)((freq-frequency_forgain_measured[0])/(frequency_forgain_measured[1]-frequency_forgain_measured[0]));

  return gainv_measured[whichbin];  
} //GetGainV

double GetGainH(double freq) {

  int whichbin=(int)((freq-frequency_forgain_measured[0])/(frequency_forgain_measured[1]-frequency_forgain_measured[0]));

  return gainh_measured[whichbin];  
} //GetGainH

double GaintoHeight(double gain,double freq) {
  // from f=4*pi*A_eff/lambda^2
  // and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
  return 2*sqrt(gain/4/PI*CLIGHT*CLIGHT/(freq*freq)*Zr/Z0*N_RECEIVER);
} //GaintoHeight

int GetBeamWidths(double flare[4][NFREQ_MAX],double gain[2][NFREQ_MAX],double freq[NFREQ_MAX]) {

  // first component is frequency
  // second component is which plane and which polarization
  // it goes  e-plane: vp/hp, h-plane: vp/hp

  // these number were read from antenna specs

  double specs[5][4];

  specs[0][0]=57.5;
  specs[0][1]=58.5;
  specs[0][2]=66;
  specs[0][3]=57;

  specs[1][0]=33.5;
  specs[1][1]=34.5;
  specs[1][2]=36.5;
  specs[1][3]=38;
  
  specs[2][0]=50.5;
  specs[2][1]=53;
  specs[2][2]=33;
  specs[2][3]=32;

  specs[3][0]=43.5;
  specs[3][1]=43;
  specs[3][2]=39;
  specs[3][3]=41.5;

  specs[4][0]=36.5;
  specs[4][1]=46.5;
  specs[4][2]=32;
  specs[4][3]=31;

  double specs2[5][2];

  specs2[0][0]=8.5;
  specs2[0][1]=8.8;

  specs2[1][0]=11.0;
  specs2[1][1]=9.2;
  
  specs2[2][0]=9.3;
  specs2[2][1]=9.6;

  specs2[3][0]=10.1;
  specs2[3][1]=11.5;

  specs2[4][0]=8.9;
  specs2[4][1]=9.0;

  

  double freq_specs[5];

  freq_specs[0]=300.E6;
  freq_specs[1]=600.E6;
  freq_specs[2]=900.E6;
  freq_specs[3]=1200.E6;
  freq_specs[4]=1500.E6;

  double scale=0;
  
  for (int k=0;k<NFREQ;k++) {
      
    if (freq[k]<freq_specs[0]) {
      for (int j=0;j<4;j++) {
	flare[j][k]=specs[0][j]*RADDEG;
      } //for
    } //if
    else if (freq[k]>=freq_specs[3]) {
      for (int j=0;j<4;j++) {
	flare[j][k]=specs[3][j]*RADDEG;  
      } //for
    } //else if
    else {
      for (int i=0;i<4;i++) {
	if (freq[k]>=freq_specs[i] && freq[k]<freq_specs[i+1]) {
	  scale = (freq[k]-freq_specs[i])/(freq_specs[i+1]-freq_specs[i]);
	  
	  for (int j=0;j<4;j++) {
	    flare[j][k]=(specs[i][j]+scale*(specs[i+1][j]-specs[i][j]))*RADDEG;
	  } //for
	  i=4; 
	} //if 
      } //for
    } //else
  } //for (frequencies)
  

  for (int k=0;k<NFREQ;k++) {
      
    if (freq[k]<freq_specs[0]) {
      for (int j=0;j<2;j++) {
	gain[j][k]=10*log10(specs2[0][j]);
      } //for
    } //if
    else if (freq[k]>=freq_specs[3]) {
      for (int j=0;j<2;j++) {
	gain[j][k]=10*log10(specs2[3][j]);
      } //for
    } //else if
    else {
      for (int i=0;i<4;i++) {
	if (freq[k]>=freq_specs[i] && freq[k]<freq_specs[i+1]) {
	  scale = (freq[k]-freq_specs[i])/(freq_specs[i+1]-freq_specs[i]);
	  
	  for (int j=0;j<2;j++) {    
	    gain[j][k]=10*log10(specs2[i][j]+scale*(specs2[i+1][j]-specs2[i][j]));  
	  } //for
	  i=4;
	} //if
      } //for
    } //else
  } //for (frequencies)
    
  return 1;
} //GetBeamWidths

double dGetTheta(double p[3]) {
  return atan2(sqrt(p[0]*p[0]+p[1]*p[1]), p[2]);
} //dGetTheta



inline int GetSigma(double pnu,double& sigma,double &len_int_kgm2) {
  // calculate cross section
  if (pnu<1.2E15) {
    cout <<  "Need a parameterization for this energy region.\n";
    return 0;
  } //if
  else {
    // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
    //sigma=(2.501E-39)*pow(pnu/1.E9,0.3076)*SIGMA_FACTOR; // 10^18 eV - 10^21 eV (use this one for ANITA)
    sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA and ARIANNA)
  } //if
  // interaction length in kg/m^2
  len_int_kgm2=M_NUCL/sigma; // kg/m^2
  
  return 1;
} //GetSigma

//  int GetSigma(double pnu,double& sigma,double &len_int_kgm2) {
//    // calculate cross section
//    if (pnu<1.2E13) {
//      cout <<  "Need a parameterization for this energy region.\n";
//      return 0;
//    } //if
//    else {
//      // from Gandhi et al.
//      sigma=(7.84E-40)*pow(pnu/1.E9,0.363)*SIGMA_FACTOR;   // total xsect (m**2)
//    } //if
  
//    // interaction length in kg/m^2
//    len_int_kgm2=M_NUCL/sigma; // km/m^2
  
//    return 1;
//  } //GetSigma

// The interaction
double Gety() {
  // THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM 
  //  Ghandhi,Reno,Quigg,Sarcevc  hep-ph/9512364
  //  (the curves are not in their later article.)
  //  There is also a slow energy dependence.
  
  float rnd;
  float x = 0;
  const double R1=0.36787944;  // 1/e
  const double R2=0.63212056;  // 1-r1
  
  // generate according to Ghandi fig. 6 
  // adjust exponent until looks like the curve
  //  and has right mean.
  //  (Note this is not the fcn, but the inverse of the integral...)
  
  rnd = Rand3.Rndm(1); // (0,1)
  x=TMath::Power(-log(R1+rnd*R2),2.5); 
  
  sum_y+=x;
  n_foraddingy++;

  return x;   
} //Gety


void GetCurrent(string& current) {
  // choose CC or NC
  //  get from ratios in Ghandi etal paper
  double rnd=Rand3.Rndm();

  if (rnd<=0.7064) 
    current="cc";
  else
    current="nc";  

} //GetCurrent
// void GetNuParticleType(AskCons::ParticleType_t flavour) {
// //   double rnd=Rand3.Rndm();
  
// //   double ntypes=6.;

// //   for (int i=0;i<(int)ntypes;i++) {
// //     if (rnd<=(double)i/ntypes)
// //       flavour=i+11;
// //   }

// }
void GetNuFlavor(string& nuflavor) {
    // pick a neutrino type, flavor ratio 1:1:1
  double rnd=Rand3.Rndm();

  if (rnd<=(1./3.)) {  
    nnu_e++;
    nuflavor="nue";
  } //if
  else if(rnd<=(2./3.)) { 
    nnu_mu++;
    nuflavor="numu";
  } //else if
  else if(rnd<=(1.)) { 
    nnu_tau++;
    nuflavor="nutau";
  } //else if
  else
    cout << "unable to pick nu flavor\n";

} //GetNuFlavor

void GetSpread(double pnu,
	       double emfrac,
	       double hadfrac,
	       double freq,
	       double n_depth,
	       double X0DEPTH,
	       int WHICHPARAMETERIZATION,

	       double& deltheta_em_max,
	       double& deltheta_had_max) {

  //  scale by how far off Cherenkov angle this viewing antenna is
  //  c.f. A-MZ  astro-ph/9706064 and astro-ph/0003315
  //  and for non-LPM (non-EM) showers from Phys.Lett.B434,396 (1998)
  //  The lengths are different hence the angular thickness of 
  //  the shower is different.  Get the angular thickness for
  //  both the EM and hadroic parts.

  double elpm=GetLPM(X0DEPTH);

  freq=freq/1.E6;  // frequency in MHz
  double showerlength=3.1;  //shower length in meters-gets a modification
                            //for em showers due to lpm effect.
  // this shower length is chosen somewhat arbitrarily, but is 
  // approximately the length of a shower in ice.
  // Then, the coefficient out front of the equations for
  // deltheta_em_max and deltheta_had_max are set so that
  // for ice, we get the equations in astro-ph/9706064
  // with the coefficient in front being 2.7 degrees.
  // I wanted to make the dependence on the shower length
  // and index of refraction explicit, so I pulled those variables
  // out of the equations for deltheta_em_max and deltheta_had_max.

  double eshower;  // shower energy
  double nu0; // reference frequency

  eshower=emfrac*pnu; // first, consider the electromagnetic shower.

  // lengthen the shower to account for the lpm effect.
  if (pnu<eshower || !LPM) 
    showerlength/=pow((eshower/1.e15),-0.03); 
  else 
    showerlength/=pow(elpm/(0.14*(eshower)+elpm),0.3);



 
  if (hadfrac>0.00001) { // if there is a hadronic component
    
    eshower=hadfrac*pnu;  // just the energy of the hadronic component of the shower

    if (WHICHPARAMETERIZATION==0) {

      nu0=500.E6/1.E6*X0ICE/X0DEPTH; // for rego (astro-ph/9706064)
      // decoherence frequency scales with radiation length
      // when X0DEPTH=X0ICE, nu0=500 MHz as in astro-ph/9706064
      
      // these equations are in astro-ph/9706064, but we have pulled
      // out the dependence on index of refraction and shower length. 
      // note that 12.32/sqrt(pow(n_depth,2)-1)*RADDEG/showerlength=2.7 degrees.
      deltheta_em_max=12.32/sqrt(pow(n_depth,2)-1)*(nu0/freq)*RADDEG/showerlength;

    // these equations are in astro-ph/9706064, but we have pulled
    // out the dependence on index of refraction and shower length.
    double epsilon=log10(eshower/1.E12);
    if (eshower>=1E12 && eshower<100.E12) 
      deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(2.07-0.33*epsilon+(7.5e-2)*epsilon*epsilon);
    else if (eshower<100.E15) 
      deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(1.744-(1.21e-2)*epsilon);
    else if (eshower<10.E18)   
      deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(4.23-0.785*epsilon+(5.5e-2)*epsilon*epsilon);
    else {
      //  beyond param, just use value at 10 EeV since slow variation
      //  and parameterization might diverge
      //  so scale from 10 EeV at 7.5% per decade (30/4=7.5)
      deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(4.23-0.785*7.+5.5e-2*49.);  // the last part in parenthesis if the previous equation evaluated at epsilon=7.
      deltheta_had_max=deltheta_had_max*(1.+(epsilon-7.)*0.075);
    } //else : beyond paramatrization
    }
    else if (WHICHPARAMETERIZATION==1) {

      nu0=500.E6/1.E6; // for rego (astro-ph/9706064)

      deltheta_em_max=12.32*(nu0/freq)*RADDEG/showerlength;

      deltheta_em_max*=1./KDELTA_MEDIUM*KDELTA_ICE
	/X0MEDIUM*X0ICE
	/sqrt(NMEDIUM*NMEDIUM-1)*sqrt(NICE*NICE-1);



      deltheta_had_max=CLIGHT*100.// speed of light in cm/s
	/(freq*1.E6)
      *1/KDELTA_MEDIUM
	/(X0MEDIUM*100.) // radiation length in cm
	/sqrt(NMEDIUM*NMEDIUM-1.);
   

    }
  } //if (hadronic component)

  else
    deltheta_had_max=1.E-10;

} //GetSpread


//  inline void GetSpread(double pnu,
//  		      double emfrac,
//  		      double hadfrac,
//  		      double freq,
//  		      double n_depth,
//  		      double X0DEPTH,
//  		      int WHICHPARAMETERIZATION,
		      
//  		      double& deltheta_em_max,
//  		      double& deltheta_had_max) {

//    //  scale by how far off Cherenkov angle this viewing antenna is
//    //  c.f. A-MZ  astro-ph/9706064 and astro-ph/0003315
//    //  and for non-LPM (non-EM) showers from Phys.Lett.B434,396 (1998)
//    //  The lengths are different hence the angular thickness of 
//    //  the shower is different.  Get the angular thickness for
//    //  both the EM and hadroic parts.

//    double elpm=GetLPM(X0DEPTH);

//    freq=freq/1.E6;  // frequency in MHz
//    double showerlength=3.1;  //shower length in meters-gets a modification
//                              //for em showers due to lpm effect.
//    // this shower length is chosen somewhat arbitrarily, but is 
//    // approximately the length of a shower in ice.
//    // Then, the coefficient out front of the equations for
//    // deltheta_em_max and deltheta_had_max are set so that
//    // for ice, we get the equations in astro-ph/9706064
//    // with the coefficient in front being 2.7 degrees.
//    // I wanted to make the dependence on the shower length
//    // and index of refraction explicit, so I pulled those variables
//    // out of the equations for deltheta_em_max and deltheta_had_max.

//    double eshower;  // shower energy

//    double nu0=500.E6/1.E6*X0ICE/X0DEPTH; // for rego (astro-ph/9706064)
//    // decoherence frequency scales with radiation length
//    // when X0DEPTH=X0ICE, nu0=500 MHz as in astro-ph/9706064

//    eshower=emfrac*pnu; // first, consider the electromagnetic shower.

//    // lengthen the shower to account for the lpm effect.
//    if (eshower<elpm || !LPM) 
//      showerlength/=pow((eshower/elpm),-0.03); 
//    else 
//      showerlength/=pow(elpm/(0.14*(eshower)+elpm),0.3);

//    // these equations are in astro-ph/9706064, but we have pulled
//    // out the dependence on index of refraction and shower length. 
//    // note that 12.32/sqrt(pow(n_depth,2)-1)*RADDEG/showerlength=2.7 degrees.

//    //cout << "frequency, ratio is " << freq << " " << X0ICE/X0DEPTH/sqrt(pow(n_depth,2)-1) << "\n";

//    deltheta_em_max=12.32/sqrt(pow(n_depth,2)-1)*(nu0/freq)*RADDEG/showerlength;
  
//    //cout << "nu0, eshower, elpm, showerlength are " << nu0 << " " << eshower << " " << elpm << " " << showerlength << "\n";

//    if (hadfrac>0.00001) { // if there is a hadronic component
    
//      eshower=hadfrac*pnu;  // just the energy of the hadronic component of the shower

//      // these equations are in astro-ph/9706064, but we have pulled
//      // out the dependence on index of refraction and shower length.
//      double epsilon=log10(eshower/1.E12);
//      if (eshower>=1E12 && eshower<100.E12) 
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(2.07-0.33*epsilon+(7.5e-2)*epsilon*epsilon);
//      else if (eshower<100.E15) 
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(1.744-(1.21e-2)*epsilon);
//      else if (eshower<10.E18)   
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(4.23-0.785*epsilon+(5.5e-2)*epsilon*epsilon);
//      else {
//        //  beyond param, just use value at 10 EeV since slow variation
//        //  and parameterization might diverge
//        //  so scale from 10 EeV at 7.5% per decade (30/4=7.5)
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(4.23-0.785*7.+5.5e-2*49.);  // the last part in parenthesis if the previous equation evaluated at epsilon=7.
//        deltheta_had_max=deltheta_had_max*(1.+(epsilon-7.)*0.075);
//      } //else : beyond paramatrization
//    } //if (hadronic component)

//    else
//      deltheta_had_max=1.E-10;

//  } //GetSpread


//  void GetSpread(double pnu,
//  	       double emfrac,
//  	       double hadfrac,
//  	       double freq,
//  	       double n_depth,
//  	       double X0DEPTH,
//  	       int WHICHPARAMETERIZATION,

//  	       double& deltheta_em_max,
//  	       double& deltheta_had_max) {

//    //  scale by how far off Cherenkov angle this viewing antenna is
//    //  c.f. A-MZ  astro-ph/9706064 and astro-ph/0003315
//    //  and for non-LPM (non-EM) showers from Phys.Lett.B434,396 (1998)
//    //  The lengths are different hence the angular thickness of 
//    //  the shower is different.  Get the angular thickness for
//    //  both the EM and hadroic parts.

//    double elpm=GetLPM(X0DEPTH);

//    freq=freq/1.E6;  // frequency in MHz
//    double showerlength=3.1;  //shower length in meters-gets a modification
//                              //for em showers due to lpm effect.
//    // this shower length is chosen somewhat arbitrarily, but is 
//    // approximately the length of a shower in ice.
//    // Then, the coefficient out front of the equations for
//    // deltheta_em_max and deltheta_had_max are set so that
//    // for ice, we get the equations in astro-ph/9706064
//    // with the coefficient in front being 2.7 degrees.
//    // I wanted to make the dependence on the shower length
//    // and index of refraction explicit, so I pulled those variables
//    // out of the equations for deltheta_em_max and deltheta_had_max.

//    double eshower;  // shower energy

//    double nu0=500.E6/1.E6*X0ICE/X0DEPTH; // for rego (astro-ph/9706064)
//    // decoherence frequency scales with radiation length
//    // when X0DEPTH=X0ICE, nu0=500 MHz as in astro-ph/9706064


//    eshower=emfrac*pnu; // first, consider the electromagnetic shower.

//    // lengthen the shower to account for the lpm effect.
//    if (pnu<eshower || !LPM) 
//      showerlength/=pow((eshower/1.e15),-0.03); 
//    else 
//      showerlength/=pow(elpm/(0.14*(eshower)+elpm),0.3);

//    // these equations are in astro-ph/9706064, but we have pulled
//    // out the dependence on index of refraction and shower length. 
//    // note that 12.32/sqrt(pow(n_depth,2)-1)*RADDEG/showerlength=2.7 degrees.

//    deltheta_em_max=12.32/sqrt(pow(n_depth,2)-1)*(nu0/freq)*RADDEG/showerlength;



 
//    if (hadfrac>0.00001) { // if there is a hadronic component
    
//      eshower=hadfrac*pnu;  // just the energy of the hadronic component of the shower

//      if (WHICHPARAMETERIZATION==0) {

//      // these equations are in astro-ph/9706064, but we have pulled
//      // out the dependence on index of refraction and shower length.
//      double epsilon=log10(eshower/1.E12);
//      if (eshower>=1E12 && eshower<100.E12) 
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(2.07-0.33*epsilon+(7.5e-2)*epsilon*epsilon);
//      else if (eshower<100.E15) 
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(1.744-(1.21e-2)*epsilon);
//      else if (eshower<10.E18)   
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(4.23-0.785*epsilon+(5.5e-2)*epsilon*epsilon);
//      else {
//        //  beyond param, just use value at 10 EeV since slow variation
//        //  and parameterization might diverge
//        //  so scale from 10 EeV at 7.5% per decade (30/4=7.5)
//        deltheta_had_max=1.473/sqrt(pow(n_depth,2)-1)*nu0/freq*RADDEG*(4.23-0.785*7.+5.5e-2*49.);  // the last part in parenthesis if the previous equation evaluated at epsilon=7.
//        deltheta_had_max=deltheta_had_max*(1.+(epsilon-7.)*0.075);
//      } //else : beyond paramatrization
//      }
//      else if (WHICHPARAMETERIZATION==1) {

//        deltheta_had_max=CLIGHT*100.// speed of light in cm/s
//  	/(freq*1.E6)
//        *1/KDELTA_MEDIUM
//  	/(X0MEDIUM*100.) // radiation length in cm
//  	/sqrt(NMEDIUM*NMEDIUM-1.);
   

//      }
//    } //if (hadronic component)

//    else
//      deltheta_had_max=1.E-10;

//  } //GetSpread

double GetLPM(double X0DEPTH) {

  // LPM
  // elpm =7.7 TeV/cm * rho * X0 in PDG, but our x0 is in meters
  // so use elpm =  7.7TeV/cm*X0 
  // X0 is radiation lengths in cm

  //double elpm=7.7E12*(X0ICE*100.);

  //  double elpm=2.E15*(X0DEPTH/X0ICE);  // this is what Jaime uses.  see caption under figure 4 of 0003315.
  double elpm=61.5E14*X0DEPTH;  // this is what Jaime uses.  see caption under figure 4 of 0003315.

  return elpm;




} //GetLPM

double dGetPhi(double p[3]) {      
  // returns phi between 0 and 2pi.
  double pt=0;
  double phi=0;
  pt=sqrt(p[0]*p[0]+p[1]*p[1]);

  if (pt==0)
    return 0.;
  else if (pt!=0) {
    if (p[1]/pt>1 || p[1]/pt<-1)
      {
	cout << "Error in GetPhi. \n";
	return 0;
      } //if
    phi=asin(p[1]/pt);
  } //else if
  if (p[1]<0. && p[0]>0) phi += TWOPI;
  //else if (phi>0 && p[0]<0.) phi = PI - phi;
  //else if (phi<0 && p[0]<0.) phi = -(PI+phi)+TWOPI;
  else if (p[0]<0.) phi = PI - phi;

  return phi;
} //dGetPhi


double GetVmMHz1m(double pnu,double freq,double X0DEPTH,double ECDEPTH,double NDEPTH,double AEXDEPTH,int WHICHPARAMETERIZATION) {

  double vmmhz1m; // V/m/MHz at 1m

  if (WHICHPARAMETERIZATION==0) {
    // parametrization from Jaime Alvarez Munhiz  
    //  here using astro-ph/0003315 
    double nu0=1150.E6/1.E6;
    double nu0_modified=nu0
      *(X0ICE/ECICE)/(X0DEPTH/ECDEPTH)
      *(1/sqrt(NDEPTH*NDEPTH-1.))/(1/sqrt(NICE*NICE-1.));
    
    freq=freq/1.E6;  // frequency in MHz
    
    double factor=
      //1/sin(changle) // should be cerenkov angle for ice
      1/sqrt(1-1/(NICE*NICE)) // sin(changle) for ice
      *1/nu0 //
      *X0DEPTH/X0ICE  // track length *** use dE/dX rho instead
      *ECICE/ECDEPTH
      *AEXDEPTH/AEX_ICE;  // to account for critical energy
    // to account for cerenkov threshold // maybe should be "a" instead

    vmmhz1m=factor*(2.53E-7)*(pnu/1.E12)*freq*(1./(1.+pow(freq/nu0_modified,ALPHAMEDIUM)))*JAIME_FACTOR;
  }
  else if (WHICHPARAMETERIZATION==1) {

    double nu_r=(RHOMEDIUM/1000.) // density in g/cm^3
      /KR_MEDIUM/RM_MEDIUM*
      CLIGHT*100./NMEDIUM/sin(acos(1/NMEDIUM));

    vmmhz1m=KE_MEDIUM/ECMEDIUM* // KE in V/cm/MHz^2, Ec in MeV
      (X0MEDIUM*100.) // radiation length in cm
      *freq/1.E6 // frequency in MHz
      *sqrt(NMEDIUM*NMEDIUM-1)/NMEDIUM
      *pnu/1.E6; // energy in MeV

    vmmhz1m*=1./(1.+pow(freq/nu_r,ALPHAMEDIUM));
    

  }


//      // this is the old version
//      double factor=
//        X0DEPTH/X0ICE  // track length
//        *(1-1/(NDEPTH*NDEPTH))/(1-1/(NICE*NICE)) // cerenkov index of refraction factor
//        *NDEPTH/NICE // to account for cerenkov threshold
//        *ECICE/ECDEPTH;  // to account for critical energy
    
//      double vmmhz1m=factor*(2.53E-7)*(pnu/1.E12)*(freq/nu0)*(1./(1.+pow(freq/nu0_modified,1.44)))*JAIME_FACTOR;


  return vmmhz1m;
} //GetVmMHz1m
// Incident neutrinos
double GetE1() {

  double energy[12];

  for (int i=0;i<12;i++) {
    energy[i]=17.5+((double)i)/2.;
  } //for
  

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=1.;
  } 

  
  double maxflux=GetMax(EdNdEdAdt,12);

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
  } //for

  // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(GetMax(energy,12)-GetMin(energy,12));
    energybin=(int)(thisenergy/0.5);
    max=EdNdEdAdt[energybin];
    thisflux=Rand3.Rndm();
  } //while
  
  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
  } //for

  return pow(10.,thisenergy+GetMin(energy,12));
	
} //GetGZK
double GetE2() {

  double energy[12];
  double E2dNdEdAdt[12]; //log(brightness)

  for (int i=0;i<12;i++) {
    energy[i]=17.5+((double)i)/2.;
  } //for
  
  for (int i=0;i<12;i++) {
    E2dNdEdAdt[i]=1.; 
  }

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
  } 

  
  double maxflux=GetMax(EdNdEdAdt,12);

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
  } //for

  // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(GetMax(energy,12)-GetMin(energy,12));
    energybin=(int)(thisenergy/0.5);
    max=EdNdEdAdt[energybin];
    thisflux=Rand3.Rndm();
  } //while
  
  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
  } //for

  return pow(10.,thisenergy+GetMin(energy,12));
	
} //GetGZK
double GetE3() {

  double energy[12];
  double E2dNdEdAdt[12]; //log(brightness)

  for (int i=0;i<12;i++) {
    energy[i]=17.5+((double)i)/2.;
  } //for
  
  for (int i=0;i<12;i++) {
    E2dNdEdAdt[i]=-(energy[i]-9); 
  }

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
  } 

  
  double maxflux=GetMax(EdNdEdAdt,12);

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
  } //for

  // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(GetMax(energy,12)-GetMin(energy,12));
    energybin=(int)(thisenergy/0.5);
    max=EdNdEdAdt[energybin];
    thisflux=Rand3.Rndm();
  } //while
  
  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
  } //for

  return pow(10.,thisenergy+GetMin(energy,12));
	
} //GetGZK
double GetE4() {

  double energy[12];
  double E2dNdEdAdt[12]; //log(brightness)

  for (int i=0;i<12;i++) {
    energy[i]=17.5+((double)i)/2.;
  } //for
  
  for (int i=0;i<12;i++) {
    E2dNdEdAdt[i]=-2.*(energy[i]-9); 
  }

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
  } 

  
  double maxflux=GetMax(EdNdEdAdt,12);

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
  } //for

  // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(GetMax(energy,12)-GetMin(energy,12));
    energybin=(int)(thisenergy/0.5);
    max=EdNdEdAdt[energybin];
    thisflux=Rand3.Rndm();
  } //while
  
  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
  } //for

  return pow(10.,thisenergy+GetMin(energy,12));
	
} //GetGZK

// // Incident neutrinos
// double GetGZK() {

//   double energy[12];

//   double Emuons[12]; // E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
//   double Eelectrons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
  
//   for (int i=0;i<12;i++) {
//     energy[i]=16.+((double)i)/2.;
//     Emuons[i]=-30.;
//     Eelectrons[i]=-30.;
//   } //for
  
// //    E2dNdEdAdt[0]=-9.6; // 16.
// //    E2dNdEdAdt[1]=-8.9; // 16.5
// //    E2dNdEdAdt[2]=-8.1; // 17.
// //    E2dNdEdAdt[3]=-7.5; // 17.5
// //    E2dNdEdAdt[4]=-7.2; // 18.
// //    E2dNdEdAdt[5]=-6.8; // 18.5
// //    E2dNdEdAdt[6]=-6.7; // 19
// //    E2dNdEdAdt[7]=-6.8; // 19.5
// //    E2dNdEdAdt[8]=-7.2; // 20.
// //    E2dNdEdAdt[9]=-7.5; // 20.5
// //    E2dNdEdAdt[10]=-8.2; // 21.0
// //    E2dNdEdAdt[11]=-9.1; // 21.5

// //    Eelectrons[0]=-17.2; // 16.
// //    Eelectrons[1]=-17.35; // 16.5
// //    Eelectrons[2]=-17.2; // 17.
// //    Eelectrons[3]=-17.1; // 17.5
// //    Eelectrons[4]=-17.2; // 18.
// //    Eelectrons[5]=-17.5; // 18.5
// //    Eelectrons[6]=-18.0; // 19
// //    Eelectrons[7]=-18.5; // 19.5
// //    Eelectrons[8]=-19.4; // 20.
// //    Eelectrons[9]=-30.; // 20.5 punt
// //    Eelectrons[10]=-30.; // 21.0 punt
// //    Eelectrons[11]=-30.; // 21.5 punt

// //    Emuons[0]=-17.8; // 16.
// //    Emuons[1]=-17.4; // 16.5
// //    Emuons[2]=-17.; // 17.
// //    Emuons[3]=-16.75; // 17.5
// //    Emuons[4]=-16.9; // 18.
// //    Emuons[5]=-17.2; // 18.5
// //    Emuons[6]=-17.7; // 19
// //    Emuons[7]=-18.3; // 19.5
// //    Emuons[8]=-19.1; // 20.
// //    Emuons[9]=-30.; // 20.5 punt
// //    Emuons[10]=-30.; // 21.0 punt
// //    Emuons[11]=-30.; // 21.5 punt

// //    Emuons[0]=-17.1;  //16.
// //    Emuons[1]=-16.6;  //16.5
// //    Emuons[2]=-16.3;  //17.
// //    Emuons[3]=-16.2; // 17.5
// //    Emuons[4]=-16.4; // 18.
// //    Emuons[5]=-16.7; // 18.5
// //    Emuons[6]=-17.3; // 19
// //    Emuons[7]=-17.95; // 19.5
// //    Emuons[8]=-18.85; // 20.
// //    Emuons[9]=-19.9; // 20.5 punt
// //    Emuons[10]=-30.; // 21.0 punt
// //    Emuons[11]=-30.; // 21.5 punt

//   Emuons[0]=-16.85;  //16.
//   Emuons[1]=-16.4;  //16.5
//   Emuons[2]=-16.05;  //17.
//   Emuons[3]=-16.; // 17.5
//   Emuons[4]=-16.15; // 18.
//   Emuons[5]=-16.5; // 18.5
//   Emuons[6]=-17.1; // 19
//   Emuons[7]=-17.7; // 19.5
//   Emuons[8]=-18.65; // 20.
//   Emuons[9]=-19.75; // 20.5 punt
//   Emuons[10]=-30.; // 21.0 punt
//   Emuons[11]=-30.; // 21.5 punt
  





//   for (int i=0;i<12;i++) {
//     EdNdEdAdt[i]=pow(10.,Eelectrons[i])+pow(10.,Emuons[i]);
    
//   }


  


//   //    for (int i=0;i<12;i++) {
// //      EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
// //    } //for

  
//   double maxflux=GetMax(EdNdEdAdt,12);

//   for (int i=0;i<12;i++) {
//     EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
//   } //for

//   // now throw at a dartboard.
  
//   double thisenergy=16.;
//   double thisflux=2.;
//   double max=1.;
//   int energybin=0;
//   while(thisflux>max) {
//     // pick an energy  
//     thisenergy=Rand3.Rndm()*(GetMax(energy,12)-GetMin(energy,12));
//     energybin=(int)(thisenergy/0.5);
//     max=EdNdEdAdt[energybin];
//     thisflux=Rand3.Rndm();
//   } //while
  
//   for (int i=0;i<12;i++) {
//     EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
//   } //for

//   return pow(10.,thisenergy+GetMin(energy,12));
	
// } //GetGZK
// Incident neutrinos
double GetGZK() {

  //double E2dNdEdAdt[12]; //log(brightness)

  double energy[12];
  double Emuons[12]; // E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
  double Eelectrons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
  
  for (int i=0;i<12;i++) {
    energy[i]=16.+((double)i)/2.;
    Emuons[i]=-30.;
    Eelectrons[i]=-30.;
  } //for
  // what I was using previously, from upper curve of ANITA proposal
//    E2dNdEdAdt[0]=-9.6; // 16.
//    E2dNdEdAdt[1]=-8.9; // 16.5
//    E2dNdEdAdt[2]=-8.1; // 17.
//    E2dNdEdAdt[3]=-7.5; // 17.5
//    E2dNdEdAdt[4]=-7.2; // 18.
//    E2dNdEdAdt[5]=-6.8; // 18.5
//    E2dNdEdAdt[6]=-6.7; // 19
//    E2dNdEdAdt[7]=-6.8; // 19.5
//    E2dNdEdAdt[8]=-7.2; // 20.
//    E2dNdEdAdt[9]=-7.5; // 20.5
//    E2dNdEdAdt[10]=-8.2; // 21.0
//    E2dNdEdAdt[11]=-9.1; // 21.5

  // electron component of Figure 4 of ES&S
  // astro-ph/0101216
  Eelectrons[0]=-17.2; // 16.
  Eelectrons[1]=-17.35; // 16.5
  Eelectrons[2]=-17.2; // 17.
  Eelectrons[3]=-17.1; // 17.5
  Eelectrons[4]=-17.2; // 18.
  Eelectrons[5]=-17.5; // 18.5
  Eelectrons[6]=-18.0; // 19
  Eelectrons[7]=-18.5; // 19.5
  Eelectrons[8]=-19.4; // 20.
  Eelectrons[9]=-30.; // 20.5 punt
  Eelectrons[10]=-30.; // 21.0 punt
  Eelectrons[11]=-30.; // 21.5 punt

  // muon component of Figure 4 of ES&S
  // astro-ph/0101216
//    Emuons[0]=-17.8; // 16.
//    Emuons[1]=-17.4; // 16.5
//    Emuons[2]=-17.; // 17.
//    Emuons[3]=-16.75; // 17.5
//    Emuons[4]=-16.9; // 18.
//    Emuons[5]=-17.2; // 18.5
//    Emuons[6]=-17.7; // 19
//    Emuons[7]=-18.3; // 19.5
//    Emuons[8]=-19.1; // 20.
//    Emuons[9]=-30.; // 20.5 punt
//    Emuons[10]=-30.; // 21.0 punt
//    Emuons[11]=-30.; // 21.5 punt

  // lower curve of Figure 9 of ES&S
  // astro-ph/0101216
//    Emuons[0]=-17.1;  //16.
//    Emuons[1]=-16.6;  //16.5
//    Emuons[2]=-16.3;  //17.
//    Emuons[3]=-16.2; // 17.5
//    Emuons[4]=-16.4; // 18.
//    Emuons[5]=-16.7; // 18.5
//    Emuons[6]=-17.3; // 19
//    Emuons[7]=-17.95; // 19.5
//    Emuons[8]=-18.85; // 20.
//    Emuons[9]=-19.9; // 20.5 punt
//    Emuons[10]=-30.; // 21.0 punt
//    Emuons[11]=-30.; // 21.5 punt


  // upper curve in Figure 9 of ES&S
  // astro-ph/0101216
  Emuons[0]=-16.85;  //16.
  Emuons[1]=-16.4;  //16.5
  Emuons[2]=-16.05;  //17.
  Emuons[3]=-16.; // 17.5
  Emuons[4]=-16.15; // 18.
  Emuons[5]=-16.5; // 18.5
  Emuons[6]=-17.1; // 19
  Emuons[7]=-17.7; // 19.5
  Emuons[8]=-18.65; // 20.
  Emuons[9]=-19.75; // 20.5 punt
  Emuons[10]=-30.; // 21.0 punt
  Emuons[11]=-30.; // 21.5 punt

  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=pow(10.,Eelectrons[i])+pow(10.,Emuons[i]);
//      cout << "EdNdEdAdt is " << EdNdEdAdt[i] << "\n";
//      cout << "energy is " << energy[i] << "\n";
//      cout << "log(EdNdEdAdt) is " << log10(EdNdEdAdt[i]) << "\n";
  }

  //  for (int i=0;i<12;i++) {
  //EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
//      cout << "EdNdEdAdt is " << EdNdEdAdt[i] << "\n";
//      cout << "energy is " << energy[i] << "\n";
//      cout << "log(EdNdEdAdt) is " << log10(EdNdEdAdt[i]) << "\n";
  //} //for

  
   maxflux=GetMax(EdNdEdAdt,12);

    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
    } //for
//   // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(GetMax(energy,12)-GetMin(energy,12));
    energybin=(int)(thisenergy/0.5);
    max=EdNdEdAdt[energybin];
    thisflux=Rand3.Rndm();
  } //while
  
  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
  } //for

  return pow(10.,thisenergy+GetMin(energy,12));

}
// Incident neutrinos
double GetGZKIron() {

  double energy[8];
 
  for (int i=0;i<8;i++) {
    energy[i]=12.+((double)i);
  } //for
  
  EdNdEdAdt[0]=-17.1; // 12.
  EdNdEdAdt[1]=-16.5; // 13.
  EdNdEdAdt[2]=-16.1; // 14.
  EdNdEdAdt[3]=-16.5; // 15.
  EdNdEdAdt[4]=-17.3; // 16.
  EdNdEdAdt[5]=-16.8; // 17.
  EdNdEdAdt[6]=-17.2; // 18.
  EdNdEdAdt[7]=-19.; // 19.


  for (int i=0;i<12;i++) {
    EdNdEdAdt[i]=pow(10,EdNdEdAdt[i]);
  } //for
  
  double maxflux=GetMax(EdNdEdAdt,8);

  for (int i=0;i<8;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]/maxflux;
  } //for

  // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(GetMax(energy,8)-GetMin(energy,8));
    energybin=(int)(thisenergy);
    max=EdNdEdAdt[energybin];
    thisflux=Rand3.Rndm(); 
  } //while
  
  for (int i=0;i<8;i++) {
    EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
  } //for 

  return pow(10.,thisenergy+GetMin(energy,8));
} //GetGZKIron

////////////////////////////////////////////////////////////////////////
// Simple tasks that I just don't want to repeat many times in the code
//
void Print(int *p,int i) {
  for (int j=0;j<i;j++) {
    cout << p[j] << " ";
  } //for
  cout << "\n";
} //Print (int*,int)

void Print(double *p,int i) {
  for (int j=0;j<i;j++) {
    cout << p[j] << " ";
  } //for
  cout << "\n";
} //Print (double*,int)

void Zero(int *anarray,int n) {
  for (int i=0;i<n;i++) {
    anarray[i]=0;
  } //for
} //Zero (int*,int)

void Zero(double *anarray,int n) {
  for (int i=0;i<n;i++) {
    anarray[i]=0;
  } //for
} //Zero (double*,int)

void GetNextNumberAsString(ifstream& fin,ofstream& fout,string& number) {
  string temp;  
  getline(fin,temp); // get next line of the input file 
 
  fout << temp << "\n"; // output this line to the summary file

  int place=0; 
  place=temp.find_first_of(" \t"); // find where the first space it

  number=temp.substr(0,place); // everything up until the first space is what we're looking for
} //GetNextNumberAsString

void GetNext2NumbersAsString(ifstream& fin,ofstream& fout,string& number1,string& number2, string& stherest) {

  string temp;  
  getline(fin,temp); // get next line of the input file 
 
  fout << temp << "\n"; // output this line to the summary file

  int place=0; 
  place=temp.find_first_of(" \t"); // find where the first space it

  number1=temp.substr(0,place); // everything up until the first space is what we're looking for

  temp=temp.substr(place+1,temp.size());
 

  number2=temp.substr(0,temp.find_first_of(" "));

  stherest=temp.substr(2,temp.size());
} //GetNext2NumbersAsString

// sums the first n elements of the array thisarray.
double dSum(double* thisarray,int n) {

  double sum=0;
  for (int i=0;i<n;i++) {
    sum+=thisarray[i];
  } //for
  return sum;
} //dSum

int iSum(int* thisarray,int n) {

  int sum=0;
  for (int i=0;i<n;i++) {
    sum+=thisarray[i];
  } //for
  return sum;
} //iSum

double dSquare(double *p) {
  return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
} //dSquare

double dDot(double *a,double *b, int n) {
  double dtemp=0;
  for (int k=0;k<3;k++) {
    dtemp+=a[k]*b[k];
  } //for
  return dtemp;
} //dDot

void dCross(double *a,double *b,double *c) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=-1*(a[0]*b[2]-a[2]*b[0]);
  c[2]=a[0]*b[1]-a[1]*b[0];
} //dCross

double GetMin(double *x,int n) {
  double min=1.E22;
  for (int i=0;i<n;i++) {
    if (x[i]<min)
      min=x[i];
  } //for
  return min;
} //GetMin

double GetMax(double *x,int n) {
  double max=-1.E22;
  for (int i=0;i<n;i++) {
    if (x[i]>max)
      max=x[i];
  } //for
  return max;
} //GetMax
int WhichIsMax(double *x,int n) {
  double max=-1.E22;
  int imax=0;
  for (int i=0;i<n;i++) {
    if (x[i]>max) {
      max=x[i];
      imax=i;
    }
  } //for
  return imax;
} //GetMax


/////////////////////////////////////////////////////////////
//Class Definitions
//

/////////////////////////////////////////////////////////////
//Methods for class Vector
/////////////////////////////////////////////////////////////

double Vector::operator [](int i) const {
  //code taken from ROOT's TVector3 class
  switch(i) {
  case 0:
    return x;
  case 1:
    return y;
  case 2:
    return z;
  default:
    printf("Vector::operator[](i) has been given a bad index: %d.  Returning zero.",i);
  } //end switch

  return 0.;
} //operator[]

Vector Vector::ChangeCoord(const Vector &new_z_axis) const {
  Vector new_vector = this->RotateY(new_z_axis.theta);
  new_vector = new_vector.RotateZ(new_z_axis.phi);

  return new_vector;
} //Vector::ChangeCoord

Vector Vector::RotateX(double angle) const {
  double new_x = x;
  double new_y = cos(angle)*y - sin(angle)*z;
  double new_z = sin(angle)*y + cos(angle)*z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateX

Vector Vector::RotateY(double angle) const {
  double new_x = cos(angle)*x + sin(angle)*z;
  double new_y = y;
  double new_z = -sin(angle)*x + cos(angle)*z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateY

Vector Vector::RotateZ(double angle) const {
  double new_x = cos(angle)*x - sin(angle)*y;
  double new_y = sin(angle)*x + cos(angle)*y;
  double new_z = z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateZ

Vector Vector::Rotate(const double angle,const Vector& axis) const {
  //Code blatently stolen from Root's TRotation::Rotate method
  //Example: If you rotate the vector (0,0,1) by 90 degrees around the vector (0,1,0), the result is (1,0,0).
  double length = axis.Mag();

  double s = sin(angle);
  double c = cos(angle);
  double dx = axis.x / length;
  double dy = axis.y / length;
  double dz = axis.z / length;

  double newx = (c+(1-c)*dx*dx) * x + ((1-c)*dx*dy-s*dz) * y + ((1-c)*dx*dz+s*dy) * z;
  double newy = ((1-c)*dy*dx+s*dz) * x + (c+(1-c)*dy*dy) * y + ((1-c)*dy*dz-s*dx) * z;
  double newz = ((1-c)*dz*dx-s*dy) * x + ((1-c)*dz*dy+s*dx) * y + (c+(1-c)*dz*dz) * z;

  return Vector(newx,newy,newz);
} //Vector::Rotate

void Vector::UpdateThetaPhi() {
  //This is a private method that will calculate values of theta and phi from x,y,z coordinates,
  //and store the results in the class variables.
  double transverse = sqrt(x*x + y*y);

  theta = atan2(transverse, z);

  if(transverse != 0)
    phi = atan2(y, x);
  else phi = 0;
  if(phi < 0) phi += TWOPI;
  return;
} //UpdateThetaPhi


/////////////////////////////////////////////////////////////
//Begin methods for class Position
/////////////////////////////////////////////////////////////

Position::Position() : Vector() 
{
  //This method intentionally left blank.
}

Position::Position(Vector vec) : Vector(vec.GetX(),vec.GetY(),vec.GetZ()) 
{
  //This method intentionally left blank.
}

Position::Position(double x_inp,double y_inp,double z_inp) : Vector(x_inp,y_inp,z_inp)
{
  //This method intentionally left blank.
}

double Position::GetLat() const {
  return theta*DEGRAD;
} //GetLat

double Position::GetLon() const {
  double phi_deg = phi*DEGRAD;

  if (phi_deg > 270)
    phi_deg = phi_deg-360;
  
  return (360.*(3./4.) - phi_deg);
} //GetLon

double GetDistance(double* p1, double* p2) 
{
  double distance;
  distance=sqrt(((p1[0]-p2[0])*(p1[0]-p2[0])) +
		((p1[1]-p2[1])*(p1[1]-p2[1])) +
		((p1[2]-p2[2])*(p1[2]-p2[2])));
  return distance;
  
}

void GetDirection(double* p1, double* p2, double* nhat) 
{
  // p1=rx p2=posnu

  double distance=GetDistance(p1,p2);
  for (int i=0;i<3;i++) {
    nhat[i]=(p1[i]-p2[i])/distance;
  }

}

double GetViewAngle(double* direction, double* nnu) {
  double rtemp=dDot(direction,nnu,3);
  if (fabs(rtemp)<=1) { 
    return acos(rtemp);
  } else {
    cout <<  "Danger trying to take acos of rtemp= " << rtemp << "\n";
    cout << "inu, Length of direction, nnu are " << inu << " " << dSquare(direction) << " " << dSquare(nnu) << "\n";
    return acos(1.0);
  }
}

double GetFac1(double* p1, double* p2)
{
  // two positions in meters and return 1/R factor
  double distance=GetDistance(p1,p2);
  return 1.0/distance;
}

double GetFac2(double* p1, double* p2, double atten)
{
  // two positions in meters and return attenuation factor
  double distance=GetDistance(p1,p2);
  if (distance/atten<=10) {    // if less than 10 attenuation lengths away
    return exp(-distance/atten); // get attenuation factor
  } else { // if more than 10 attenuation lengths away, then return 0
    return 0;
  }
}
	
double GetFac3(double viewangle, double changle, 
	       double pnu,       double emfrac,  double hadfrac,
	       double thisfreq,
	       double& deltheta_em,double& deltheta_had)
{
  double facs[2];
  static int init=0;
  if (init==0) {
    cout << "GetFac3 init : Global Variables are:" << endl;
    cout << "NMedium " << NMEDIUM << endl;
    cout << "X0Medium " << X0MEDIUM << endl;

    init=1;
  }
 
  GetSpread(pnu,emfrac,hadfrac,thisfreq,N_DEPTH,X0MEDIUM,WHICHPARAMETERIZATION,
	    deltheta_em,deltheta_had);

  double taperfactor=1.0;
  TaperVmMHz(viewangle,deltheta_em,deltheta_had,emfrac,hadfrac,WHICHPARAMETERIZATION,
	     taperfactor);
  // why is this here
  //  return taperfactor;

  // TaperVmMHz replaces below;

  double va=viewangle;
  double ch=changle;

  double term1=(va-ch)/deltheta_em;
  double rtemp1=ALOG2*(term1*term1);  // number of cone widths away from cherenkov angle you are
  if (rtemp1<=20) {
    facs[0]=exp(-rtemp1)*taperfactor;
  } else { // if you are more than 20 cone widths away, set this factor to 0
    facs[0]=0.0;
  }

  double term2=(va-ch)/deltheta_had;
  double rtemp2=ALOG2*(term2*term2); // number of cone widths away from cherenkov angle you are
  if (rtemp2<=20) {
    facs[1]=exp(-rtemp2)*taperfactor;
  } else {
    facs[1]=0.0; // if you are more than 20 cone widths away, set this factor to 0
  }
  return facs[0]*emfrac+facs[1]*hadfrac;

  // I include old comments here


		    
		    //		  vmmhz1m=GetVmMHz1m(pnu,freq,X0SALT,ECSALT,NMEDIUM);
		    
		    
		    //		  cout << "vmmhz1m is " << vmmhz1m << "\n";
		    
		    // but this factor freq/nu0 is the factor by which the shower is
		    //  shorter which should reduce the power by that factor but
		    //  not the e-field, so correct it back.
		    // same param for EM and HAD showers.  Only difference will
		    //  be thickness of cone            
		    // **** why the sqrt?
		    //vmmhz1m=vmmhz1m*sqrt(X0H20/X0SALT);
		    
		    // scale by how far off Cherenkov angle this viewing antenna is
		    //  c.f. A-MZ  astro-ph/9706064 and astro-ph/0003315
		    // and for non-LPM (non-EM) showers from Phys.Lett.B434,396 (1998)
		    //  The lengths are different hence the angular thickness of 
		    //  the shower is different.  Get the angular thickness for
		    //   both the EM and hadroic parts.
		    
		    // scale these widths for track length
		    //deltheta_em*=X0ICE/X0SALT;
		    //deltheta_had*=X0ICE/X0SALT;
		    
		    //cout << "deltheta_em is " << deltheta_em << " " << DEGRAD*deltheta_em << "\n";
		    //cout << "deltheta_had is " << deltheta_had << " " << DEGRAD*deltheta_had << "\n";
		    
		    //--EM
		    //		  cout << "deltheta_em, deltheta_had are " << deltheta_em << " " << deltheta_had << "\n";
		    //cout << "viewangle is " << viewangle << "\n";
		    //		    cout << "watch it.\n";

}

double GetFac4(double* direction, 
	       double hitangle_e,double hitangle_h,
	       double e_component,double h_component,
	       double* nnu, int cantid, int cantpolar, 
	       double thisfreq, int ifreq,
	       double& costh_pol, double& costh_dir, double& heff)
{


  static int init2=0;
  if (init2==0) {
    cout << "GetFac4 init : Global Variables are:" << endl;
    cout << "CLIGHT " << CLIGHT << endl;
    cout << "NFREQ " << NFREQ << endl;
    cout << "NMedium " << NMEDIUM << endl;
    cout << "X0Medium " << X0MEDIUM << endl;
    cout << "NFREQ " << NFREQ << endl;

    init2=1;
  }


  double polarfactor;
  double thislambda=CLIGHT/(NMEDIUM*thisfreq);

  double pol1[3]={0.,0.,1.};      // polarization of vert. pol. antennas
  double pol2[3]={1.,0.,0.};// polarization of horiz.pol. antennas
  double pol3[3]={0.,1.,0.};// polarization of horiz.pol. antennas

  double down[3]={0.,0.,-1.}; // normal of seavey antennas


  //double pol1[3]={0.,1.,0.};      // polarization of vert. pol. antennas
  //double pol2[3]={0.,0.,1.};// polarization of horiz.pol. antennas
  //double pol3[3]={1.,0.,0.};// polarization of horiz.pol. antennas
  //  cout << "dipoles and slots turned on their side!!!\n";

  if (dSquare(direction)>1.00001)
    cout << "Warning:  In GetFac4, direction vector is not a unit vector.  Length is " << dSquare(direction) << "\n";
  if (dSquare(nnu)>1.00001)
    cout << "Warning:  In GetFac4, nnu vector is not a unit vector.\n";


  double ctemp[3],signal_pol[3]; // ctemp is a temporary vector that is perpendicular to both the neutrino direction and the rf
  TMath::Cross(direction,nnu,ctemp); 
  Double_t tvlen=TMath::Normalize(ctemp);
 


  if (dSquare(ctemp)>1.00001)
    cout << "Warning:  In GetFac4, ctemp vector is not a unit vector.\n";


  TMath::Cross(direction,ctemp,signal_pol); // signal_pol is the polarization of the signal.  It is perp. to the line of sight in the plane of the neutrino path

  if (dSquare(signal_pol)>1.00001)
    cout << "Warning:  In GetFac4, signal_pol vector is not a unit vector.\n";


  costh_pol=dDot(signal_pol,pol1,3); // cos theta of polarization
  //cout << "watch it.\n";


  costh_dir=dDot(direction,pol1,3); // cos theta of neutrino direction

  if (cantid==0) { // cantid=0 means dipoles and slots

    heff= thislambda*sqrt(2*Zr*gain_dipole/Z0/4/PI*NMEDIUM);   // effective height of dipole, using formula from Ped's note
    if (thisfreq>(FREQ_LOW_DIPOLESSLOTS+BW_DIPOLESSLOTS)/N_RECEIVER ||
	thisfreq<FREQ_LOW_DIPOLESSLOTS/N_RECEIVER)  // if this is outside the bandwidth of the dipoles and slots then set the effective height to zero..
      heff=0.;
    GetPolarizationFactorDipolesSlots(pol1, pol2, pol3,
				      signal_pol, direction,tvlen,
				      cantpolar,polarfactor,
				      ifreq,thisfreq); // get polarization factor for dipoles or slots.  whichpol determines which antenna/polarization we are interested in here.        
  }
  if (cantid==1) { // 
    GetHitAngles(down,direction,signal_pol,
		 hitangle_e,hitangle_h,e_component,h_component); // get hit angle of the ray wrt planes of polarization
    

    //    kfreq=(int)(100.*(double)ifreq/(double)NFREQ); // frequency bin
    
 
    //cout << "ifreq, NFREQ, kfreq are " << ifreq << " " << NFREQ << " " << kfreq << "\n";



    if (cantpolar==0) { // horizontal polarization of seaveys


      //     cout << "thisfreq, gain, height is " << thisfreq << " " << GetGainH(thisfreq) << " " << GaintoHeight(GetGainH(thisfreq),thisfreq) << "\n";
      heff=GaintoHeight(GetGainH(thisfreq*N_RECEIVER),thisfreq); // horizontal polarization of quad ridged gain horn
      // uses e- and h-components instead of lcp,rcp, for Anita-lite
      
 

      polarfactor=sqrt(pow(e_component*exp(-1*ALOG2*(hitangle_e/flare[0][ifreq])*(hitangle_e/flare[0][ifreq])),2)+
		       // e_component*exp(-1*ALOG2*(hitangle_e/flare[0][ifreq])*(hitangle_e/flare[0][ifreq]))+
		       0.01*pow(h_component*exp(-1*ALOG2*(hitangle_h/flare[1][ifreq])*(hitangle_h/flare[1][ifreq])),2));


      //      cout << "inside GetFac4, cantpolar==0, polarfactor is " << polarfactor << "\n";
      
    }
    else if (cantpolar==1) { // vertical polarization of seaveys
      heff=GaintoHeight(GetGainV(thisfreq*N_RECEIVER),thisfreq); // vertical polarization of quad ridged gain horn
      polarfactor=sqrt(pow(h_component*exp(-1*ALOG2*(hitangle_h/flare[2][ifreq])*(hitangle_h/flare[2][ifreq])),2)+
		       // h_component*exp(-1*ALOG2*(hitangle_e/flare[2][ifreq])*(hitangle_e/flare[2][ifreq]))+
		       //e_component*exp(-1*ALOG2*(hitangle_h/flare[3][ifreq])*(hitangle_h/flare[3][ifreq]))*
		       0.01*pow(e_component*exp(-1*ALOG2*(hitangle_e/flare[3][ifreq])*(hitangle_e/flare[3][ifreq])),2));
      //      cout << "inside GetFac4, cantpolar==1, polarfactor is " << polarfactor << "\n";
    }
    else if (cantpolar==2) { // dipole
      heff= thislambda*sqrt(2*Zr*gain_dipole/Z0/4/PI*NMEDIUM);   // effective height of dipole, using formula from Ped's note
      polarfactor=(double)fabs(dDot(signal_pol,pol1,3));

      //cout << "inside GetFac4, cantpolar==2, polarfactor is " << polarfactor << "\n";
      if (thisfreq>(FREQ_LOW_DIPOLESSLOTS+BW_DIPOLESSLOTS)/N_RECEIVER ||
	  thisfreq<FREQ_LOW_DIPOLESSLOTS/N_RECEIVER)  // if this is outside the bandwidth of the dipoles and slots then set the effective height to zero..
	heff=0.;

    } // dipole antenna
    
  } // the seavey-dipole configuration
  if (cantid==3) { // octagon configuration
    // the antenna's axis is vertical, pointing downward.
    // GetHitAngles takes the antenna normal and calculates the hitangles and components of the signal along the e plane and h plane.  when down is vertical, the e and h planes will be along the x and y axes.  
    // for this octagon configuration, we want to have the antenna polarization always in the horizontal plane, but rotated in phi.  Instead of changing GetHitAngles, we just rotate the direction and signal_pol in phi temporarily

    double angle=(double)cantpolar*2.*PI/8.;

    double signal_pol_temp[3];
    double direction_temp[3];

    for (int i=0;i<3;i++) {
      signal_pol_temp[i]=signal_pol[i];
      direction_temp[i]=direction[i];
    }

    signal_pol_temp[0]=cos(angle)*signal_pol[0]-sin(angle)*signal_pol[1];
    signal_pol_temp[1]=sin(angle)*signal_pol[0]+cos(angle)*signal_pol[1];
    signal_pol_temp[2]=signal_pol[2];

    direction_temp[0]=cos(angle)*direction[0]-sin(angle)*direction[1];
    direction_temp[1]=sin(angle)*direction[0]+cos(angle)*direction[1];
    direction_temp[2]=direction[2];


    GetHitAngles(down,direction_temp,signal_pol_temp,
		 hitangle_e,hitangle_h,e_component,h_component); // get hit angle of the ray wrt planes of polarization

    heff=GaintoHeight(GetGainV(thisfreq*N_RECEIVER),thisfreq); // vertical polarization of quad ridged gain horn

    polarfactor=sqrt(pow(h_component*exp(-1*ALOG2*(hitangle_h/flare[2][ifreq])*(hitangle_h/flare[2][ifreq])),2)+
		     // h_component*exp(-1*ALOG2*(hitangle_e/flare[2][ifreq])*(hitangle_e/flare[2][ifreq]))+
		     //e_component*exp(-1*ALOG2*(hitangle_h/flare[3][ifreq])*(hitangle_h/flare[3][ifreq]))*
		     0.01*pow(e_component*exp(-1*ALOG2*(hitangle_e/flare[3][ifreq])*(hitangle_e/flare[3][ifreq])),2));


  }

  return polarfactor;
	    
}


void GetPolarizationFactorDipolesSlots(double *pol1, double *pol2, double *pol3,
				       double *signal_pol, double *direction,double tvlen,
				       int cantpolar,double& polarfactor,
				       int ifreq,double freq) {
  //  double opposite_pol1[3]={0.,0.,-1.};// polarization of horiz.pol. antennas
  double opposite_pol2[3]={-1.,0.,0.};// polarization of horiz.pol. antennas
  double opposite_pol3[3]={0.,-1.,0.};// polarization of horiz.pol. antennas
  double ctemp[3];
  
  
  
  polarfactor=1;
  
  if (tvlen>0) {
    if (cantpolar==0) { // if antenna number is even
      // even are dipoles

      if (SIDEWAYS==0) {
	  polarfactor*=fabs(dDot(signal_pol,pol1,3));
       
      }
      else {
	// if we are laying antennas sideways (like for a mine), turn half of them in the +x direction and half of them in the +y direction.
	if (gRandom->Rndm()>0.5)
	  polarfactor*=fabs(dDot(signal_pol,pol2,3));
	else
	  polarfactor*=fabs(dDot(signal_pol,pol3,3));
      }
    } // end if it is a dipole
    else if (abs(cantpolar)==1) { // if antenna number is 1 or -1
      // slot antennas with beam patterns pointed in +/- y direction
      //polarfactor*=dDot(signal_pol,pol2,3);

      if (SIDEWAYS==0)
	dCross(signal_pol,pol1,temp); // slots pick out electric field in a circle around orientation of the antenna.  This on is oriented in the +z direction
      else
	dCross(signal_pol,pol2,temp); // slots pick out electric field in a circle around orientation of the antenna.  This one is oriented in the +x direction
      if (SIDEWAYS==0) {
	if (cantpolar>0) 
	  // We chop off 1/2 of a dipole and multiply the power by a factor of 2.
	  // The length of temp is the component of the signal's polarization that is perp to orientation of the antenna.
	  polarfactor*=sqrt(2.)*Step(dDot(direction,pol3,3))*sqrt(dSquare(temp));// this slot's beam pattern prefers the +y direction
	else 
	  polarfactor*=sqrt(2.)*Step(dDot(direction,opposite_pol3,3))*sqrt(dSquare(temp)); // this slot's beam pattern prefers the -y direction 
      }
      else
	polarfactor*=sqrt(2.)*Step(dDot(direction,pol1,3))*sqrt(dSquare(temp));
      // this slot's beam pattern prefers the +z direction
      
      //polarfactor*=0;

      //polarfactor*=(freqlist[ifreq]/freq);
      // Second *= is ratio of lambda for dipole vs slot
      // Pretty sure (97%) it is this ratio and not reciprocal
      // PG uses (used) 150 MHz for dipole and 220 for slot
      // Also, not sure about VmMhz1m at 2 freq...
    } // end slot antennas with y orientation 
    else if (abs(cantpolar)==2) {
      // slot antennas with beam patterns pointed in +/- x direction
      if (SIDEWAYS==0)
	dCross(signal_pol,pol1,ctemp);// this slot antenna is oriented in the +z direction
      else 
	dCross(signal_pol,pol3,ctemp); // this slot antenna is oriented in the +y direction
      if (SIDEWAYS==0) {
	if (cantpolar>0) 
	  polarfactor*=sqrt(2.)*Step(dDot(direction,pol2,3))*sqrt(dSquare(temp)); // this antenna's beam pattern prefers the +x direction
	else
	  polarfactor*=sqrt(2.)*Step(dDot(direction,opposite_pol2,3))*sqrt(dSquare(temp)); // this antenna's beam pattern prefers the -x direction
      }
      else
	  polarfactor*=sqrt(2.)*Step(dDot(direction,pol1,3))*sqrt(dSquare(temp));
    
      //polarfactor*=0;
      //polarfactor*=(freqlist[ifreq]/freq);
      
    }
    else {
      cerr << "You should never see this!";
      cerr << __LINE__ << endl;
    }
  }
  else {
    cerr << "Setting polarization (hitangle) to zero!";
    cerr << endl;
    polarfactor=0;
  }
  
  // absolute value
  if (polarfactor<0) polarfactor*=-1.0;
  if (freq>(FREQ_LOW_DIPOLESSLOTS+BW_DIPOLESSLOTS)/N_RECEIVER ||
	  freq<FREQ_LOW_DIPOLESSLOTS/N_RECEIVER)  // if this is outside the bandwidth of the dipoles and slots then set the effective height to zero..

    polarfactor=0.;
}
void GetHitAngles(double *ant_normal,double *ray,double *n_pol,
	  double& hitangle_e,double& hitangle_h,double& e_component,double& h_component) {


//void GetHitAngles(double *n_exit2bn,double *n_pol,int ilayer,int ifold,double theta_bn,double phi_bn,double phi_spin,double& hitangle_e,double& hitangle_h,double& e_component,double& h_component,double* ant_normal) {

  double n_eplane[3]={0,0,1};
  double n_hplane[3]={0,-1,0};
  double n_normal[3]={1,0,0};


  double test[3];
  double test2[3];
  double test3[3];
  double temp[3];
  for (int i=0;i<3;i++) {
    test[i]=n_eplane[i];
    test2[i]=n_hplane[i];
    test3[i]=n_normal[i];
  }
  
  double theta=dGetTheta(ant_normal);
  double phi=dGetPhi(ant_normal);


  //  cout << "theta, phi are " << theta << " " << phi << "\n";

  // now rotate for antenna's orientation on the payload.

  n_eplane[0]=cos(theta-PI/2)*test[0]+sin(theta-PI/2)*test[2];
  n_eplane[1]=test[1];
  n_eplane[2]=-sin(theta-PI/2)*test[0]+cos(theta-PI/2)*test[2];

  n_hplane[0]=cos(theta-PI/2)*test2[0]+sin(theta-PI/2)*test2[2];
  n_hplane[1]=test2[1];
  n_hplane[2]=-sin(theta-PI/2)*test2[0]+cos(theta-PI/2)*test2[2];

  n_normal[0]=cos(theta-PI/2)*test3[0]+sin(theta-PI/2)*test3[2];
  n_normal[1]=test3[1];
  n_normal[2]=-sin(theta-PI/2)*test3[0]+cos(theta-PI/2)*test3[2];




  for (int i=0;i<3;i++) {
    test[i]=n_eplane[i];
    test2[i]=n_hplane[i];
    test3[i]=n_normal[i];
  }


  n_eplane[0]=cos(phi)*test[0]-sin(phi)*test[1];
  n_eplane[1]=sin(phi)*test[0]+cos(phi)*test[1];
  n_eplane[2]=test[2];


  n_hplane[0]=cos(phi)*test2[0]-sin(phi)*test2[1];
  n_hplane[1]=sin(phi)*test2[0]+cos(phi)*test2[1];
  n_hplane[2]=test2[2];

  n_normal[0]=cos(phi)*test3[0]-sin(phi)*test3[1];
  n_normal[1]=sin(phi)*test3[0]+cos(phi)*test3[1];
  n_normal[2]=test3[2];   


//    for (int i=0;i<3;i++) {
//      test[i]=n_eplane[i];
//      test2[i]=n_hplane[i];
//      test3[i]=n_normal[i];
//    }

//   //   now rotate to balloon's position on the continent


//    n_eplane[0]=cos(PI-theta_bn)*test[0]+sin(PI-theta_bn)*test[2];
//    n_eplane[1]=test[1];
//    n_eplane[2]=-sin(PI-theta_bn)*test[0]+cos(PI-theta_bn)*test[2];

//    n_hplane[0]=cos(PI-theta_bn)*test2[0]+sin(PI-theta_bn)*test2[2];
//    n_hplane[1]=test2[1];
//    n_hplane[2]=-sin(PI-theta_bn)*test2[0]+cos(PI-theta_bn)*test2[2];

//    n_normal[0]=cos(PI-theta_bn)*test3[0]+sin(PI-theta_bn)*test3[2];
//    n_normal[1]=test3[1];
//    n_normal[2]=-sin(PI-theta_bn)*test3[0]+cos(PI-theta_bn)*test3[2];


//    for (int i=0;i<3;i++) {
//      test[i]=n_eplane[i];   
//      test2[i]=n_hplane[i];
//      test3[i]=n_normal[i];

//    }

//    n_eplane[0]=cos(phi_bn)*test[0]-sin(phi_bn)*test[1];
//    n_eplane[1]=sin(phi_bn)*test[0]+cos(phi_bn)*test[1];
//    n_eplane[2]=test[2];


//    n_hplane[0]=cos(phi_bn)*test2[0]-sin(phi_bn)*test2[1];
//    n_hplane[1]=sin(phi_bn)*test2[0]+cos(phi_bn)*test2[1];
//    n_hplane[2]=test2[2];

//    n_normal[0]=cos(phi_bn)*test3[0]-sin(phi_bn)*test3[1];
//    n_normal[1]=sin(phi_bn)*test3[0]+cos(phi_bn)*test3[1];
//    n_normal[2]=test3[2];   


//    cout << "n_eplane is ";Print(n_eplane,3);
//    cout << "n_hplane is ";Print(n_hplane,3);
//    cout << "n_normal is ";Print(n_normal,3);

  for (int i=0;i<3;i++) {
    temp[i]=-1*ray[i];    
  }

  e_component=dDot(temp,n_eplane,3);
  
  // find component along the h-plane
  
  h_component=dDot(temp,n_hplane,3);
  
  // find the component normal
  
  double n_component=dDot(temp,n_normal,3);
  
  hitangle_e=atan2(h_component, n_component);

  //    if (temp[0]>0.65 && temp[0]<0.75 && temp[1]>0.65 && temp[1]<0.75) {
  //cout << "hitangle_e is " << hitangle_e << "\n";
  //}  
  //  if (inu==1341) {
  //cout << "e_component, n_component,phi,phi(temp),n_eplane are " << e_component << " " << n_component << " " << phi << " ";Print(temp,3);cout << " ";Print(n_eplane,3);cout << " ";Print(n_normal,3);
  //}
  
  hitangle_h=atan2(e_component, n_component);

  //if (inu==120 || inu==877 || inu==1927)
  //cout << "and now it's " << hitangle_h << "\n";

  //if (temp[0]>0.65 && temp[0]<0.75 && temp[1]>0.65 && temp[1]<0.75) {
  //cout << "hitangle_h is " << hitangle_h << "\n";
  //}  

  //  if (ilayer==3 && inu==382)
  //cout << "theta, phi for e_plane is " << dGetTheta(n_eplane) << " " << dGetPhi(n_eplane) << "\n";
  //if (inu==382)
  //cout << "n_component,e_component,hitangle_e, hitangle_h is " << n_component << " " << e_component << " " << hitangle_e << " " << hitangle_h << "\n";
  


  e_component=dDot(n_pol,n_eplane,3);
  
  // find component along the h-plane



  
  h_component=dDot(n_pol,n_hplane,3);


  
} //end void GetHitAngles



