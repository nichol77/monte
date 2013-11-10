#ifndef ASKGEOM_H
#define ASKGEOM_H


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AskGeom                                                            //
//                                                                      //
// Encapsulate geometry and coordinate related thingies                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#include "TMath.h"
#include "TVector3.h"
#include <iostream>

namespace crust2 {
    const int N_LON=180; // number of bins in longitude for crust 2.0
    const int N_LAT=90;  // number of bins in latitude
    const double MAX_THETA=180; // maximum value of theta in degrees in polar coordinates
}

namespace AskGeom {
    const double R_EARTH=6.378137E6;
    const double GEOID_MAX=6.378137E6; // parameters of geoid model
    const double GEOID_MIN=6.356752E6;
    
/*     Double_t getSinFromCos(Double_t cosTheta); */
/*     Double_t getCosThetaBewteenVectors(TVector3 &vec1, TVector3 &vec2); */
/*     void fastRotateX(TVector3 &veccy, Double_t sinAng, Double_t cosAng); */
/*     void fastRotateY(TVector3 &veccy, Double_t sinAng, Double_t cosAng); */
/*     void fastRotateZ(TVector3 &veccy, Double_t sinAng, Double_t cosAng); */


    inline Double_t getSinFromCos(Double_t cosTheta) {
      return TMath::Sqrt(1. - cosTheta*cosTheta);
    }

    inline Double_t getCosThetaBewteenVectors(TVector3 &vec1, TVector3 &vec2) {
      Double_t ptot2 = vec1.Mag2()*vec2.Mag2();
      if(ptot2 <= 0) {
	return 0.0;
      } else {
      Double_t arg = vec1.Dot(vec2)/TMath::Sqrt(ptot2);
      if(arg >  1.0) arg =  1.0;
      if(arg < -1.0) arg = -1.0;
      return arg;
      }
    }
    
    inline void fastRotateX(TVector3 &veccy, Double_t sinAng, Double_t cosAng) {
      Double_t zz=veccy.Z();
      Double_t yy=veccy.Y();
      veccy.SetY(cosAng*yy - sinAng*zz);
      veccy.SetZ(sinAng*yy + cosAng*zz);
    }
    
   inline  void fastRotateY(TVector3 &veccy, Double_t sinAng, Double_t cosAng) {
      Double_t zz=veccy.Z();
      Double_t xx=veccy.X();
      veccy.SetZ(cosAng*zz - sinAng*xx);
      veccy.SetX(sinAng*zz + cosAng*xx);
    }
    
    inline void fastRotateZ(TVector3 &veccy, Double_t sinAng, Double_t cosAng) {
      Double_t yy=veccy.Y();
      Double_t xx=veccy.X();
      veccy.SetX(cosAng*xx - sinAng*yy);
      veccy.SetY(sinAng*xx + cosAng*yy);
    }
    
    
    //Need to check these to make sure that they use the correct variables
    inline Double_t getTheta(Int_t lat)  {        
      return (((Double_t)lat+0.5)/(Double_t)crust2::N_LAT*crust2::MAX_THETA)*TMath::DegToRad(); 
    }
    inline Double_t getPhi(Int_t lon) { 
      return (double)(-1*((double)lon+0.5)+(double)crust2::N_LON)*2*TMath::Pi()/(double)crust2::N_LON-TMath::Pi()/2;
    }
    inline Int_t getLat(Double_t theta) { 
       
      return (int)(90.-((theta*TMath::RadToDeg()))); 
    }
    
    inline Int_t getLon(Double_t phi){ 
      //Need to fix this somehow
      double phi_deg = phi*TMath::RadToDeg();
      if (phi_deg>270)
	phi_deg-=360;	
      return Int_t(90.-phi_deg);
    }
    
    inline Double_t getPhi(Double_t p[3]){    
      // returns phi between 0 and 2pi.
	double pt=0;
	double phi=0;
	pt=sqrt(p[0]*p[0]+p[1]*p[1]);    
	if (pt==0)
	  return 0.;
	else if (pt!=0) {
	  if (p[1]/pt>1 || p[1]/pt<-1) {
		std::cerr << "Error in getPhi. \n";
		return 0;
	    }
	    phi=asin(p[1]/pt);
	}
	if (p[1]<0. && p[0]>0) phi += 2*TMath::Pi();
	else if (phi>0 && p[0]<0.) phi = TMath::Pi() - phi;
	else if (phi<0 && p[0]<0.) phi = -(TMath::Pi()+phi)+2*TMath::Pi();
	return phi;
    }

     inline Double_t getPhi(TVector3 &thePos) {
       double p[3]={thePos.X(),thePos.Y(),thePos.Z()};
       return getPhi(p);
       //return thePos.Theta();
    }
    
    
     inline Double_t getTheta(Double_t p[3]) {
      double pz,pt;
      double tantheta1=0;
      double theta=0;
      
      pz=p[2];
      pt=sqrt(p[0]*p[0]+p[1]*p[1]);
      tantheta1=pt/pz;
      theta=atan(tantheta1);
      
      if (pz<0)
	theta += TMath::Pi();  
      return theta;  
    }


     inline Double_t getTheta(TVector3 &thePos) {
       double p[3]={thePos.X(),thePos.Y(),thePos.Z()};
       thePos.GetXYZ(p);
       return getTheta(p);
    }

      inline void getLonLat(Double_t p[3],Int_t& lon,Int_t& lat) {
	lon=AskGeom::getLon(AskGeom::getPhi(p));
	lat=AskGeom::getLat(AskGeom::getTheta(p));
/* 	if(p[2]<0) { */
/* 	   lat*=-1; */
/* 	} */
    }

     inline void getLonLat(TVector3 &thePos,Int_t& lon,Int_t& lat) {
	lon=AskGeom::getLon(AskGeom::getPhi(thePos));
	lat=AskGeom::getLat(AskGeom::getTheta(thePos));
    }

    inline Double_t getGeoidFromCosTheta(Double_t c) {
      return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-
					     (GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*c*c);    
    }

    inline Double_t getGeoid(Double_t p[3]) {
      Double_t c=cos(AskGeom::getTheta(p));
      return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-
					     (GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*c*c);    
    }
    inline Double_t getGeoid(TVector3 &thePos) {
      Double_t p[3]={thePos.X(),thePos.Y(),thePos.Z()};
      return getGeoid(p);
    }


    inline Double_t getGeoid(Int_t lat) {
	return GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2)-(pow(GEOID_MIN,2)-pow(GEOID_MAX,2))*pow(cos(AskGeom::getTheta(lat)),2));    
    }

    inline Double_t getGeoidFromTheta(Double_t theta) {
      Double_t c=cos(theta);
      return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-
					     (GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*c*c);    
  
    }

    //Crust 2 thingies
    inline Int_t getILat(Double_t theta) {
        return (int)((theta*TMath::RadToDeg())/2.);
    }

    inline Int_t getILon(Double_t phi){
        double phi_deg = phi*TMath::RadToDeg();
        if (phi_deg>270)
            phi_deg-=360;
        return (int)((360.*0.75-phi_deg)*180./360.);
    }

    inline void getILonILat(Double_t p[3],Int_t& lon,Int_t& lat) {
        lon=getILon(getPhi(p));
        lat=getILat(getTheta(p));
    }



}


#endif //ASKGEOM_H
