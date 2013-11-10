///////////////////////////////////////////////////////////////////////////////
/////   NeutrinoCrossSectionModel.h
/////  Namespace for holding the cross section related code and calls
//////////////////////////////////////////////////////////////////////////////

#ifndef NEUTRINOCROSSSECTIONMODEL_H
#define NEUTRINOCROSSSECTIONMODEL_H

#include "TObject.h"
#include "TRandom3.h"
#include "TF1.h"

#include "AskConventions.h"


namespace NeutrinoCrossSectionModel  {
   //class or namespace?
   
   //All lengths in metres unless otherwise stated
   //All energies in GeV unless otherwise stated


    enum WhichCrossection_t {
	kDefault    =  0,
	kRenoCteq6 = 1
    };
    
    void setWhichCrossection(NeutrinoCrossSectionModel::WhichCrossection_t whichCross);
    void setScaleFactor(Double_t scaleFactor);
         

    //General Functions
    Double_t getCrossSection(Double_t energy, AskCons::ParticleType_t particle);
    Double_t getInteractionLength(Double_t energy, AskCons::ParticleType_t particle, AskCons::MaterialType_t mat);
    Double_t getDistanceToInteraction(double intLength);
    

    //Reno CTEQ6 -- Add ref
    Double_t getRenoCteq6CrossSection(double energy) ;
    
    //Gandhi et al 1998 -- Add Ref
    Double_t funcCCCrossSec(Double_t *energy, Double_t *par);
    Double_t funcNCCrossSec(Double_t *energy, Double_t *par);
    Double_t funcTotCrossSec(Double_t *energy, Double_t *par);
    Double_t funcNuCCInteractionLength(Double_t *energy, Double_t *par);
    Double_t funcNuNCInteractionLength(Double_t *energy, Double_t *par);
    Double_t funcNuTotInteractionLength(Double_t *energy, Double_t *par);
    Double_t funcNuCCInteractionLengthInCM(Double_t *energy, Double_t *par);
    Double_t funcNuNCInteractionLengthInCM(Double_t *energy, Double_t *par);
    Double_t funcNuTotInteractionLengthInCM(Double_t *energy, Double_t *par);  


    //Utility Functions
    double intLengthFromProb(double *p, double *par);

};


#endif //NEUTRINOCROSSSECTIONMODEL_H
