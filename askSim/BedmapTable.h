///////////////////////////////////////////////////////////////////////////////
/////   BedmapTable.h
/////   Class to hold Bedmap data
/////   rjn@mps.ohio-state.edu
/////   4th May 2005
//////////////////////////////////////////////////////////////////////////////

#ifndef BEDMAPTABLE_H
#define BEDMAPTABLE_H

#include "TObject.h"

#include "AskConventions.h"
#include "AskGeom.h"

class BedmapTable : public TObject {
public:
    BedmapTable();    
    BedmapTable(AskCons::BedmapType_t whichTable);    
    ~BedmapTable();
   
    Float_t getValue(Double_t x, Double_t y, Int_t &goodFlag);

    Float_t fXLowerLeft;
    Float_t fYLowerLeft;
    Float_t fXUpperRight;
    Float_t fYUpperRight;
    
private:    
    Int_t index(Int_t rowNum, Int_t colNum);
    AskCons::BedmapType_t fType;
    void getCellCentre(Int_t rowNum, Int_t colNum, Double_t &xCentre, Double_t &yCentre);
    Int_t fNullValue;
    Int_t fNumberRows;
    Int_t fNumberCols;
    Float_t fCellSize;
    Float_t *fArray;//[fNumberRows*fNumberCols]
    void LoadTable();    
    ClassDef(BedmapTable,2);
};


#endif //BEDMAPTABLE_H
