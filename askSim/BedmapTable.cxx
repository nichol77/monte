///////////////////////////////////////////////////////////////////////////////
/////   BedmapTable.cxx
/////   Class to wrap round Bedmap table
///////////////////////////////////////////////////////////////////////////////
#include "BedmapTable.h"

#include <iostream>
#include <fstream>

using namespace std;

ClassImp(BedmapTable)



BedmapTable::BedmapTable() {
    //Default constructor
}

BedmapTable::BedmapTable(AskCons::BedmapType_t whichTable)
    :fType(whichTable),fArray(0) 
{
    LoadTable();
}

BedmapTable::~BedmapTable() {
    //Default destructor
    if(fArray) {
	
	delete [] fArray;
    }
}


void BedmapTable::LoadTable() {

    ifstream DataFile;
    switch(fType) {
	case AskCons::kBathymetry:
	    DataFile.open("geomData/bathy.asc");
	    break;
	case AskCons::kBedElevation:
	    DataFile.open("geomData/bedelev.asc");
	    break;
	case AskCons::kGroundBed:
	    DataFile.open("geomData/groundbed.asc");
	    break;
	case AskCons::kIceThickness:
	    DataFile.open("geomData/icethic.asc");
	    break;
	case AskCons::kSurfaceHeight:
	    DataFile.open("geomData/surface.asc");
	    break;
	case AskCons::kWaterDepth:
	    DataFile.open("geomData/water.asc");
	    break;
	default:
	    cerr << "Unknown BedmapType: " << fType << endl;
	    exit(1);
    }
    if(!DataFile) {
	cerr << "Couldn't open file" << endl;
	exit(1);
    }
    string tempBuf1;
    string tempBuf2;
    string tempBuf3;
    string tempBuf4;
    string tempBuf5;
    string tempBuf6;
    Int_t temp1,temp2,temp6;
    Float_t temp3,temp4,temp5;

    if(!(DataFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2
       >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
	 >> tempBuf5 >> temp5 >> tempBuf6 >> temp6)) {       
	cerr << "Bugger couldn't get header data" << endl;
 	exit(1);
    }
    if(tempBuf1 == string("ncols")) {
	fNumberCols=temp1;
    }
    if(tempBuf2 == string("nrows")) {
	fNumberRows=temp2;
    }
    if(tempBuf3 == string("xllcorner")) {
	fXLowerLeft=temp3;
    }
    if(tempBuf4 == string("yllcorner")) {
	fYLowerLeft=temp4;
    }
    if(tempBuf5 == string("cellsize")) {
	fCellSize=temp5;
    }
    if(tempBuf6 == string("NODATA_value")) {
	fNullValue=temp6;
    }
  
    fArray= new Float_t[fNumberRows*fNumberCols];

    fXUpperRight=fXLowerLeft+(fCellSize*fNumberCols);
    fYUpperRight=fYLowerLeft+(fCellSize*fNumberRows);
    
    Float_t theValue;
    for(int rowNum=0;rowNum<fNumberRows;rowNum++) {
	for(int colNum=0;colNum<fNumberCols;colNum++) {
	    if(!(DataFile >> theValue)) {		
		cout << "BedmapType: " << fType << endl;
		cout << fNumberCols << "\t" << fNumberRows << "\t"
		     << fXLowerLeft << "\t" << fXUpperRight << "\t"
		     << fCellSize << "\t" << fNullValue << endl;
		    
		cerr << "Bugger couldn't get data" << endl;
		exit(1);
	    }
	    fArray[index(rowNum,colNum)]=theValue;
	}
    }

    
}

Int_t BedmapTable::index(Int_t rowNum, Int_t colNum)
{
    if(rowNum>=fNumberRows || colNum>=fNumberCols) return -1;
    if(rowNum<0 || colNum<0) return -1;
    return colNum + fNumberCols*rowNum;    
}

void BedmapTable::getCellCentre(Int_t rowNum, Int_t colNum, Double_t &xCentre, Double_t &yCentre) {
    xCentre=double(colNum)*fCellSize+fXLowerLeft+fCellSize/2.;
    yCentre=double(rowNum)*fCellSize+fYLowerLeft+fCellSize/2.;    
}

Float_t BedmapTable::getValue(Double_t x, Double_t y, Int_t &goodFlag)
{
    //Will need to do something more clever
    y*=-1; //Silly hack as I don't understand the projection

    goodFlag=0;
    if(x<fXLowerLeft || x>fXUpperRight) return 0;
    if(y<fYLowerLeft || y>fYUpperRight) return 0;
 
    int colNum=int((x-fXLowerLeft)/fCellSize);
    int rowNum=int((y-fYLowerLeft)/fCellSize);
//    std::cout << "Here:\t" << rowNum << "\t" << colNum << std::endl;

    Double_t thisCellX,thisCellY;
    getCellCentre(rowNum,colNum,thisCellX,thisCellY);

    Int_t leftOrRight=0;
    if(x<thisCellX) leftOrRight=-1;
    else leftOrRight=1;
    
    Int_t upOrDown=0;
    if(y<thisCellY) upOrDown=-1;
    else upOrDown=1;

    Int_t indexes[4];
    Double_t xCentres[4];
    Double_t yCentres[4];
    Double_t values[4];

    Int_t rowNums[4]={rowNum,rowNum,rowNum+upOrDown,rowNum+upOrDown};
    Int_t colNums[4]={colNum,colNum+leftOrRight,colNum,colNum+leftOrRight};

    for(int i=0;i<4;i++) {
	indexes[i]=index(rowNums[i],colNums[i]);
	getCellCentre(rowNums[i],colNums[i],xCentres[i],yCentres[i]);
	if(indexes[i]!=-1) values[i]=fArray[indexes[i]];
    }

    Double_t mean=0;
    Double_t weight=0;
    for(int i=0;i<4;i++) {
	if(indexes[i]!=-1 && ((int)values[i]!=fNullValue)) {
	    Double_t temp=sqrt(pow(x-xCentres[i],2)+pow(y-yCentres[i],2));
	    mean+=values[i]*temp;
	    weight+=temp;
	}
    } 
    if(weight>0) {
	mean/=weight;
	goodFlag=1;
	return mean;
    }
    return 0;
   
}

    
