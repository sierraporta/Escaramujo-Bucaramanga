#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <bitset>
#include <cstddef> 
using namespace std;

#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TNtuple.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TApplication.h"
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TROOT.h>
#include "TStyle.h"
#include "TTreeViewer.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>

const int gVerbosity   = 1;
const int kMaxLine     = 1024;
const int kDataCharacters = 72; // This number given by the format from QuarkNet
const int kClockTSize = 8;
const int kNumTotalEdges = 8; //4 channels with 2 edges each (?)
const int kNumChannels = 4; //Given by the acquisition board
const int kNumEdgeBytes = 3; //2 digits for edges plus 1 for space

const float kMaxTime = 30e3;
const float kMinTime = 0;

const int kMaxNHitsTDC = 16;

const double kClockRes = 40;    // ns/count
const double kBoardRes = kClockRes/32;  // ns/count (1.25 for 6000, 0.75 for 5000 series

struct event_t
{
  vector<vector<Double_t> > tdcLE;
  vector<vector<Double_t> > tdcTE;  
  double tGMT;
  int dGMT;
  int okGPS;
  int nSat;
};

struct run_t
{
  vector<int> chHits;
  vector<event_t> event;
};


void printHelp(){
  cout << "\nUsage:\n";
  cout << " * Save root file:\n";
  cout << "    ./qNet2root.exe <input qNet file> <output root file>\n";
}


std::bitset<4> hex2bitset(const string& hexStr)
{
  stringstream hexStrSS;
  hexStrSS << hex << hexStr;
  unsigned n;
  hexStrSS >> n;
  std::bitset<4> b(n);
  return b;
}

int hexData2Time(string edgeStrLSB, string edgeStrMSB){
  int intedHex;
  edgeStrMSB.append(edgeStrLSB);
  stringstream hexStrSS;
  hexStrSS << hex << edgeStrMSB;
  unsigned n;
  hexStrSS >> n;
  std::bitset<5> b(n);
  return intedHex = b.to_ulong();
}

bool evStart(string hexStr)
{  // decides if the event start maker is present. MSB in MSByte
  std::bitset<4> b = hex2bitset(hexStr);
  return (b[3]==1);
}    

bool hit(const string& hexStr)
{ 
  /* decides if there is a hit or not in an hex number 
     (if the bit 1 is on of MSByte) */
  bool hit = 0;
  std::bitset<4> b = hex2bitset(hexStr);
  if(b[1]==1) 
    hit = 1;
  return hit;
}

double clockStr2double(const string& clockStr)
{ // transform the clock string into a time in [ns]
  stringstream clockTSS;
// 	  cout << " >>" << clockStr;
  clockTSS << hex << clockStr;
  unsigned unsigCl;
  clockTSS >> unsigCl;
  std::bitset<32> bCl(unsigCl);
// 	  cout << " " << bCl.to_ulong() << " ";
  double clockTi = bCl.to_ulong()*kClockRes;
//   cout << " " << std::setprecision(0) << std::fixed << clockTi << "<< ";
  return clockTi;
}

bool saneLine(string lineS)
{
  bool saneLineRet = 1;

  //Check if the number of characters makes sense
  if(lineS.size()<kDataCharacters)
  {
    if(gVerbosity>2)
      cout<<"SANE LINE excluding "<<lineS<<" because it has only "<<lineS.size()<<" characters, as opposed to the expected "<<kDataCharacters<<endl;
    return false;
  }

  { // check if there is trash in the line //
    size_t found = lineS.find_first_not_of(" 0123456789ABCDEFV\t\v\r\n#+-."); 
    //V marks an error in the GPS; + is needed for final column
    if(found != string::npos || lineS.empty())
    {
      saneLineRet = 0;              
      if(gVerbosity>2)
	cout<<"SANE LINE Excluding "<<lineS<<" because of bad characters"<<endl;
    }
  } 
  
  {  // check if the clock tick is an hex number 
    string clockStr = lineS.substr(0,kClockTSize);
    size_t found = clockStr.find_first_not_of("0123456789ABCDEF");
    if(found != string::npos){
      saneLineRet = 0;
      if(gVerbosity>2) 
	cout << "SANE LINE Clock Tick with characters other than hex\n";
    } 
  }
  
  { // check if there are spaces between the data
    //int kFirstSpace = 8; //This is unnecessary; equal to kClockTSize
    bool pars = 0;
    for(int i=0;i<8;++i){
      std::locale loc;
      //char s = lineS[kFirstSpace+3*i];
      char s = lineS[kClockTSize+kNumEdgeBytes*i];
      if(! std::isspace(s,loc) ) pars=1;
    }
    if(pars)
    {
      saneLineRet = 0;
      if(gVerbosity>2)
	cout<<"SANE LINE excluding "<<lineS<<" because of bad spacing"<<endl;
    } 
  }  
  /*Probably we should check the spacing of the remaining columns, 
    but for now we will assume that they are all the same, 
    since so far so good*/
//   cout << "saneLineRet = " << saneLineRet << endl;
  if(saneLineRet) return true;
  else return false;
}
    
    

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    cerr <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

int fillEvInfo(const string lineS, event_t& ev)
{
  stringstream iss(lineS);
  string dummy="";
  iss>>dummy; //Initial clock
  for(int i=0; i<kNumTotalEdges; ++i)
    iss>>dummy; //Clocks for edges
  iss>>dummy; //Clock for 1-second pulser
  iss>>ev.tGMT; //GPS time
  iss>>ev.dGMT; //GPS date
  string GPS_label="";
  iss>>GPS_label;
  ev.okGPS = (int)(GPS_label=="A");
  iss>>ev.nSat;
  return 0;
}

bool readFile(const char *inFile, run_t &run){
  
  vector<event_t> &evV = run.event;
  run.chHits.resize(4);
  
  if(gVerbosity>0){
    cout << "Reading input file: " << inFile << endl;
  }
  
  ifstream in(inFile);
  
  if(!in.is_open()){
    cerr << "\nCan't open file: " << inFile << endl << endl;
    return false;
  }
  
  unsigned int n=0;
  unsigned int nErr=0;
  unsigned int ln=0;
  
  // aks the first line
  char line[kMaxLine];
  in.getline(line,kMaxLine);
  ++ln;
  string lineS(line,kDataCharacters);
  
  while(!saneLine(lineS) && !(in.eof()))
  {
    in.getline(line,kMaxLine);
    ++ln;
    lineS.assign(line,kDataCharacters);
  } 
  if(in.eof())
  {
    cout<<"ERROR: could not find the first line containing data"<<endl;
    return -1;
  }
  if(gVerbosity>2) 
    cout << ln << ": " << lineS << "\tFirst line" << endl;
    
  while(!in.eof()){ // start to read file
    
//     cout << "evV.size()=" << evV.size() << endl;
//     if(evV.size()>20) break;

    const string evStartString = lineS.substr(kClockTSize+1,1);
    bool bEvStart = evStart(evStartString);  // test the event start flag
    if(bEvStart)
    {  // event start
      ++n;
      if(gVerbosity>2) 
	cout << "\t\t * * * event Start n = " << n << "\t\t * * * " << endl;
      event_t ev;
      ev.tdcLE.resize(4);
      ev.tdcTE.resize(4);
      if(fillEvInfo(lineS,ev))
      {
	cout<<"ERROR: failed at filling event info"<<endl;
	return -1;
      }
      bEvStart = 0; // put the event start flag OFF
      bool evSane = 1;
      const string clockStrEvStart = lineS.substr(0,kClockTSize);
      const double clockTEvStart = clockStr2double(clockStrEvStart);
      do
      {
	
	const string clockStr = lineS.substr(0,kClockTSize); //Not start time
	const double clockT = clockStr2double(clockStr) - clockTEvStart;
	if(gVerbosity>2) 
	  cout << ln << ": " << lineS << "\tIN of event" << endl;
	
	for(int i=0;i<kNumTotalEdges;++i)
	{
	  const string edgeStrMSB=lineS.substr(kClockTSize+1+i*kNumEdgeBytes,1);
	  const string edgeStrLSB=lineS.substr(kClockTSize+2+i*kNumEdgeBytes,1);
	  if(hit(edgeStrMSB))
	  {
	    const Double_t t=
	      hexData2Time(edgeStrLSB,edgeStrMSB)* kBoardRes + clockT; 
	    if(t<kMinTime || t>kMaxTime) 
	    {
	      if(gVerbosity>2)
		cout<<"EV SANE throwing out ev because there are hits beyond "
		    <<(int)(kMaxTime/1000.)<<"us"<<endl<<lineS<<endl;
	      evSane = 0;
	    }
	    const unsigned int chTag = (unsigned int) ((i/2)%kNumChannels);
	    const int edgeTypeInt = (i%2);
	    if     (edgeTypeInt==0) 
	      ev.tdcLE[(int)chTag].push_back(t);
	    else if(edgeTypeInt==1) 
	      ev.tdcTE[(int)chTag].push_back(t);
	    else 
	    {
	      if(gVerbosity>2)
		cout<<"EV SANE edge is neither leading nor trailing"<<endl;
	      evSane = 0;
	    }
	    if(gVerbosity>2)
	    {
	      cout << ln << ": "; 
	      cout << edgeStrMSB << edgeStrLSB << " ";
	      cout << "ch" << chTag << "_" << edgeTypeInt << " " 
		   << std::setprecision(2) << std::fixed << t << " " 
		   << clockStr << " " << clockStrEvStart;
	      cout << " " << clockStr2double(clockStr) << " " 
		   << clockStr2double(clockStrEvStart) << "\n";
	    }
	  }
	}
	/*char line[kMaxLine];
	  in.getline(line,kMaxLine);*/
	string line="";
	std::getline(in,line);
	if(line.size()<kDataCharacters)
	{
	  //This line is too small, probably it's not data
	  if(line=="CD" || line=="")
	  {
	    ++nErr;
	    if(gVerbosity>2)
	      cout<<"LINE SKIP: Removing small line "<<endl<<line<<endl;
	    break;
	  }
	  else
	  {
	    cout<<"ERROR: line "<<line<<" is not acceptable"<<endl;
	    return -1;
	  }
	}
	++ln;
	lineS.clear();
	lineS.assign(line.c_str(),kDataCharacters);
	if(!saneLine(lineS))
	{
	  ++nErr;
	  if(gVerbosity>2)
	    cout<<"LINE SKIP error in line"<<endl;
	  break;
	}
	
	const string evStartStringAgain = lineS.substr(kClockTSize+1,1);
	bEvStart = evStart(evStartStringAgain);  // test the event start flag
	if(bEvStart)
	{
	  for(unsigned int c=0;c<kNumChannels;++c)
	  {
	    run.chHits[c] += ev.tdcLE[c].size();
	    if(ev.tdcLE[c].size()>kMaxNHitsTDC-2) 
	    {
	      if(gVerbosity>2)
		cout<<"EV SANE throwing out event because there are more than "
		    <<kMaxNHitsTDC-2<<" hits"<<endl<<lineS<<endl;
	      evSane = 0;
	    }
	    if(ev.tdcLE[c].size() != ev.tdcTE[c].size() ) 
	    {
	      if(gVerbosity>2)
		cout<<"EV SANE diff nums of leading and trailing edges"<<endl;
//	      cout<<"DEBUG: "<<ev.tdcLE[c].size()<<" "<<ev.tdcTE[c].size()<<" "<<c<<endl;
	      evSane = 0;
	    }
	  }
	  if(!evSane)
	  {
	    ++nErr;
	    if(gVerbosity>2)
	      cout<<"LINE SKIP line excluded because of event error:\n  "
		  <<lineS<<endl;
	    continue;
	  } 
	  evV.push_back(ev);
	}
	
      } while(!bEvStart); // do this while there the new event flag is off
      
      
    } // event ends
    else //I believe this is not really doing anything other than getting the next line
    {
      // aks the first line
      char line[kMaxLine];
      in.getline(line,kMaxLine);
      ++ln;
      lineS.assign(line,kDataCharacters);
      if(gVerbosity>2) 
	cout << ln << ": " << lineS << "\tOUT of event" << endl;
      /*if(!saneLine(lineS)) //I'm pretty sure this is unnecessary
      {
	if(gVerbosity>2)
	  cout<<"DEBUG: unsane line "<<lineS<<endl;
	continue;
	}*/
    }
    
    
  } // end of while that reads the file until de eof
  
  in.close();
  
  if(gVerbosity>0){
    cout << "=============================\n";
    cout << "Events read:\t" << n
    << "\nParsing errors:\t" << nErr
    << "\nGood events:\t"    << evV.size() << endl;

    cout << "\nChannels present (TDC)\n";
    cout << " CH\tCounts\n";
    for(unsigned int c=0;c<kNumChannels;++c){
      cout << " " << c << "\t" << run.chHits[c] << endl;
    }
//     for(std::map<int, int>::iterator it = run.chMap.begin(); it != run.chMap.end(); it++){
//       cout << " " << it->first << "\t" << it->second << endl;
//     }
    cout << "=============================\n";
  }

  return true;
}  

void data2tree(const run_t &run, TTree *dataTree){

  const vector<event_t> &evV = run.event;
  const unsigned int nCh = kNumChannels;
  
  Double_t **tdcLEV = new Double_t *[nCh];
  Double_t **tdcTEV = new Double_t *[nCh];
  Int_t *tdcSizeV = new Int_t[nCh];
  Int_t GMTDate = -1, GPSOK = -1, NumSat = -1;
  Double_t GMTTime = -1;
  
  //Channel name map
  map<string, string> channel2chName;
  vector<string> tdcLEName;
  vector<string> tdcTEName;
  vector<string> tdcSizeName;
  
  for(unsigned int c=0;c<nCh;++c){

    ostringstream varName;
    varName << "TDC_LE_" << c;
    tdcLEName.push_back(varName.str());
    
    varName.str("");
    varName << "TDC_TE_" << c;
    tdcTEName.push_back(varName.str());
    
    ostringstream varSizeName;
    varSizeName << "sizeTDC_" << c;
    tdcSizeName.push_back(varSizeName.str());
    
  }
  
  cout << "\nInitialize TTree structure\n";
  //Init tree structure
  for(unsigned int c=0;c<nCh;++c){
    tdcLEV[c] = new  Double_t[kMaxNHitsTDC];
    tdcTEV[c] = new  Double_t[kMaxNHitsTDC];
    dataTree->Branch(tdcSizeName[c].c_str(), &(tdcSizeV[c]),(tdcSizeName[c]+"/I").c_str() );
    dataTree->Branch(tdcLEName[c].c_str(),tdcLEV[c],(tdcLEName[c]+"["+tdcSizeName[c]+"]/D").c_str() );
    dataTree->Branch(tdcTEName[c].c_str(),tdcTEV[c],(tdcTEName[c]+"["+tdcSizeName[c]+"]/D").c_str() );
  }
  dataTree->Branch("tGPS",&GMTTime);
  dataTree->Branch("dGPS",&GMTDate);
  dataTree->Branch("sGPS",&GPSOK);
  dataTree->Branch("nGPS",&NumSat);

  //Fill tree
  cout << "\nFilling TTree \n";
  for(unsigned int i=0;i<evV.size();++i){ // loop on run events
    
//     cout << "i=" << i << endl;
//     if(i>4000 && i<5000) continue;

    for(unsigned int c=0;c<nCh;++c){// loop on tdc channels
      
      const unsigned int nEntr = evV[i].tdcLE[c].size();
      tdcSizeV[c] = nEntr;
//       cout << "ch=" << c << "\tnEntr=" << nEntr << endl;
      for(unsigned int e=0;e<nEntr;++e){
	
	tdcLEV[c][e]         = evV[i].tdcLE[c][e];
// 	cout << "llega\n";
	tdcTEV[c][e]         = evV[i].tdcTE[c][e];
	
// 	cout << i << "," << c << "," << e << " = " << tdcV[c][e] << " " << tdcEdgeTypeV[c][e] << endl;
      }
        
      
//       std::map<int, vector<double> >::const_iterator tdcIt=evV[i].tdc.find(tdcCh[c]);
//       std::map<int, vector<int> >::const_iterator tdcEdgeTypeIt=evV[i].tdcEdgeType.find(tdcCh[c]);
//       if(tdcIt==evV[i].tdc.end() || tdcEdgeTypeIt==evV[i].tdcEdgeType.end()){
//         tdcSizeV[c]=0;
//       }
//       else{
//         const unsigned int nEntr = tdcIt->second.size();
//         tdcSizeV[c] = nEntr;
//         for(unsigned int e=0;e<nEntr;++e){
//           tdcV[c][e]=(tdcIt->second)[e];
//           tdcEdgeTypeV[c][e]=(tdcEdgeTypeIt->second)[e];
//         }
//       }
      
    }
    GMTTime = evV[i].tGMT;
    GMTDate = evV[i].dGMT;
    GPSOK = evV[i].okGPS;
    NumSat = evV[i].nSat;

    dataTree->Fill();
  }

  //dataTree->Print(); 
}

int main(int argc,char *argv[]){

  if( argc!=3 ){
    printHelp();
    return 0;
  }
  
  if(!fileExist(argv[1])){
    return 1;
  }

  run_t run;
  readFile(argv[1], run);

  TFile *outF = new TFile(argv[2],"RECREATE");
  TTree *dataTree = new TTree("dataTree","dataTree");
  data2tree(run, dataTree);
  cout << "\nWriting TTree on disc\n";
  dataTree->Write();
  cout << "\n Closing ROOT file \n";
//   outF->Close();  
  delete outF;
  cout << "\n * * * All done * * * \n";
  
 return 0;
}
