// auto-veto.cc
// Runs during production, creates ROOT files of MJD veto data.
// C. Wiseman, A. Lopez, 10/24/16.
//
// NOTE: The scans are split across a few different loops over the events in the run.
// This is done to increase the flexibility of the code, since it checks many different
// quantities.  The performance hit should be minimal, since the size of the veto trees
// is relatively small.

#include <iostream>
#include <fstream>
#include "TTreeReader.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TFile.h"
#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

vector<int> VetoThreshFinder(TChain *vetoChain, string outputDir, bool makePlots=false);
void ProcessVetoData(TChain *vetoChain, vector<int> thresholds, string outputDir, bool errorCheckOnly=false);

void SetCardNumbers(int runNum, int &card1, int &card2);
int FindPanelThreshold(TH1D *qdcHist, int threshVal, int panel, int runNum);
double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler);
int PanelMap(int i, int runNum=0);
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry, vector<int> &ErrorVec);
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry);

int main(int argc, char** argv)
{
	// get command line args
	if (argc < 1) {
		cout << "Usage: ./auto-veto [run number]"
			  << "[-d (optional: draws QDC & multiplicity plots)]"
			  << "[-e (optional: error check only)]\n";
		return 1;
	}
	int run = stoi(argv[1]);
	if (run > 60000000 && run < 70000000) {
		cout << "Veto data not present in Module 2 runs.  Exiting ...\n";
		return 1;
	}
	bool draw = false, errorCheckOnly = false;
	string opt1 = "", opt2 = "";
	if (argc > 2) opt1 = argv[2];
	if (argc > 3) opt2 = argv[3];
	if (argc > 2 && (opt1 == "-d" || opt2 == "-d"))
		draw = true;
	if (argc > 2 && (opt1 == "-e" || opt2 == "-e"))
		errorCheckOnly = true;

	// output file directory
	string outputDir = "./output";

	// Only get the run path (so we can run with veto-only runs too)
	GATDataSet ds;
	string runPath = ds.GetPathToRun(run,GATDataSet::kBuilt);
	TChain *vetoChain = new TChain("VetoTree");

	// Verify that the veto data exists in the given run
	if (!vetoChain->Add(runPath.c_str())){
		cout << "File doesn't exist.  Exiting ...\n";
		return 1;
	}

	printf("\n========= Processing run %i ... %lli entries. =========\n",run,vetoChain->GetEntries());
	cout << "Path: " << runPath << endl;

	// Find the QDC pedestal location in each channel.
	// Set a software threshold value above this location,
	// and optionally output plots that confirm this choice.
	vector<int> thresholds = VetoThreshFinder(vetoChain, outputDir, draw);

	// Check for data quality errors,
	// tag muon and LED events in veto data,
	// and output a ROOT file for further analysis.
	ProcessVetoData(vetoChain,thresholds,outputDir,errorCheckOnly);

	printf("=================== Done processing. ====================\n\n");
	return 0;
}


vector<int> VetoThreshFinder(TChain *vetoChain, string outputDir, bool makePlots)
{
	// format: (panel 1) (threshold 1) (panel 2) (threshold 2) ...
	vector<int> thresholds;
	int threshVal = 35;	// how many QDC above the pedestal we set the threshold at

	long vEntries = vetoChain->GetEntries();
	TTreeReader reader(vetoChain);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");
	reader.Next();
	int runNum = vRun->GetRunNumber();
	reader.SetTree(vetoChain);  // resets the reader

	gStyle->SetOptStat(0);
	int bins=500, lower=0, upper=500;
	TH1D *hLowQDC[32];
	TH1D *hFullQDC[32];
	char hname[50];
	for (int i = 0; i < 32; i++) {
		sprintf(hname,"hLowQDC%d",i);
		hLowQDC[i] = new TH1D(hname,hname,bins,lower,upper);
		sprintf(hname,"hFullQDC%d",i);
		hFullQDC[i] = new TH1D(hname,hname,420,0,4200);
	}
	sprintf(hname,"Run %i Hit Multiplicity",runNum);
	TH1D *hMultip = new TH1D("hMultip",hname,32,0,32);

	// Set all thresholds to 1, causing all entries to have a multiplicity of 32
	int def[32] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	long skippedEvents = 0, prevGoodEntry=0;

	// MJVetoEvent variables, with run-based card numbers
	int card1=0, card2=0;
	SetCardNumbers(runNum,card1,card2);
	MJVetoEvent veto(card1,card2);
	MJVetoEvent prev, first;
	cout << "QDC 1 in slot " << card1 << ", QDC 2 in slot " << card2 << endl;

	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(def);
    	veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (CheckEventErrors(veto,prev,first,prevGoodEntry))
		{
    		skippedEvents++;
			// do the end of event resets before continuing
			prev = veto;
			prevGoodEntry = i;
			veto.Clear();
    		continue;
    	}
    	for (int q = 0; q < 32; q++) {
    		hLowQDC[q]->Fill(veto.GetQDC(q));
    		hFullQDC[q]->Fill(veto.GetQDC(q));
		}
		// save previous entries for the event error check
		prev = veto;
		prevGoodEntry = i;
		veto.Clear();
	}
	if (skippedEvents > 0) printf("VetoThreshFinder skipped %li of %li entries.\n",skippedEvents,vEntries);

	int thresh[32] = {9999};
	for (int i = 0; i < 32; i++)
	{
		thresh[i] = FindPanelThreshold(hLowQDC[i],threshVal,i,runNum);
		thresholds.push_back(i);
		thresholds.push_back(thresh[i]);
	}
	// cout << "vtf thresholds: " << endl;
	// for (int i = 0; i < 32; i++)
		// cout << i << " " << thresh[i] << endl;

	// re-scan with the found thresholds to make a multiplicity plot
	reader.SetTree(vetoChain);  // resets the reader
	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();
		veto.SetSWThresh(thresh);
    	veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (CheckEventErrors(veto,prev,first,prevGoodEntry))
		{
    		skippedEvents++;
			// do the end of event resets before continuing
			prev = veto;
			prevGoodEntry = i;
			veto.Clear();
    		continue;
    	}
		hMultip->Fill(veto.GetMultip());
		// save previous entries for the event error check
		prev = veto;
		prevGoodEntry = i;
		veto.Clear();
	}
	if (makePlots)
	{
		TCanvas *c1 = new TCanvas("c1","full QDC",1600,1200);
		c1->Divide(8,4,0,0);
		for (int i=0; i<32; i++)
		{
			c1->cd(i+1);
			TVirtualPad *vpad1 = c1->cd(i+1); vpad1->SetLogy();
			hFullQDC[i]->Draw();
		}
		TCanvas *c2 = new TCanvas("c2","QDC thresholds",1600,1200);
		c2->Divide(8,4,0,0);
		for (int i=0; i<32; i++)
		{
			c2->cd(i+1);
			TVirtualPad *vpad2 = c2->cd(i+1); vpad2->SetLogy();
			hLowQDC[i]->GetXaxis()->SetRange(lower,upper);
			hLowQDC[i]->Draw();
			double ymax = hLowQDC[i]->GetMaximum();
			TLine *line = new TLine(thresh[i],0,thresh[i],ymax+10);
			line->SetLineColor(kRed);
			line->SetLineWidth(2.0);
			line->Draw();
		}
		TCanvas *c3 = new TCanvas("c3","multiplicity",800,600);
		c3->cd();
		c3->SetLogy();
		hMultip->Draw();

		char plotName[200];
		sprintf(plotName,"%s/veto-%i-qdc.pdf",outputDir.c_str(),runNum);
		c1->Print(plotName);

		sprintf(plotName,"%s/veto-%i-qdcThresh.pdf",outputDir.c_str(),runNum);
		c2->Print(plotName);

		sprintf(plotName,"%s/veto-%i-multip.pdf",outputDir.c_str(),runNum);
		c3->Print(plotName);

	}
	return thresholds;
}

void ProcessVetoData(TChain *vetoChain, vector<int> thresholds, string outputDir, bool errorCheckOnly)
{
	// QDC software threshold (obtained from VetoThreshFinder)
	int swThresh[32] = {0};
	for (int i = 0; i < (int)thresholds.size(); i+=2)
		swThresh[thresholds[i]] = thresholds[i+1];
	// cout << "ProcessVetoData is using these QDC thresholds: " << endl;
	// for (int i = 0; i < 32; i++)
	// 	cout << i << " " << swThresh[i] << endl;

	// LED variables
	int LEDMultipThreshold = 5;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	int LEDSimpleThreshold = 10;  // used when LED frequency measurement is bad.
	int highestMultip=0, multipThreshold=0;
	double LEDWindow = 0.1;
	double LEDfreq=0, LEDperiod=0, LEDrms=0;
	bool badLEDFreq=false;
	bool IsLED = false, IsLEDPrev = false;
	bool LEDTurnedOff = false;
	int simpleLEDCount=0;
	bool useSimpleThreshold=false;
	bool firstLED=false;

	// Error checks
	const int nErrs = 29; // error 0 is unused
	int SeriousErrorCount = 0;
	int TotalErrorCount = 0;
	vector<double> EntryTime;
	vector<double> EntryNum;
	vector<bool> BadScalers;
	vector<int> EventError(nErrs), ErrorCount(nErrs);
	bool badEvent = false;
	// Specify which error types to print during the loop over events (event-level errors)
	vector<int> SeriousErrors = {1, 13, 14, 18, 19, 20, 21, 22, 23, 24};

	// muon ID variables
	bool TimeCut = true;	// actually an LED cut
	bool EnergyCut = false;
	int PlaneHitCount=0;
	vector<int> CoinType(32), CutType(32), PlaneHits(32), PlaneTrue(32);

	// time variables
	double SBCOffset=0;
	bool ApproxTime = false;
	bool foundScalerJump = false;	// this is related to ApproxTime but should be kept separate.
	double deltaTadj = 0;
	double timeSBC=0, prevtimeSBC=0;
	long skippedEvents=0;
	long start=0, stop=0;
	double duration=0, livetime=0, xTime=0, x_deltaT=0, x_LEDDeltaT=0;
	double xTimePrev=0, x_deltaTPrev=0, xTimePrevLED=0, xTimePrevLEDSimple=0, TSdifference=0;
	double prevGoodTime=0, firstGoodScaler=0;
	int prevGoodEntry=0;
	bool foundFirst=false, foundFirstScaler=false;

	// initialize input data
	long vEntries = vetoChain->GetEntries();
	TTreeReader reader(vetoChain);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");
	reader.Next();
	int runNum = vRun->GetRunNumber();
	start = vRun->GetStartTime();
	stop = vRun->GetStopTime();
	duration = (double)(stop - start);
	reader.SetTree(vetoChain);  // resets the reader

	// MJVetoEvent variables, with run-based card numbers
	int card1=0, card2=0;
	SetCardNumbers(runNum,card1,card2);
	MJVetoEvent veto(card1,card2);
	MJVetoEvent first, prev, last, prevLED, out;

	// initialize output file
	char outputFile[200];
	sprintf(outputFile,"%s/veto_run%i.root",outputDir.c_str(),runNum);
	TFile *RootFile = new TFile(outputFile, "RECREATE");
	TTree *vetoTree = new TTree("vetoTree","MJD Veto Events");
	// event info
	vetoTree->Branch("run",&runNum);
	vetoTree->Branch("vetoEvent","MJVetoEvent",&out,32000,1);	// TODO: organize these better
	// LED variables
	vetoTree->Branch("LEDfreq",&LEDfreq);
	vetoTree->Branch("LEDrms",&LEDrms);
	vetoTree->Branch("multipThreshold",&multipThreshold);
	vetoTree->Branch("highestMultip",&highestMultip);
	vetoTree->Branch("LEDWindow",&LEDWindow);
	vetoTree->Branch("LEDMultipThreshold",&LEDMultipThreshold);
	vetoTree->Branch("LEDSimpleThreshold",&LEDSimpleThreshold); // TODO: add a bool signifiying this is in use
	vetoTree->Branch("useSimpleThreshold",&useSimpleThreshold);
	// time variables
	vetoTree->Branch("start",&start,"start/L");
	vetoTree->Branch("stop",&stop,"stop/L");
	vetoTree->Branch("duration",&duration);
	vetoTree->Branch("livetime",&livetime);
	vetoTree->Branch("xTime",&xTime);
	vetoTree->Branch("timeSBC",&timeSBC);
	vetoTree->Branch("x_deltaT",&x_deltaT);
	vetoTree->Branch("x_LEDDeltaT",&x_LEDDeltaT);
	// muon ID variables
	vetoTree->Branch("CoinType",&CoinType);
	vetoTree->Branch("CutType",&CutType);
	vetoTree->Branch("PlaneHits",&PlaneHits);
	vetoTree->Branch("PlaneTrue",&PlaneTrue);
	vetoTree->Branch("PlaneHitCount",&PlaneHitCount);
	// error variables
	vetoTree->Branch("badEvent",&badEvent);
	vetoTree->Branch("EventErrors",&EventError);
	vetoTree->Branch("ErrorCount",&ErrorCount);

	// ================= 1st loop over veto entries  ==============
	// Measure the LED frequency, find the first good entry,
	// highest multiplicity, and SBC offset.
	TH1D *LEDDeltaT = new TH1D("LEDDeltaT","LEDDeltaT",100000,0,100); // 0.001 sec/bin

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (!veto.GetBadScaler()) {
			BadScalers.push_back(0);
			xTime = veto.GetTimeSec();
		}
		else {
			BadScalers.push_back(1);
			xTime = ((double)i / vEntries) * duration; // this breaks if we have corrupted duration
		}
		EntryNum.push_back(i);
		EntryTime.push_back(xTime);

		if (foundFirst && veto.GetError(1)) {
			foundFirst = false;
		}
		if (!foundFirstScaler && !veto.GetError(4)){
			foundFirstScaler = true;
			firstGoodScaler = veto.GetTimeSec();
		}
		if (CheckEventErrors(veto,prev,first,prevGoodEntry)){
			skippedEvents++;
			// do end of loop resets before continuing
			prev = veto;
			prevGoodTime = xTime;
			prevGoodEntry = i;
			veto.Clear();
			continue;
		}
		if (!foundFirst && veto.GetTimeSBC() > 0 && veto.GetTimeSec() > 0 && !veto.GetError(4)) {
			first = veto;
			foundFirst = true;
		}
		if (veto.GetMultip() > highestMultip) {
			highestMultip = veto.GetMultip();
		}
		// cout << veto.GetMultip() << endl;
		if (veto.GetMultip() > LEDSimpleThreshold) {
			LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
			simpleLEDCount++;
		}
		// end of loop
		prev = veto;
		prevGoodTime = xTime;
		prevGoodEntry = i;
		veto.Clear();
	}

	// =====================  Run-level checks =====================

	SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
	if (duration <= 0) {
		cout << "Corrupted duration.  Did we get a stop packet?\n";
		cout << "   Raw duration is " << duration << "  start: " << start << " stop: " << stop << endl;
		cout << "   Last good timestamp: " << prevGoodTime-firstGoodScaler << endl;
		duration = prevGoodTime-firstGoodScaler;
		cout << "   Set duration to " << duration << endl;
	}
	livetime = duration - (first.GetTimeSec() - firstGoodScaler);
	cout << "Veto livetime: " << livetime << " seconds\n";

	// set LED multiplicity threshold
	multipThreshold = highestMultip - LEDMultipThreshold;
	if (multipThreshold < 0) multipThreshold = 0;

	// find LED frequency, and try to adjust if we have a short run (only a few LED events)
	int dtEntries = LEDDeltaT->GetEntries();
	if (dtEntries > 0) {
		int maxbin = LEDDeltaT->GetMaximumBin();
		LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
		LEDrms = LEDDeltaT->GetRMS();
		LEDfreq = 1/LEDDeltaT->GetMean();
	}
	else {
		cout << "Warning! No multiplicity > " << LEDSimpleThreshold << " events.  LED may be off.\n";
		LEDrms = 9999;
		LEDfreq = 9999;
		badLEDFreq = true;
	}
	LEDperiod = 1/LEDfreq;
	delete LEDDeltaT;
	if (LEDperiod > 9 || vEntries < 100)
	{
		cout << "Warning: Short run.\n";
		if (simpleLEDCount > 3) {
			cout << "   From histo method, LED freq is " << LEDfreq
				 << "  Using approximate rate: " << simpleLEDCount/duration << endl;
			LEDperiod = duration/simpleLEDCount;
		}
		else {
			LEDperiod = 9999;
			badLEDFreq = true;
		}
	}
	if (LEDperiod > 20 || LEDperiod < 0 || badLEDFreq) {
		ErrorCount[25]++;
		EventError[25] = true;
	}

	// ========= 2nd loop over entries - Error checks =========

	reader.SetTree(vetoChain); // reset the reader
	std::fill(EventError.begin(), EventError.end(), 0); // reset error bools
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		// this time we don't skip anything, and count up each type of error.
		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		CheckEventErrors(veto,prev,first,prevGoodEntry,EventError);
		for (int j=0; j<nErrs; j++) if (EventError[j]==1) ErrorCount[j]++;

		// find event time
		if (!veto.GetBadScaler())
		{
			xTime = veto.GetTimeSec();
			if(run > 8557 && veto.GetTimeSBC() < 2000000000)
				timeSBC = veto.GetTimeSBC() - SBCOffset;
		}
		else if (run > 8557 && veto.GetTimeSBC() < 2000000000)
			xTime = veto.GetTimeSBC() - SBCOffset;
		else
			xTime = InterpTime(i,EntryTime,EntryNum,BadScalers);
		EntryTime[i] = xTime;	// replace entry with the more accurate one

		// Print errors to screen
		bool PrintError = false;
		for (auto i : SeriousErrors){
			if (EventError[i]) {
				PrintError = true;
				break;
			}
		}
		if (PrintError)
		{
			cout << "\nSerious errors found in entry " << i << ":\n";

			// debug block (don't delete!) used to compare with original vetoCheck code
			// cout << veto.GetTimeSec() << "\n" 	// STime
			// 	  << prev.GetTimeSec() << "\n" 	// STimePrev
			// 	  << timeSBC << "\n"					// SBCTime
			// 	  << SBCOffset << "\n"				// SBCOffset
			// 	  << prevtimeSBC << "\n"			// SBCTimePrev
			// 	  << veto.GetError(1) << "\n"
			// 	  << i  << "\n"						// entry
			// 	  << first.GetEntry() << "\n"		// firstGoodEntry
			// 	  << "error 18: " << EventError[18] << "\n"			// Error[18]
			// 	  << veto.GetSEC() << "\n"
			// 	  << prev.GetSEC() << "\n"
			// 	  << prevGoodEntry << "\n"			//  EventNumPrev_good
			// 	  << "error 20: " << EventError[20] << "\n";

			if (EventError[1]) {
				cout << "EventError[1] Missing Packet."
					 << "  Scaler index " << veto.GetScalerIndex()
					 << "  Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << veto.GetTimeSBC() << "\n";
			}
			if (EventError[13]) {
				cout << "EventError[13] ORCA packet indexes of QDC1 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC1 Index " << veto.GetQDC1Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC1 Index " << prev.GetQDC1Index() << endl;
			}
			if (EventError[14]) {
				cout << "EventError[14] ORCA packet indexes of QDC2 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC2 Index " << veto.GetQDC2Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC2 Index " << prev.GetQDC2Index() << endl;
			}
			if (EventError[18]) {
				cout << "EventError[18] Scaler/SBC Desynch."
					 << "\n    DeltaT (adjusted) " << veto.GetTimeSec() - timeSBC - TSdifference
					 << "  DeltaT " << veto.GetTimeSec() - timeSBC
					 << "\n    Prev TSdifference " << TSdifference
					 << "  Scaler DeltaT " << veto.GetTimeSec()-prev.GetTimeSec()
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  Previous Scaler Index " << prev.GetScalerIndex()
					 << "\n    Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << timeSBC << "\n";
				TSdifference = veto.GetTimeSec() - timeSBC;
			}
			if (EventError[19]) {
				cout << "EventError[19] Scaler Event Count Reset. "
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
			}
			if (EventError[20]) {
				cout << "EventError[20] Scaler Event Count Jump."
					 << "    xTime " << xTime
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "\n    SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
			}
			if (EventError[21]) {
				cout << "EventError[21] QDC1 Event Count Reset."
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  QEC1 " << veto.GetQEC()
					 << "  Previous QEC1 " << prev.GetQEC() << "\n";
			}
			if(EventError[22]) {
				cout << "EventError[22] QDC 1 Event Count Jump."
					 << "  xTime " << xTime
					 << "  QDC 1 Index " << veto.GetQDC1Index()
					 << "  QEC 1 " << veto.GetQEC()
					 << "  Previous QEC 1 " << prev.GetQEC() << "\n";
			}
			if (EventError[23]) {
				cout << "EventError[23] QDC2 Event Count Reset."
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  QEC2 " << veto.GetQEC2()
					 << "  Previous QEC2 " << prev.GetQEC2() << "\n";
			}
			if(EventError[24]) {
				cout << "EventError[24] QDC 2 Event Count Jump."
					 << "  xTime " << xTime
					 << "  QDC 2 Index " << veto.GetQDC2Index()
					 << "  QEC 2 " << veto.GetQEC2()
					 << "  Previous QEC 2 " << prev.GetQEC2() << "\n";
			}
		}

		TSdifference = veto.GetTimeSec() - timeSBC;
		prevtimeSBC = timeSBC;
		timeSBC = 0;
		prev = veto;
		last = prev;
		prevGoodEntry = i;
		veto.Clear();

		// Reset error bools each entry
		std::fill(EventError.begin(), EventError.end(), 0);
	}

	// Calculate total errors and total serious errors
	for (int i = 1; i < nErrs; i++)
	{
			// Ignore 10 and 11, these will always be present
			// as long as the veto counters are not reset at the beginning of runs.
			if (i != 10 && i != 11) TotalErrorCount += ErrorCount[i];

			// make sure LED being off is counted as serious
			if (i == 25 && ErrorCount[i] > 0) SeriousErrorCount += ErrorCount[i];

			// count up the serious errors in the vector
			for (auto j : SeriousErrors)
				if (i == j)
					SeriousErrorCount += ErrorCount[i];
	}
	cout << "=================== Veto Error Report ===================\n";
	cout << "Serious errors found :: " << SeriousErrorCount << endl;
	if (SeriousErrorCount > 0)
	{
		cout << "Total Errors : " << TotalErrorCount << endl;
		if (duration != livetime)
			cout << "Run duration (" << duration << " sec) doesn't match live time: " << livetime << endl;

		for (int i = 1; i < nErrs; i++)
		{
			if (ErrorCount[i] > 0)
			{
				if (i != 25)
					cout << "  Error[" << i <<"]: " << ErrorCount[i] << " events ("
						 << 100*(double)ErrorCount[i]/vEntries << " %)\n";
		 		else if (i == 25)
				{
					cout << "  EventError[25]: Bad LED rate: " << LEDfreq << "  Period: " << LEDperiod << endl;
					if (LEDperiod > 0.1 && (abs(duration/LEDperiod) - simpleLEDCount) > 5)
					{
						cout << "   Simple LED count: " << simpleLEDCount
							 << "  Expected: " << (int)(duration/LEDperiod) << endl;
					}
				}
			}
		}
		cout << "For reference, \"serious\" error types are: ";
		for (auto i : SeriousErrors) cout << i << " ";
		cout << "\nPlease report these to the veto group.\n";
	}
	if (errorCheckOnly) return;

	// ========= 3rd loop over veto entries - Find muons! Write ROOT output! =========

	cout << "================= Scanning for muons ... ================\n";
	reader.SetTree(vetoChain); // reset the reader
	prev.Clear();
	skippedEvents = 0;
	printf("Highest multiplicity found: %i.  Using LED threshold: %i\n",highestMultip,multipThreshold);

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		CheckEventErrors(veto,prev,first,prevGoodEntry,EventError);
		for (int j=0; j<nErrs; j++) if (EventError[j]==1) ErrorCount[j]++;

		// Calculate time of event.
		// TODO: implement an estimate of the error when alternate methods are used.
		ApproxTime = false;
		if (!veto.GetBadScaler())
		{
			xTime = veto.GetTimeSec();
			if(run > 8557 && veto.GetTimeSBC() < 2000000000)
				timeSBC = veto.GetTimeSBC() - SBCOffset;
		}
		else if (run > 8557 && veto.GetTimeSBC() < 2000000000) {
			timeSBC = veto.GetTimeSBC() - SBCOffset;
			xTime = timeSBC;
		}
		else {
			xTime = InterpTime(i,EntryTime,EntryNum,BadScalers);
			ApproxTime = true;
		}
		// Scaler jump handling:
		// When the error is initally found, we adjust xTime by "DeltaT (adjusted)" from EventError[18].
		// If we find a scaler jump in one event, we force the scaler and SBC times to match
		// for the rest of the run, or until they come back in sync on their own.
		deltaTadj = veto.GetTimeSec() - timeSBC - TSdifference;
		if (EventError[18]) {
			foundScalerJump=true;
			xTime = xTime - deltaTadj;
		}
		if (foundScalerJump) {
			ApproxTime = true;
			xTime = xTime - TSdifference;
		}
		if (fabs(TSdifference) < 0.001)	// SBC is accurate to microseconds
			foundScalerJump = false;
		// debug block for a run w/ jump at i==248 (don't delete!)
		// if (i==247||i==248||i==249||i==250)
			// printf("%i  e18: %i  %-6.3f  %-6.3f  %-6.3f  %-6.3f  %-6.3f\n", 		i,EventError[18],veto.GetTimeSec(),timeSBC,xTime,TSdifference,deltaTadj);

		// Skip events after the event time is calculated,
		// so we still properly reset for the next event.

		// These events are tagged as "badEvent==1" in the ROOT output.
		badEvent = CheckEventErrors(veto,prev,first,prevGoodEntry,EventError);
		if (badEvent)
		{
			skippedEvents++;
			// do the end-of-event reset
			vetoTree->Fill();
			TSdifference = veto.GetTimeSec() - timeSBC;
			prevtimeSBC = timeSBC;
			timeSBC = 0;
			prev = veto;
			last = prev;
			prevGoodEntry = i;
			veto.Clear();
			continue;
		}

		// LED / Time Cut

		// IsLED = true;
		TimeCut = false;
		LEDTurnedOff = ErrorCount[25];
		x_deltaT = xTime - xTimePrevLED;

		// TODO: right now this enforces a hard multiplicity cut.
		// Need to evaluate if the frequency-based multiplicity cut will
		// recover high-multiplicity muon events without mis-identifying LED's as muons.
		if (veto.GetMultip() < multipThreshold && !LEDTurnedOff) {
			// IsLED = false;
			TimeCut = true;
		}
		else if (LEDTurnedOff)
			TimeCut = true;

		/*
		// if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) < LEDWindow && veto.GetMultip() > multipThreshold)
		// {
		// 	TimeCut = false;
		// 	IsLED = true;
		// }
		// // almost missed a high-multiplicity event somehow ...
		// // often due to skipping previous events.
		// else if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) >= (LEDperiod - LEDWindow)
		// 				&& veto.GetMultip() > multipThreshold && i > 1)
		// {
		// 	TimeCut = false;
		// 	IsLED = true;
		// 	cout << "Almost missed LED:\n";
		// 	printf("Current: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
		// 		,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);
		// }
		// else TimeCut = true;
		// // Grab first LED
		// if (!LEDTurnedOff && !firstLED && veto.GetMultip() > multipThreshold) {
		// 	printf("Found first LED.  i %-2li m %-2i t %-5.2f thresh:%i  LEDoff:%i\n\n",i,veto.GetMultip(),xTime,multipThreshold,LEDTurnedOff);
		// 	IsLED=true;
		// 	firstLED=true;
		// 	TimeCut=false;
		// 	x_deltaT = -1;
		// }
		// // If frequency measurement is bad, revert to standard multiplicity cut
		// if (badLEDFreq && veto.GetMultip() >= LEDSimpleThreshold){
		// 	IsLED = true;
		// 	TimeCut = false;
		// }
		// // Simple x_LEDDeltaT uses the multiplicity-only threshold, veto.GetMultip() > multipThreshold.
		// x_LEDDeltaT = xTime - xTimePrevLEDSimple;
		//
		// // If LED is off, all events pass time cut.
		// if (LEDTurnedOff) {
		// 	IsLED = false;
		// 	TimeCut = true;
		// }
		// // Check output
		// printf("%-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
		// 	,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);
		*/

    	// Energy (Gamma) Cut
    	// The measured muon energy threshold is QDC = 500.
    	// Set TRUE if at least TWO panels are over 500.
		EnergyCut = false;
    	int over500Count = 0;
    	for (int q = 0; q < 32; q++) {
    		if (veto.GetQDC(q) > 500)
    			over500Count++;
    	}
    	if (over500Count >= 2) EnergyCut = true;

		// debug block (don't delete): check entries
		// printf("Entry %li  Time %-6.2f  QDC %-5i  Mult %i  hits>500 %i  Loff? %i  T/LCut %i  ECut %i  \n",i,xTime,veto.GetTotE(),veto.GetMultip(),over500Count,LEDTurnedOff,TimeCut,EnergyCut);

		// Hit Pattern "Cut": Map hits above SW threshold to planes and count the hits.
		PlaneHitCount = 0;
		for (int k = 0; k < 12; k++) {
			PlaneTrue[k] = 0;
			PlaneHits[k]=0;
		}
		for (int k = 0; k < 32; k++) {
			if (veto.GetQDC(k) > veto.GetSWThresh(k)) {
				if (PanelMap(k)==0) { PlaneTrue[0]=1; PlaneHits[0]++; }			// 0: Lower Bottom
				else if (PanelMap(k)==1) { PlaneTrue[1]=1; PlaneHits[1]++; }		// 1: Upper Bottom
				else if (PanelMap(k)==2) { PlaneTrue[2]=1; PlaneHits[2]++; }		// 3: Inner Top
				else if (PanelMap(k)==3) { PlaneTrue[3]=1; PlaneHits[3]++; }		// 4: Outer Top
				else if (PanelMap(k)==4) { PlaneTrue[4]=1; PlaneHits[4]++; }		// 5: Inner North
				else if (PanelMap(k)==5) { PlaneTrue[5]=1; PlaneHits[5]++; }		// 6: Outer North
				else if (PanelMap(k)==6) { PlaneTrue[6]=1; PlaneHits[6]++; }		// 7: Inner South
				else if (PanelMap(k)==7) { PlaneTrue[7]=1; PlaneHits[7]++; }		// 8: Outer South
				else if (PanelMap(k)==8) { PlaneTrue[8]=1; PlaneHits[8]++; }		// 9: Inner West
				else if (PanelMap(k)==9) { PlaneTrue[9]=1; PlaneHits[9]++; }		// 10: Outer West
				else if (PanelMap(k)==10) { PlaneTrue[10]=1; PlaneHits[10]++; }	// 11: Inner East
				else if (PanelMap(k)==11) { PlaneTrue[11]=1; PlaneHits[11]++; }	// 12: Outer East
			}
		}
		for (int k = 0; k < 12; k++) if (PlaneTrue[k]) PlaneHitCount++;

		// Muon Identification
		// Use EnergyCut, TimeCut, and the Hit Pattern to identify them sumbitches.
		for (int r = 0; r < 32; r++) {CoinType[r]=0; CutType[r]=0;}
		if (TimeCut && EnergyCut)
		{
			int type = 0;
			bool a=0,b=0,c=0;
			CoinType[0] = true;

			if (PlaneTrue[0] && PlaneTrue[1] && PlaneTrue[2] && PlaneTrue[3]) {
				CoinType[1]=true;
				a=true;
				type=1;
			}
			if ((PlaneTrue[0] && PlaneTrue[1]) && ((PlaneTrue[2] && PlaneTrue[3]) || (PlaneTrue[4] && PlaneTrue[5])
					|| (PlaneTrue[6] && PlaneTrue[7]) || (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))){
				CoinType[2] = true;
				b=true;
				type = 2;
			}
			if ((PlaneTrue[2] && PlaneTrue[3]) && ((PlaneTrue[4] && PlaneTrue[5]) || (PlaneTrue[6] && PlaneTrue[7])
					|| (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))) {
				CoinType[3] = true;
				c=true;
				type = 3;
			}
			if ((a && b)||(a && c)||(b && c)) type = 4;

			char hitType[200];
			if (type==0) sprintf(hitType,"2+ panels");
			if (type==1) sprintf(hitType,"vertical");
			if (type==2) sprintf(hitType,"side+bottom");
			if (type==3) sprintf(hitType,"top+sides");
			if (type==4) sprintf(hitType,"compound");

			// print the details of the hit
			printf("Hit: %-12s Entry %-4li Time %-6.2f  QDC %-5i  Mult %i  LEDoff %i  ApxT %i\n", hitType,i,xTime,veto.GetTotE(),veto.GetMultip(),LEDTurnedOff,ApproxTime);
		}

		// Output

		CutType[0] = LEDTurnedOff;
		CutType[1] = EnergyCut;
		CutType[2] = ApproxTime;
		CutType[3] = TimeCut;
		CutType[4] = IsLED;
		CutType[5] = firstLED;
		CutType[6] = badLEDFreq;

		out = veto;
		vetoTree->Fill();

		// Reset for next entry

		// resets used by error check (don't delete!)
		TSdifference = veto.GetTimeSec() - timeSBC;
		prevtimeSBC = timeSBC;
		timeSBC = 0;
		prev = veto;
		last = prev;
		prevGoodEntry = i;

		// resets used by muon finder (may need revision)
		if (IsLED) {
			prevLED = veto;
			xTimePrevLED = xTime;
		}
		if (veto.GetMultip() > multipThreshold) {
			xTimePrevLEDSimple = xTime;
		}
		IsLEDPrev = IsLED;
		xTimePrev = xTime;
		x_deltaTPrev = x_deltaT;

		veto.Clear();
	}
	if (skippedEvents > 0) printf("ProcessVetoData skipped %li of %li entries.\n",skippedEvents,vEntries);

	vetoTree->Write("",TObject::kOverwrite);
	RootFile->Close();
	cout << "Wrote ROOT file: " << outputFile << endl;
}
// ====================================================================================
// ====================================================================================

void SetCardNumbers(int runNum, int &card1, int &card2)
{
	if (runNum > 45000000){
		card1 = 11;
		card2 = 18;
	}
	else {
		card1 = 13;
		card2 = 18;
	}
}

int FindPanelThreshold(TH1D *qdcHist, int threshVal, int panel, int runNum)
{
	if (runNum > 45000000 && panel > 23)
		return 9999;

	int firstNonzeroBin = qdcHist->FindFirstBinAbove(1,1);
	qdcHist->GetXaxis()->SetRange(firstNonzeroBin-10,firstNonzeroBin+50);
	//qdcHist->GetXaxis()->SetRangeUser(0,500); //alternate method of finding pedestal
	int bin = qdcHist->GetMaximumBin();
	if (firstNonzeroBin == -1) return -1;
	double xval = qdcHist->GetXaxis()->GetBinCenter(bin);
	return xval+threshVal;
}

double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler)
{
	if ((times.size() != entries.size()) || times.size() != badScaler.size())
	{
		cout << "Vectors are different sizes!\n";
		if (entry >= (int)times.size())
			cout << "Entry is larger than number of entries in vector!\n";
		return -1;
	}

	double iTime = 0;
	double lower = 0;
	double upper = 0;
	if (!badScaler[entry]) iTime = times[entry];
	else
	{
		for (int i = entry; i < (int)entries.size(); i++)
			if (badScaler[i] == 0) { upper = times[i]; break; }
		for (int j = entry; j > 0; j--)
			if (badScaler[j] == 0) { lower = times[j]; break; }
		iTime = (upper + lower)/2.0;
	}
	return iTime;
}

int PanelMap(int i, int runNum)
{
	// TODO: update this for different run numbers / configurations!

	// For tagging plane-based coincidences.
	// This uses a zero-indexed map.

	// 0: Lower Bottom
	// 1: Upper Bottom
	// 2: Inner Top
	// 3: Outer Top
	// 4: Inner North
	// 5: Outer North
	// 6: Inner South
	// 7: Outer South
	// 8: Inner West
	// 9: Outer West
	// 10: Inner East
	// 11: Outer East

	if 		(i == 0) return 0;  // L-bot 1
	else if (i == 1) return 0;
	else if (i == 2) return 0;
	else if (i == 3) return 0;
	else if (i == 4) return 0;
	else if (i == 5) return 0;  // L-bot 6

	else if (i == 6) return 1;  // U-bot 1
	else if (i == 7) return 1;
	else if (i == 8) return 1;
	else if (i == 9) return 1;
	else if (i == 10) return 1;
	else if (i == 11) return 1; // U-bot 6

	else if (i == 17) return 3; // Top outer
	else if (i == 18) return 3; // Top outer
	else if (i == 20) return 2; // Top inner
	else if (i == 21) return 2; // Top inner

	else if (i == 15) return 5; // North outer
	else if (i == 16) return 5; // North outer
	else if (i == 19) return 4; // North inner
	else if (i == 23) return 4; // North inner

	else if (i == 24) return 6; // South inner
	else if (i == 25) return 7; // South outer
	else if (i == 26) return 6; // South inner
	else if (i == 27) return 7; // South outer

	else if (i == 12) return 8; // West inner
	else if (i == 13) return 8; // West inner
	else if (i == 14) return 9; // West outer
	else if (i == 22) return 9; // West outer

	else if (i == 28) return 10; // East inner
	else if (i == 29) return 11; // East outer
	else if (i == 30) return 10; // East inner
	else if (i == 31) return 11; // East outer

	else return -1;
}

// This is overloaded so we don't have to use the error vector if we don't need it
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry)
{
	vector<int> ErrorVec(29);
	std::fill(ErrorVec.begin(), ErrorVec.end(), 0);
	return CheckEventErrors(veto,prev,first,prevGoodEntry,ErrorVec);
}
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry, vector<int> &ErrorVec)
{
	// Returns skip = false if the event is analyzable (either clean, or a workaround exists)
	bool skip = false;

	/*
	Recoverable:
	Actionable:
		LED Off -> Output string Dave can grep for, send veto experts an email
		No QDC events -> Need to have people check if a panel went down
		High error rate (desynchs, bad scalers, etc) -> Need to try and restart data taking
	Diagnostic:

	// Event-level error checks ('s' denotes setting skip=true)
	// s 1. Missing channels (< 32 veto datas in event)
	// s 2. Extra Channels (> 32 veto datas in event)
	// s 3. Scaler only (no QDC data)
	//   4. Bad Timestamp: FFFF FFFF FFFF FFFF
	// s 5. QDCIndex - ScalerIndex != 1 or 2
	// s 6. Duplicate channels (channel shows up multiple times)
	//   7. HW Count Mismatch (SEC - QEC != 1 or 2)
	//   8. MJTRun run number doesn't match input file
	// s 9. MJTVetoData cast failed (missing QDC data)
	//   10. Scaler EventCount doesn't match ROOT entry
	//   11. Scaler EventCount doesn't match QDC1 EventCount
	//   12. QDC1 EventCount doesn't match QDC2 EventCount
	// s 13. Indexes of QDC1 and Scaler differ by more than 2
	// s 14. Indexes of QDC2 and Scaler differ by more than 2
	//   15. Indexes of either QDC1 or QDC2 PRECEDE the scaler index
	//   16. Indexes of either QDC1 or QDC2 EQUAL the scaler index
	//   17. Unknown Card is present.
	// s 18. Scaler & SBC Timestamp Desynch.
	// s 19. Scaler Event Count reset.
	// s 20. Scaler Event Count increment by > +1.
	// s 21. QDC1 Event Count reset.
	// s 22. QDC1 Event Count increment by > +1.
	// s 23. QDC2 Event Count reset.
	// s 24. QDC2 Event Count increment > +1.
	//   25. Used interpolated time

	// Run-level error checks (not checked in this function)
	// 26. LED frequency very low/high, corrupted, or LED's off.
	// 27. QDC threshold not found
	// 28. No events above QDC threshold
	*/


	vector<int> EventErrors(29);
	std::fill(EventErrors.begin(), EventErrors.end(), 0);

	// Errors 1-18 are checked automatically when we call MJVetoEvent::WriteEvent
	for (int i=0; i < 18; i++) EventErrors[i] = veto.GetError(i);

	for (int q=0; q<18; q++) {
		if (EventErrors[q]==1 && (q==1||q==2||q==3||q==5||q==6||q==9||q==13||q==14))
			skip = true;
	}

	if (veto.GetBadScaler() && (veto.GetRun() < 8557 || veto.GetTimeSBC() > 2000000000))
		EventErrors[25] = true;

	int entry = veto.GetEntry();
	int firstGoodEntry = first.GetEntry();

	// assign EventErrors to the output vetor
	ErrorVec = EventErrors;

	// debug: print error vector
	// for (int i=0; i < 26; i++) cout << i << ":" << EventErrors[i] << "  ";
	// cout << endl;

	// return if we haven't found the first good entry yet
	if (firstGoodEntry == -1) return skip;

	double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
	double timeSBC = veto.GetTimeSBC() - SBCOffset;
	double prevtimeSBC = prev.GetTimeSBC() - SBCOffset;

	if (veto.GetTimeSec() > 0 && timeSBC > 0 && SBCOffset !=0
			&& !veto.GetError(1) && entry > firstGoodEntry
			&& fabs((veto.GetTimeSec() - prev.GetTimeSec())-(timeSBC - prevtimeSBC)) > 2)
				EventErrors[18] = true;

	if (veto.GetSEC() == 0 && entry != 0 && entry > firstGoodEntry)
		EventErrors[19] = true;

	if (abs(veto.GetSEC() - prev.GetSEC()) > entry-prevGoodEntry
			&& entry > firstGoodEntry && veto.GetSEC()!=0)
		EventErrors[20] = true;

	if (veto.GetQEC() == 0 && entry != 0 && entry > firstGoodEntry && !veto.GetError(1))
		EventErrors[21] = true;

	if (abs(veto.GetQEC() - prev.GetQEC()) > entry-prevGoodEntry
			&& entry > firstGoodEntry && veto.GetQEC() != 0)
		EventErrors[22] = true;

	if (veto.GetQEC2() == 0 && entry != 0 && entry > firstGoodEntry && !veto.GetError(1))
		EventErrors[23] = true;

	if (abs(veto.GetQEC2() - prev.GetQEC2()) > entry-prevGoodEntry
			&& entry > firstGoodEntry && veto.GetQEC2() != 0)
		EventErrors[24] = true;

	// Check errors 18-25
	for (int q=18; q<26; q++) {
		if (EventErrors[q]==1 && (q==18||q==19||q==20||q==21||q==22||q==23||q==24))
			skip = true;
	}
	// assign EventErrors to the output vetor
	ErrorVec = EventErrors;

	return skip;
}
