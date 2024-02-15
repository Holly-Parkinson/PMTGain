//This macro was made by abullock and aszelc for SBND, July 2023
//For a given channel, scan waveforms for single photoelectron (SPE) events and perform an analysis on them
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include <iostream>
#include "TStopwatch.h"
#include "TMath.h"

using namespace std;


void speana_noTSpec_single() {
  TStopwatch timer;
  timer.Start();


  //CONFIG
  Int_t opc = 85; //opchannel as an integer
  Int_t eventid = 1; //event id as an integer
  bool all_events = true; //do all events or just the selected one specified by eventid: RIGHT NOW, JUST ONE CHANNEL
  bool cut = false; //make the 100ns cut (disregard-we don't use this anymore)
  bool save = true; // do we want to save the analysis histograms? (probably)
  Option_t *filemode = "RECREATE"; //how to initialize the file, see TFile() reference
  bool do_avgspe=true, do_amp=true, do_integ=true; // which analysis do we want to do?
  bool draw_avgspe=false, draw_amp=false, draw_integ=false; // which analysis do we want to draw?
  int lowbin=200,hibin=200; //sample range around the peak
  Int_t spacethresh = 100; //average spacing in ns above which a microsecond range qualifies for SPEs
  int nstdev=3; //number of noise stdevs to set the threshold to (I used nstdev=2, you might get better results with 3)
  int spe_region_start = 1500;
  double* gamma; //was advised to keep outside of main loop
  double nbmax_factor=0.9,nbmin_factor=0.1;  //set noise analysis range 1 (pre peaks): noisebinmax = highestbin*nbmax_factor
  double n2bmin_factor=0.9; //set noise analysis range 2 (post peaks), upper bound is just wvf_nbins: noisebin2min = n2bmin_factor*wvf_nbins

  double pulse_peak;
  double pulse_t_start;
  double pulse_t_end;
  double pulse_t_peak;
  vector<double> pulse_t_peak_v;

  //VARIABLE CREATION
  int NBINS=hibin+lowbin+1; //number of total bins around the peak
  int navspes=0; //number of SPEs in the analysis histograms
  int success=0; //number of successful waveform analyses
  int failed=0; //number of failed waveform analyses


  //TSpectrum *s = new TSpectrum(200); //create new TSpectrum object with 200 peaks (should be plenty) TAKEN OUT
  TH1D *avgspe = new TH1D(Form("avgspe_opchannel_%d", opc), "Average SPE Shape;Samples from peak;Count", NBINS, -lowbin, hibin); // create histogram for average shape
  TH1D *amp = new TH1D(Form("amp_opchannel_%d", opc), "Amplitudes of SPEs;Amplitude [ADC];Count", 50, 0, 200); // create histogram for amplitude
  TH1D *integ0 = new TH1D(Form("integ_opchannel_%d_zeromode", opc), "Integral of SPEs;Integral value [ADC*samples];Count", 50, 0, 500); // create histogram for integral
  TH1D *integ1 = new TH1D(Form("integ_opchannel_%d_threshmode", opc), "Integral of SPEs;Integral value [ADC*samples];Count", 50, 0, 500); // create histogram for integral
  TH1D *integ2 = new TH1D(Form("integ_opchannel_%d_manualmode", opc), "Integral of SPEs;Integral value [ADC*samples];Count", 50, 0, 500); // create histogram for integral
  TH1D *integ3 = new TH1D(Form("integ_opchannel_%d_zeromodeB", opc), "Integral of SPEs;Integral value [ADC*samples];Count", 50, 0, 500); // create histogram for integral
  TH1D *integ4 = new TH1D(Form("integ_opchannel_%d_threshmodeB", opc), "Integral of SPEs;Integral value [ADC*samples];Count", 50, 0, 500); // create histogram for integral
  TH1D *integ5 = new TH1D(Form("integ_opchannel_%d_manualmodeB", opc), "Integral of SPEs;Integral value [ADC*samples];Count", 50, 0, 500); // create histogram for integral

  //INITIAL PRINTOUT
  cout << "======" << "SPE ANALYSIS" << "======" << endl << "Developed by abullock for SBND, July 2023." << endl << "Channel selected: " << opc << endl;
  if (!all_events) {cout << "Event selected: " << eventid << endl;} else {cout << "All events selected." << endl;}
  cout << "Launching..." << endl;

  //OPEN THE FILE, GET KEYS FOR RECO WAVEFORMS
  TFile *myWvf = new TFile("/exp/sbnd/data/users/hollyp/PMTGain/wvf_ana.root"); // open the waveform file and store under pointer myWvf
  TFile *wvfana = (TFile*) myWvf->Get("wvfana"); // get subdirectory with waveforms in it
  TList *keys = wvfana->GetListOfKeys(); // get list of all the waveforms in wvfana
  int nkeys = keys->GetEntries();

  //OPEN THE FILE, GET TTREE FOR TRUTH
  TFile *myTruth = new TFile("/exp/sbnd/data/users/hollyp/PMTGain/opdet_ana.root"); // open the file and store under pointer myTruth
  TFile *opdetana = (TFile*) myTruth->Get("opdetana"); // get subdirectory
  TTree *a = (TTree*) opdetana->Get("AllPhotons"); // get the ttree out of the file  

  //SETUP TTREE FOR READING
  Long64_t ne = a->GetEntries(); // get number of entries in ttree
  Float_t time; a->SetBranchAddress("Time", &time); // select time branch
  Int_t channel; a->SetBranchAddress("OpChannel", &channel); // channel branch
  Int_t eventID; a->SetBranchAddress("EventID", &eventID); //event branch

  //SELECT TRUTH TIMES TO INCLUDE (only for all_events=true, this is redone later if all_events=false)
  vector <Double_t> gammaV (0); // the truth values of photon times
  for (Long64_t i=0; i<ne; i++) { // fill histogram
    a->GetEntry(i);
    if (channel==opc) { // only add entries from the right channel and event
      gammaV.push_back(time);
    }
  }
  long long gammasize = gammaV.size(); //get the size of the vector
  gamma = &gammaV[0]; //copy the vector into an array

//LOOP THROUGH KEYS
bool analysis_active=false; //is this key being analyzed?
for (int q=0; q<nkeys; q++) { //loop through keys
  TKey *key = (TKey*) keys->At(q); //get a key
  TH1D *wvf = (TH1D*) key->ReadObj(); //get the waveform from the key
  TString *wvf_name = new TString(wvf->GetName()); //get the name of the waveform  

  //SELECT ONLY THE RIGHT OPTICAL CHANNEL (aszelc wrote this and I don't know how it works -abullock)
  Int_t start_index=wvf_name->Index("opchannel")+10; 
  TString subs=(*wvf_name)(start_index,4).Data();
  Int_t length=subs.First('_');
  Int_t start_index2=wvf_name->Index("event")+6;
  TString subs2=(*wvf_name)(start_index,4).Data();
  Int_t length2=subs2.First('_');
  bool right_channel = false; bool right_event = false;
  if (atoi((*wvf_name)(start_index,length).Data())==opc) {right_channel=true;}
  if (atoi((*wvf_name)(start_index2,length2).Data())==eventid) {right_event=true;}
  if (all_events) {right_event=true;}
if (right_channel && right_event) { //only continue if it's the right opc
  analysis_active=true; //this key is being analyzed
  cout << "Analyzing " << *wvf_name << endl;
  //Int_t new_index=wvf_name->Index("event_5_opchannel_85_pmt_coated_1380");
  //std::cout << new_index << std::endl;
  //std::cout << "test 1" << std::endl;
  //if (new_index == -1) continue;

  //GET INFO ON THE WAVEFORM
  Int_t wvf_nbins = wvf->GetNbinsX(); // get number of bins
  Double_t wvf_tlo = wvf->GetXaxis()->GetBinLowEdge(1); // get lower bound of t on histogram
  Double_t wvf_thi = wvf->GetXaxis()->GetBinUpEdge(wvf_nbins); // get upper bound of t on histogram
  

//SELECT TRUTH TIMES TO INCLUDE (only for all_events=false)
if (!all_events) {
  gammaV.clear(); // the truth values of photon times
  for (Long64_t i=0; i<ne; i++) { // fill histogram
    a->GetEntry(i);
    if (channel==opc && eventID==eventid) { // only add entries from the right channel and event
      gammaV.push_back(time);
    }
  }
  gammasize = gammaV.size(); //get the size of the vector
  double *gamma = &gammaV[0]; //copy the vector into an array
}



  //SORT TRUTH TIMES
  Long64_t index[gammasize]; //array for sorted indeces, use gamma[index[i]] to get values in time order!
  TMath::Sort(gammasize, gamma, index, false); //sort them

  //FIND SPE SPACING PER MICROSECOND FROM TRUTH
  Int_t tcropmin = wvf_tlo*1000; // lowest time of wvf's range in ns
  Double_t spe_spaces[10]; Double_t avg; //array of average spacing and a variable to hold an average
  Double_t lasttime; //last time value recorded
  bool active, just_one; //is the scan active, and have we only found one photon
  for (int us=0; us<10; us++) { //look at each microsecond of waveform
    Double_t dtotal=0; Double_t count=0; //sum of distances, and total number of distances counted  
    active = false; just_one=false;
    for (int i=0; i<gammasize; i++) { // look through all truth times
      Double_t timevalue = gamma[index[i]];
      if ( timevalue >= (tcropmin+(us*1000)) && timevalue < (tcropmin+((us+1)*1000)) && active==false ){
        lasttime = timevalue; //the first time is recorded
        active = true; //scan now active for this microsecond
        just_one = true; //we've found a photon
      }
      else if ( timevalue >= (tcropmin+(us*1000)) && timevalue < (tcropmin+((us+1)*1000)) && active==true ){
        dtotal += (timevalue - lasttime); //record the distance between this value and the last one
        count++; //we measured a distance
        lasttime = timevalue; //set this value as the last time
        just_one = false; //there is no longer just one
      }
    }
    //once a whole microsecond is done
    if (count==0 && just_one) {avg=1000;} //set to 1000 if only one is found
    else if (count==0) {avg=0;} //set to zero if none found
    else {avg=dtotal/count;} //find average distance if found
    spe_spaces[us] =  avg; // add to array
    //cout << us << " to " << us+1 << " µs: avg is " << avg << endl;
  }



  //DETERMINE THE SPE RANGE FROM TRUTH 
  Double_t sperange_tlo=0; Double_t sperange_thi=0; //whether or not we have entered the range, and bounds of range
  for (int i=0; i<10; i++) { //look through bins of spacing histogram
    if (spe_spaces[i] >= spacethresh) { //find the first bin that exceeds the threshold
      sperange_tlo = i;
      break;
    }
  }
  for (int i=9; i>=0; i--) { //look through bins of spacing histogram again but backwards
    if (spe_spaces[i] >= spacethresh) { //find the last bin that exceeds the threshold
      sperange_thi = i+1;
      break;
    }
  }
  //cout << "The SPE range is from " << sperange_tlo << " µs to " << sperange_thi << " µs." << endl;




  //NOISE ANALYSIS
  double thresh; //tspectrum threshold will be determined by noise analysis
//  wvf->Add(wvf, -2); //flip histogram
  Int_t highestbin = wvf->GetMinimumBin(); //which bin has the highest peak?
// this needs to be set by numbers/.fcl parameters (Andrzej)
  Int_t noisebinmin = nbmin_factor*highestbin; Int_t noisebinmax = nbmax_factor*highestbin; //set noise analysis range 1 (pre peaks)
  Int_t noisebin2min = n2bmin_factor*wvf_nbins; //set noise analysis range 2 (post peaks), upper bound is just wvf_nbins
  int m = 0; //number of noise values recorded
  Double_t total = 0; Double_t totalvar = 0; //for calculation
  //total across noise region 1
  for (int i=noisebinmin; i<=noisebinmax; i++) { //iterate across noise region
   Double_t val = wvf->GetBinContent(i); //get value in bin
   total += val; //sum all values
   m++;
  }
  //total across noise region 2
  for (int i=noisebin2min; i<=wvf_nbins; i++) { //iterate across noise region
    Double_t val = wvf->GetBinContent(i); //get value in bin
    total += val; //sum all values
    m++;
  }
  //calculating mean
  Double_t mean = total / m; //this mean describes the baseline and the average value over the noise
cout << "mean = " << mean << endl;

 vector<double> wvfm(wvf_nbins, 0);
 for (int i=0; i<wvf_nbins; i++) { 
	wvfm.at(i)=mean-wvf->GetBinContent(i);
	}


  //total variance across noise region 1
  for (int i=noisebinmin; i<=noisebinmax; i++) { //iterate across noise region
   Double_t val = wvfm.at(i); //get value in bin
   totalvar += (val)*(val); //sum together the variances of each bin
  }
  //total variance across noise region 2
  for (int i=noisebin2min-1; i<wvf_nbins; i++) { //iterate across noise region
    Double_t val = wvfm.at(i); //get value in bin
    totalvar += (val)*(val); //sum together the variances of each bin
  }
  //calculating stdev
  Double_t stdev = sqrt(totalvar/m); //stdev of noise
  Double_t highestval = wvfm.at(highestbin); //height of the highest peak, from which the threshold will be determined
  //thresh = stdev * nstdev / (highestval); //threshold defined as some number of stdevs above mean, as a fraction of the highest peak height
  thresh = stdev * nstdev ; //threshold defined as some number of stdevs above mean, as a fraction of the highest peak height
  cout << "Threshold set to:   "  << thresh << endl;

// Andrzej: why would you base your threshold on the highest value? 

  //LYNN'S CODE (pasted in from line 490 at https://github.com/SBNSoftware/sbndcode/blob/5bf42ecca1bf9e8ec6477ec7b5390aa25ffab94f/sbndcode/Trigger/PMT/pmtSoftwareTriggerProducer_module.cc#L490 )

   
  bool fire = false; // bool for if pulse has been detected
  int counter = 0; // counts the bin of the waveform
  int peak_count = 0;
  
  double start_threshold = thresh; //maybe put all of these outside the loop? 
  double end_threshold = thresh; 

 

  pulse_t_peak_v.clear();
  pulse_peak = 0; pulse_t_start = 0; pulse_t_end = 0; pulse_t_peak = 0;
  counter = spe_region_start;
  for (int i=spe_region_start; i<wvfm.size(); i++){
    //for (auto const &adc : wvfm){ 
    auto const &adc = wvfm[counter];
//   std:: cout << " adc/thresh " << adc <<  " " << start_threshold << std::endl;

    if ( !fire && ((double)adc) > start_threshold ){ // if its a new pulse 
      fire = true;
      //vic: i move t_start back one, this helps with porch
      pulse_t_start = counter - 1 > 0 ? counter - 1 : counter;    
      //std::cout << counter << " " << (double)adc << std::endl;
     }
  
     else if( fire && ((double)adc) < end_threshold ){ // found end of a pulse
       fire = false;
       //vic: i move t_start forward one, this helps with tail
       pulse_t_end = counter < ((int)wvfm.size())  ? counter : counter - 1;
       //std::cout << "pulse length = " << (pulse_t_end - pulse_t_start) << std::endl;
       if (pulse_t_end - pulse_t_start > 2){
	 //std::cout << counter << " " << (double)adc << std::endl;
	 //std::cout << "pulse length = " << (pulse_t_end - pulse_t_start) << std::endl;
        cout << "peak time = " << pulse_t_peak << endl; //this is the bin number 
	 pulse_t_peak_v.push_back(pulse_t_peak);
	 peak_count++;
       }
       pulse_peak = 0; pulse_t_start = 0; pulse_t_end = 0; pulse_t_peak = 0;
 
     }   
  
     else if(fire){ // if we're in a pulse 
       // Found a new maximum:
       if ((double)adc > pulse_peak){
         pulse_peak = ((double)adc);
         pulse_t_peak = counter;

      }
       //std::cout << counter << " " << (double)adc << std::endl;
     }
    //std::cout << counter << " " << (double)adc << std::endl;
    counter++;
     if (peak_count==200) {cout << "  Analysis Failure: Threshold setting unsuccessful." << endl; failed++; continue;}
   }

  
   if(fire){ // Take care of a pulse that did not finish within the readout window.
     fire = false;
     pulse_t_end = counter - 1;
     pulse_peak = 0; pulse_t_start = 0; pulse_t_end = 0; pulse_t_peak = 0;
   }
 

  //TSPECTRUM USED TO BE HERE
  

  auto wvf_pt = pulse_t_peak_v;

 
  //ISOLATE SPE PEAKS FROM ALL THE PEAKS -- potentially edit/remove this block because our threshold window for peak searching is already the SPE range
  int nspe = 0; //how many peaks are within the spe region
  for (int i=0; i<peak_count; i++) { //loop through all peak positions
    nspe++;
  }
  Double_t wvf_spet[nspe]; int j=0; //create array for them
  for (int i=0; i<peak_count; i++) { //loop through all peak positions
  TAxis *xaxis = wvf->GetXaxis();
  // Get the x value corresponding to the bin number
  double xValue = xaxis->GetBinCenter(wvf_pt[i]);
  cout << xValue << endl;
   // if (wvf_pt[i] >= (wvf_tlo + sperange_tlo) && wvf_pt[i] <= (wvf_tlo + sperange_thi)) { COMMENTED OUT FOR TEST
      cout << "passed next stage " << wvf_pt[i] << endl; 
      wvf_spet[j] = wvf_pt[i]; //add SPEs to the array
      cout << "peak x position in SPE array " << wvf_spet[j] << endl;

      j++;
    //} COMMENTED OUT FOR TEST
  }
  if (nspe==0) {cout << "  Analysis Failure: No SPEs found in this waveform." << endl; failed++; continue;} //if no SPEs are found
  




  //AVERAGE SPE SHAPE HISTOGRAM
if (do_avgspe) {
  for (int i=0; i<nspe; i++) { //iterate through the found spe positions
      Int_t peakbin = wvf_spet[i]; // get bin associated with peak time
       cout << "peak x position " << wvf_spet[i] << endl;
   if(peakbin-lowbin <0 || peakbin+hibin>wvf_nbins) //aszelc addition—check to make sure the bin numbers do not cause an error
      {continue;}

//selection addition——prevent any peak with another peak 100ns before it from being added
    bool selected = true;
    for (int l=0; l<nspe; l++){
      Double_t separation = wvf_spet[i] - wvf_spet[l];
      if (separation<0.1 && separation>0) {selected=false;}
    }
    if (cut==false) {selected = true;} //bypasses the selection above
    if (selected) {


    for (int j=1; j<=NBINS; j++) { //loop over range surrounding peak
      Double_t le_bin = wvfm[peakbin-lowbin+j]; //add the values to the histogram
      avgspe->AddBinContent(j,le_bin);
    }
    navspes++; //added one!
    } //close if(selected)
  }
  //this will get normalized after the key loop
}

  //AMPLITUDES OF SPES
if (do_amp) {
  for (int i=0; i<navspes; i++) { //iterate through the found spe positions
    Int_t peakbin = wvf_spet[i]; //get bin associated with peak time
    Double_t peakheight =  wvfm.at(i); //get the value of the highest bin (the amplitude) WHY IS THIS GIVING NEGATIVE?? //THIS COULD BE OUR ERROR
    cout << peakheight  << " = peakheight" << endl;
    
    amp->Fill(peakheight); //add to histogram
    if (!do_avgspe) {navspes++;} //count total spes here if we aren't already
  }
}

  //INTEGRALS OF SPES
if (do_integ) {
  //without baseline subtraction
  for (int i=0; i<nspe; i++) { //iterate through the found spe positions
    Int_t peakbin = wvf_spet[i]; //get bin associated with peak time
    if(peakbin-lowbin <0 || peakbin+hibin>wvf_nbins) //aszelc addition—check to make sure the bin numbers do not cause an error
      continue;
    //zeromode bounds
    int ilo=0, ihi=0; //set bounds initially
    Double_t val = mean - (wvf->GetBinContent(peakbin)); //edited for baseline subtraction
    while (val>(0)) { //check whether a bin is below noise threshold (indicating a bound on the peak)
      ilo++;
      val = mean - (wvf->GetBinContent(peakbin-ilo)); //go one bin left
    }
    val = mean - (wvf->GetBinContent(peakbin)); //reset for other bound
    while (val>(0)) { //repeat the process for the other bound
      ihi++;
      val = mean - (wvf->GetBinContent(peakbin+ihi)); //go one bin right
    }
    ilo--; ihi--; //the bounds search stops at one greater than the true bounds
    Double_t integral = 0;
    for (int j=ilo*-1; j<=ihi; j++) { //loop over range surrounding peak
      Double_t le_bin = mean - (wvf->GetBinContent(peakbin+j)); //add the values to the histogram
      integral += le_bin;
    }
    integ0->Fill(integral); //add integral
    //threshmode bounds
    ilo=0, ihi=0; //set bounds initially
    val = mean - (wvf->GetBinContent(peakbin));
    while (val>(stdev*nstdev)) { //check whether a bin is below noise threshold (indicating a bound on the peak)
      ilo++;       
      val = mean - (wvf->GetBinContent(peakbin-ilo)); //go one bin left
    }
    val = mean - (wvf->GetBinContent(peakbin)); //reset for other bound
    while (val>(stdev*nstdev)) { //repeat the process for the other bound
      ihi++;
      val = mean - (wvf->GetBinContent(peakbin+ihi)); //go one bin right
    }
    ilo--; ihi--; //the bounds search stops at one greater than the true bounds
    integral = 0;
    for (int j=ilo*-1; j<=ihi; j++) { //loop over range surrounding peak
      Double_t le_bin = mean - (wvf->GetBinContent(peakbin+j)); //add the values to the histogram
      integral += le_bin;
    }
    integ1->Fill(integral); //add integral 
    //manualmode bounds
    ilo=7, ihi=10; //set bounds manually
    integral = 0;
    for (int j=ilo*-1; j<=ihi; j++) { //loop over range surrounding peak
      Double_t le_bin = mean - (wvf->GetBinContent(peakbin+j)); //add the values to the histogram
      integral += le_bin;
    }
    integ2->Fill(integral); //add integral
    if (!do_avgspe && !do_amp) {navspes++;} //count total spes here if we aren't already
  }

  //with baseline subtraction
  for (int i=0; i<nspe; i++) { //iterate through the found spe positions
    Int_t peakbin = wvf_spet[i]; //get bin associated with peak time
    if(peakbin-lowbin <0 || peakbin+hibin>wvf_nbins) //aszelc addition—check to make sure the bin numbers do not cause an error
      continue;
    //local baseline calculation
    Double_t vallow = mean - (wvf->GetBinContent(peakbin-50));
    Double_t valhi = mean - (wvf->GetBinContent(peakbin+50));
    Double_t bsl = (vallow + valhi) / 2; //average the two
    //zeromode bounds
    int ilo=0, ihi=0; //set bounds initially
    Double_t val = mean - (wvf->GetBinContent(peakbin) - bsl);
    while (val>(0)) { //check whether a bin is below noise threshold (indicating a bound on the peak)
      ilo++;
      val = mean - (wvf->GetBinContent(peakbin-ilo) - bsl); //go one bin left
      if (ilo==50) {break;} //keep it from going too wide
    }
    val = mean - (wvf->GetBinContent(peakbin) - bsl); //reset for other bound
    while (val>(0)) { //repeat the process for the other bound
      ihi++;
      val = mean - (wvf->GetBinContent(peakbin+ihi) - bsl); //go one bin right
      if (ihi==50) {break;} //keep it from going too wide
    }
    ilo--; ihi--; //the bounds search stops at one greater than the true bounds
    Double_t integral = 0;
    for (int j=ilo*-1; j<=ihi; j++) { //loop over range surrounding peak
      Double_t le_bin = mean - (wvf->GetBinContent(peakbin+j) - bsl); //add the values to the histogram
      integral += le_bin;
    }
    integ3->Fill(integral); //add integral
    //threshmode bounds
    ilo=0, ihi=0; //set bounds initially
    val = mean - (wvf->GetBinContent(peakbin) - bsl);
    while (val>(stdev*nstdev)) { //check whether a bin is below noise threshold (indicating a bound on the peak)
      ilo++;
      val = mean - (wvf->GetBinContent(peakbin-ilo) - bsl); //go one bin left
      if (ilo==50) {break;} //keep it from going too wide
    }
    val = mean - (wvf->GetBinContent(peakbin) - bsl); //reset for other bound
    while (val>(stdev*nstdev)) { //repeat the process for the other bound
      ihi++;
      val = mean - (wvf->GetBinContent(peakbin+ihi) - bsl); //go one bin right
      if (ihi==50) {break;} //keep it from going too wide
    }
    ilo--; ihi--; //the bounds search stops at one greater than the true bounds
    integral = 0;
    for (int j=ilo*-1; j<=ihi; j++) { //loop over range surrounding peak
      Double_t le_bin = mean - (wvf->GetBinContent(peakbin+j) - bsl); //add the values to the histogram
      integral += le_bin;
    }
    integ4->Fill(integral); //add integral 
    //manualmode bounds
    ilo=7, ihi=10; //set bounds manually
    integral = 0;
    for (int j=ilo*-1; j<=ihi; j++) { //loop over range surrounding peak
      Double_t le_bin = mean - (wvf->GetBinContent(peakbin+j) - bsl); //add the values to the histogram
      integral += le_bin;
    }
    integ5->Fill(integral); //add integral
  }
}

  success++;
  cout << "  Analysis successful. " << nspe << " SPEs found." << endl;

} //end right channel/event if statement

else { //if it wasn't the right channel
  if (analysis_active && all_events) {eventid++;} //if we just finished an event (i.e. we got through all the waveforms of opc in eventid), switch to the next event
  analysis_active = false;
  //this lets us do all the events without needing to loop through the keys more than once!
}

} //end key loop

//NORMALIZE AVGSPE
for (int j=1;j<=NBINS;j++){
  Double_t le_bin=avgspe->GetBinContent(j); //get a bin
  Double_t norm = -1 * le_bin * (navspes-1) / navspes; //amount to subtract to reduce le_bin to le_bin/navspes
  avgspe->AddBinContent(j, norm); //add it
}

//FINAL PRINTOUT
cout << "======" << endl << "Analyses complete." << endl << navspes << " SPEs analyzed from " << success << " waveforms. Analysis failed on " << failed << " waveforms." << endl; 

//DRAW STATEMENTS
if (draw_avgspe) {avgspe->Draw();}
if (draw_amp) {amp->Draw();}
if (draw_integ) {integ1->Draw();}



//SAVE
if (save) {
  cout << "Saving histograms..." << endl;
  TFile *saved = new TFile(Form("ser_opchannel_%d.root", opc), filemode);
  //save the histograms!
  if (do_avgspe) {saved->WriteObject(avgspe, Form("avgspe_opchannel_%d", opc));}
  if (do_amp) {saved->WriteObject(amp, Form("amp_opchannel_%d", opc));}
  if (do_integ) {
    saved->WriteObject(integ0, Form("integ_opchannel_%d_zeromode", opc));
    saved->WriteObject(integ1, Form("integ_opchannel_%d_threshmode", opc));
    saved->WriteObject(integ2, Form("integ_opchannel_%d_manualmode", opc));
    saved->WriteObject(integ3, Form("integ_opchannel_%d_zeromodeB", opc));
    saved->WriteObject(integ4, Form("integ_opchannel_%d_threshmodeB", opc));
    saved->WriteObject(integ5, Form("integ_opchannel_%d_manualmodeB", opc));
  }
  cout << "Histograms saved!" << endl;
  myWvf->Close();
  myTruth->Close(); 
}
//SHOULD I BE CLOSING THE FILE SOMEWHERE?
timer.Stop();
printf("Elapsed time: %f seconds\n", timer.RealTime());
} //end function
