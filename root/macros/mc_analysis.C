// Macro to create histograms of MC events
// Author: Eric Pooser, pooser@jlab.org

// Declare data files and directories
TFile *mcRawFile;
TDirectory *mcRawDir;
// Declare trees
TTree *mcRawTree;
// Number of entries
Int_t    nEntriesMCRaw;
// Input MC variables
Float_t  xFocalMCRaw, xpFocalMCRaw, yFocalMCRaw, ypFocalMCRaw;
Float_t  yTarMCRaw, xpTarMCRaw, ypTarMCRaw, deltaMCRaw;
// Data variables
// Input MC histos
TH1F *h_xFocalMCRaw, *h_xpFocalMCRaw, *h_yFocalMCRaw, *h_ypFocalMCRaw;
TH1F *h_yTarMCRaw, *h_xpTarMCRaw, *h_ypTarMCRaw, *h_deltaMCRaw;
TH2F *h2_xVxpFocalMCRaw, *h2_xVyFocalMCRaw, *h2_xVypFocalMCRaw;
TH2F *h2_xpVyFocalMCRaw, *h2_xpVypFocalMCRaw, *h2_yVypFocalMCRaw;
TH2F *h2_yVxpTarMCRaw, *h2_yVypTarMCRaw, *h2_xpVypTarMCRaw;
// Declare constants
static const Double_t protonMass    = 0.938272;  // GeV
static const Double_t rad2mrad      = 1000.0;
// Define functions to calculate kinematic variables
Double_t calc_q2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle);
Double_t calc_w2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle);

void mc_analysis() {
  TH1::SetDefaultSumw2(false);
  // Open the ROOT files
  mcRawFile = new TFile("../input/kpp_shms_488.root");
  mcOutFile  = new TFile("../output/mc_488.root", "RECREATE");
  // Obtain the ROOT trees
  mcRawTree = dynamic_cast <TTree*> (mcRawFile->Get("h1411"));
    
  // Acquire the number of entries for each tree, and calculate scale factor
  nEntriesMCRaw = mcRawTree->GetEntries();
  
  // Acquire leafs of interest
  // Input MC leafs
  mcRawTree->SetBranchAddress("hsxfp",   &xFocalMCRaw);
  mcRawTree->SetBranchAddress("hsxpfp",  &xpFocalMCRaw);
  mcRawTree->SetBranchAddress("hsyfp",   &yFocalMCRaw);
  mcRawTree->SetBranchAddress("hsypfp",  &ypFocalMCRaw);
  mcRawTree->SetBranchAddress("hsytar",  &yTarMCRaw);
  mcRawTree->SetBranchAddress("hsxptar", &xpTarMCRaw);
  mcRawTree->SetBranchAddress("hsyptar", &ypTarMCRaw);
  mcRawTree->SetBranchAddress("hsdelta", &deltaMCRaw);
    
  // Create input MC directory and descend into it
  mcRawDir = dynamic_cast <TDirectory*> (mcOutFile->Get("mcRawDir"));
  if(!mcRawDir) {mcRawDir = mcOutFile->mkdir("mcRawDir"); mcRawDir->cd();}
  else mcOutFile->cd("mcRawDir");
  // Book input 1D MC histos
  h_xFocalMCRaw  = new TH1F("h_xFocalMCRaw",  "Input Monte-Carlo: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm",   160, -40, 40);
  h_xpFocalMCRaw = new TH1F("h_xpFocalMCRaw", "Input Monte-Carlo: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad",    60, -60.0, 60.0);
  h_yFocalMCRaw  = new TH1F("h_yFocalMCRaw",  "Input Monte-Carlo: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm",   160, -40, 40);
  h_ypFocalMCRaw = new TH1F("h_ypFocalMCRaw", "Input Monte-Carlo: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad",    60, -60.0, 60.0);
  h_yTarMCRaw    = new TH1F("h_yTarMCRaw",    "Input Monte-Carlo: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 100, -5, 5);
  h_xpTarMCRaw   = new TH1F("h_xpTarMCRaw",   "Input Monte-Carlo: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad",  60, -60.0, 60.0);
  h_ypTarMCRaw   = new TH1F("h_ypTarMCRaw",   "Input Monte-Carlo: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad",  60, -60.0, 60.0);
  h_deltaMCRaw   = new TH1F("h_deltaMCRaw",   "Input Monte-Carlo: #delta; #delta; Number of Entries",             80, -40, 40);
  // Book input 2D MC histos
  h2_xVxpFocalMCRaw  = new TH2F("h2_xVxpFocalMCRaw",  "Input Monte-Carlo: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    100, -100.0, 100.0, 160, -40, 40);
  h2_xVyFocalMCRaw   = new TH2F("h2_xVyFocalMCRaw",   "Input Monte-Carlo: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm",   160, -40, 40, 160, -40, 40);
  h2_xVypFocalMCRaw  = new TH2F("h2_xVypFocalMCRaw",  "Input Monte-Carlo: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_xpVyFocalMCRaw  = new TH2F("h2_xpVyFocalMCRaw",  "Input Monte-Carlo: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad",    160, -40, 40, 100, -100.0, 100.0);
  h2_xpVypFocalMCRaw = new TH2F("h2_xpVypFocalMCRaw", "Input Monte-Carlo: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad",     60, -60.0, 60.0, 100, -100.0, 100.0);
  h2_yVypFocalMCRaw  = new TH2F("h2_yVypFocalMCRaw",  "Input Monte-Carlo: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_yVxpTarMCRaw    = new TH2F("h2_yVxpTarMCRaw",    "Input Monte-Carlo: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_yVypTarMCRaw    = new TH2F("h2_yVypTarMCRaw",    "Input Monte-Carlo: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_xpVypTarMCRaw   = new TH2F("h2_xpVypTarMCRaw",   "Input Monte-Carlo: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad", 60, -60.0, 60.0, 100, -100.0, 100.0);
  mcOutFile->cd("../");

  UInt_t mcRawEventCntr = 0;
  // Loop over raw MC data
  for (UInt_t imcRaw = 0; imcRaw < nEntriesMCRaw; imcRaw++) {
    // Obtain the data entry
    mcRawTree->GetEntry(imcRaw);
    if (mcRawEventCntr % 10000 == 0 && mcRawEventCntr != 0) 
      cout << mcRawEventCntr << " Input Monte-Carlo events have been processed..." << endl;
    mcRawEventCntr++;
    // Fill histos
    // 1D histos
    h_xFocalMCRaw->Fill(xFocalMCRaw);
    h_xpFocalMCRaw->Fill(xpFocalMCRaw*rad2mrad);
    h_yFocalMCRaw->Fill(yFocalMCRaw);
    h_ypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad);
    h_yTarMCRaw->Fill(yTarMCRaw);
    h_xpTarMCRaw->Fill(xpTarMCRaw*rad2mrad);
    h_ypTarMCRaw->Fill(ypTarMCRaw*rad2mrad);
    h_deltaMCRaw->Fill(deltaMCRaw);
    // 2D histos
    h2_xVxpFocalMCRaw->Fill(xpFocalMCRaw*rad2mrad, xFocalMCRaw);
    h2_xVyFocalMCRaw->Fill(yFocalMCRaw, xFocalMCRaw);
    h2_xVypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad, xFocalMCRaw);
    h2_xpVyFocalMCRaw->Fill(yFocalMCRaw, xpFocalMCRaw*rad2mrad);
    h2_xpVypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad, xpFocalMCRaw*rad2mrad);
    h2_yVypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad, yFocalMCRaw);
    h2_yVxpTarMCRaw->Fill(xpTarMCRaw*rad2mrad, yTarMCRaw);
    h2_yVypTarMCRaw->Fill(ypTarMCRaw*rad2mrad, yTarMCRaw);
    h2_xpVypTarMCRaw->Fill(ypTarMCRaw*rad2mrad, xpTarMCRaw*rad2mrad);
  }  // Raw MC loop

  // Write comparison ROOT file
  mcOutFile->Write();
  // Close comparison ROOT file
  mcOutFile->Close();

}

// User specific funtions
Double_t calc_q2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle) {
  return 2.0*beamEnergy*scatMom*(1.0 - TMath::Cos(scatAngle));
}

Double_t calc_w2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle) {
  return protonMass*protonMass + 2.0*protonMass*(beamEnergy - scatMom) - 2.0*beamEnergy*scatMom*(1.0 - TMath::Cos(scatAngle));
}
