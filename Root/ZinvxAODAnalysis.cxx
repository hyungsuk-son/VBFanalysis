#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ZinvAnalysis/ZinvxAODAnalysis.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODBase/IParticleHelpers.h"
#include "AthContainers/ConstDataVector.h"

#include <TSystem.h>
#include <TFile.h>

#include "xAODRootAccess/tools/Message.h"

#include "PATInterfaces/CorrectionCode.h" // to check the return correction code status of tools
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/AuxContainerBase.h"

// header files for systematics:
#include "PATInterfaces/SystematicVariation.h" 
#include "PATInterfaces/SystematicsUtil.h"


static std::string jetType = "AntiKt4EMTopoJets";

// Global accessors and decorators
static SG::AuxElement::Decorator<char> dec_baseline("baseline");
static SG::AuxElement::Decorator<char> dec_signal("signal");
static SG::AuxElement::Decorator<char> dec_bad("bad");
static SG::AuxElement::Accessor<float>  acc_jvt("Jvt");
static SG::AuxElement::ConstAccessor<float> cacc_jvt("Jvt");
// For ORTools
static const std::string inputLabel = "selected";
static const std::string outputLabel = "overlaps";
const bool outputPassValue = false;
static const std::string bJetLabel = "";
//static SG::AuxElement::Accessor<char> overlapAcc("overlaps");
ort::inputAccessor_t selectAcc(inputLabel);
ort::inputDecorator_t selectDec(inputLabel);
ort::outputAccessor_t overlapAcc(outputLabel);
ort::inputDecorator_t bJetDec(bJetLabel);
ort::objLinkAccessor_t objLinkAcc("overlapObject");
//static const bool outputPassValue = false;
//static const std::string outputLabel = outputPassValue? "passOR" : "overlaps";
//static SG::AuxElement::Decorator<char> dec_overlap(outputLabel);

// Scale Factor decorators
static SG::AuxElement::Decorator<double> dec_scalefactor("scalefactor");

struct DescendingPt:std::function<bool(const xAOD::IParticle*, const xAOD::IParticle*)> {
  bool operator()(const xAOD::IParticle* l, const xAOD::IParticle* r)  const {
    return l->pt() > r->pt();
  }
};


// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
  do {                                                     \
    if( ! EXP.isSuccess() ) {                             \
      Error( CONTEXT,                                    \
          XAOD_MESSAGE( "Failed to execute: %s" ),    \
#EXP );                                     \
      return EL::StatusCode::FAILURE;                    \
    }                                                     \
  } while( false )


// this is needed to distribute the algorithm to the workers
ClassImp(ZinvxAODAnalysis)



ZinvxAODAnalysis :: ZinvxAODAnalysis ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode ZinvxAODAnalysis :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  job.useXAOD ();
  xAOD::Init().ignore(); // call before opening first file
  //  CP::CorrectionCode::enableFailure();
  EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.


  h_sumOfWeights = new TH1D("h_sumOfWeights", "MetaData_EventCount", 6, 0.5, 6.5);
  h_sumOfWeights -> GetXaxis() -> SetBinLabel(1, "sumOfWeights DxAOD");
  h_sumOfWeights -> GetXaxis() -> SetBinLabel(2, "sumOfWeightsSquared DxAOD");
  h_sumOfWeights -> GetXaxis() -> SetBinLabel(3, "nEvents DxAOD");
  h_sumOfWeights -> GetXaxis() -> SetBinLabel(4, "sumOfWeights initial");
  h_sumOfWeights -> GetXaxis() -> SetBinLabel(5, "sumOfWeightsSquared initial");
  h_sumOfWeights -> GetXaxis() -> SetBinLabel(6, "nEvents initial");
  wk()->addOutput (h_sumOfWeights);


  /*

  TH1::SetDefaultSumw2(kTRUE);


  // Zinv study
  // Znunu
  // MET
  h_znunu_met = new TH1F("h_znunu_met", "Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  wk()->addOutput (h_znunu_met);
  // Jets
  h_znunu_njet = new TH1F("h_znunu_njet", "Number of Jets;Numbers", 40, 0, 40);
  h_znunu_jet1_pt = new TH1F("h_znunu_jet1_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_znunu_jet2_pt = new TH1F("h_znunu_jet2_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_znunu_jet3_pt = new TH1F("h_znunu_jet3_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_znunu_jet1_phi = new TH1F("h_znunu_jet1_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_znunu_jet2_phi = new TH1F("h_znunu_jet2_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_znunu_jet3_phi = new TH1F("h_znunu_jet3_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_znunu_jet1_eta = new TH1F("h_znunu_jet1_eta", "Jet eta;#eta", 100, -5, 5);
  h_znunu_jet2_eta = new TH1F("h_znunu_jet2_eta", "Jet eta;#eta", 100, -5, 5);
  h_znunu_jet3_eta = new TH1F("h_znunu_jet3_eta", "Jet eta;#eta", 100, -5, 5);
  h_znunu_jet1_rap = new TH1F("h_znunu_jet1_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_znunu_jet2_rap = new TH1F("h_znunu_jet2_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_znunu_jet3_rap = new TH1F("h_znunu_jet3_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_znunu_mjj = new TH1F("h_znunu_mjj", "Invariant mass;m_{jj} (GeV)", 80, 0, 4000);
  h_znunu_dPhijj = new TH1F("h_znunu_dPhijj", "Delta Phi between J1 and J2;#Delta#phi(j_{1},j_{2}) (rad)", 32, 0, 3.2);
  h_znunu_dRjj = new TH1F("h_znunu_dRjj", "Delta R between J1 and J2;#DeltaR(j_{1},j_{2}) (rad)", 50, 0, 5);
  h_znunu_dPhimetj1 = new TH1F("h_znunu_dPhimetj1", "Delta Phi between MET and J1;#Delta#phi(MET,j_{1}) (rad)", 32, 0, 3.2);
  h_znunu_dPhimetj2 = new TH1F("h_znunu_dPhimetj2", "Delta Phi between MET and J2;#Delta#phi(MET,j_{2}) (rad)", 32, 0, 3.2);
  h_znunu_dPhimetj3 = new TH1F("h_znunu_dPhimetj3", "Delta Phi between MET and J3;#Delta#phi(MET,j_{3}) (rad)", 32, 0, 3.2);
  h_znunu_dPhiMinmetjet = new TH1F("h_znunu_dPhiMinmetjet", "Minimum Delta Phi between MET and Jet;#Delta#phi_{min}(MET,jet) (rad)", 32, 0, 3.2);
  wk()->addOutput (h_znunu_njet);
  wk()->addOutput (h_znunu_jet1_pt);
  wk()->addOutput (h_znunu_jet2_pt);
  wk()->addOutput (h_znunu_jet3_pt);
  wk()->addOutput (h_znunu_jet1_phi);
  wk()->addOutput (h_znunu_jet2_phi);
  wk()->addOutput (h_znunu_jet3_phi);
  wk()->addOutput (h_znunu_jet1_eta);
  wk()->addOutput (h_znunu_jet2_eta);
  wk()->addOutput (h_znunu_jet3_eta);
  wk()->addOutput (h_znunu_jet1_rap);
  wk()->addOutput (h_znunu_jet2_rap);
  wk()->addOutput (h_znunu_jet3_rap);
  wk()->addOutput (h_znunu_mjj);
  wk()->addOutput (h_znunu_dPhijj);
  wk()->addOutput (h_znunu_dRjj);
  wk()->addOutput (h_znunu_dPhimetj1);
  wk()->addOutput (h_znunu_dPhimetj2);
  wk()->addOutput (h_znunu_dPhimetj3);
  wk()->addOutput (h_znunu_dPhiMinmetjet);


  // Zmumu
  // MET
  h_zmumu_met = new TH1F("h_zmumu_met", "Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  h_zmumu_emulmet_nomu = new TH1F("h_zmumu_emulmet_nomu", "Emulated Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  wk()->addOutput (h_zmumu_met);
  wk()->addOutput (h_zmumu_emulmet_nomu);
  // Jets
  h_zmumu_njet = new TH1F("h_zmumu_njet", "Number of Jets;Numbers", 40, 0, 40);
  h_zmumu_jet1_pt = new TH1F("h_zmumu_jet1_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_zmumu_jet2_pt = new TH1F("h_zmumu_jet2_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_zmumu_jet3_pt = new TH1F("h_zmumu_jet3_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_zmumu_jet1_phi = new TH1F("h_zmumu_jet1_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zmumu_jet2_phi = new TH1F("h_zmumu_jet2_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zmumu_jet3_phi = new TH1F("h_zmumu_jet3_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zmumu_jet1_eta = new TH1F("h_zmumu_jet1_eta", "Jet eta;#eta", 100, -5, 5);
  h_zmumu_jet2_eta = new TH1F("h_zmumu_jet2_eta", "Jet eta;#eta", 100, -5, 5);
  h_zmumu_jet3_eta = new TH1F("h_zmumu_jet3_eta", "Jet eta;#eta", 100, -5, 5);
  h_zmumu_jet1_rap = new TH1F("h_zmumu_jet1_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_zmumu_jet2_rap = new TH1F("h_zmumu_jet2_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_zmumu_jet3_rap = new TH1F("h_zmumu_jet3_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_zmumu_mjj = new TH1F("h_zmumu_mjj", "Invariant mass;m_{jj} (GeV)", 80, 0, 4000);
  h_zmumu_dPhijj = new TH1F("h_zmumu_dPhijj", "Delta Phi between J1 and J2;#Delta#phi(j_{1},j_{2}) (rad)", 32, 0, 3.2);
  h_zmumu_dRjj = new TH1F("h_zmumu_dRjj", "Delta R between J1 and J2;#DeltaR(j_{1},j_{2}) (rad)", 50, 0, 5);
  h_zmumu_dPhimetj1 = new TH1F("h_zmumu_dPhimetj1", "Delta Phi between MET and J1;#Delta#phi(MET,j_{1}) (rad)", 32, 0, 3.2);
  h_zmumu_dPhimetj2 = new TH1F("h_zmumu_dPhimetj2", "Delta Phi between MET and J2;#Delta#phi(MET,j_{2}) (rad)", 32, 0, 3.2);
  h_zmumu_dPhimetj3 = new TH1F("h_zmumu_dPhimetj3", "Delta Phi between MET and J3;#Delta#phi(MET,j_{3}) (rad)", 32, 0, 3.2);
  h_zmumu_dPhiMinmetjet = new TH1F("h_zmumu_dPhiMinmetjet", "Minimum Delta Phi between MET and Jet;#Delta#phi_{min}(MET,jet) (rad)", 32, 0, 3.2);
  wk()->addOutput (h_zmumu_njet);
  wk()->addOutput (h_zmumu_jet1_pt);
  wk()->addOutput (h_zmumu_jet2_pt);
  wk()->addOutput (h_zmumu_jet3_pt);
  wk()->addOutput (h_zmumu_jet1_phi);
  wk()->addOutput (h_zmumu_jet2_phi);
  wk()->addOutput (h_zmumu_jet3_phi);
  wk()->addOutput (h_zmumu_jet1_eta);
  wk()->addOutput (h_zmumu_jet2_eta);
  wk()->addOutput (h_zmumu_jet3_eta);
  wk()->addOutput (h_zmumu_jet1_rap);
  wk()->addOutput (h_zmumu_jet2_rap);
  wk()->addOutput (h_zmumu_jet3_rap);
  wk()->addOutput (h_zmumu_mjj);
  wk()->addOutput (h_zmumu_dPhijj);
  wk()->addOutput (h_zmumu_dRjj);
  wk()->addOutput (h_zmumu_dPhimetj1);
  wk()->addOutput (h_zmumu_dPhimetj2);
  wk()->addOutput (h_zmumu_dPhimetj3);
  wk()->addOutput (h_zmumu_dPhiMinmetjet);
  // Leptons
  h_zmumu_muon1_pt = new TH1F("h_zmumu_muon1_pt", "Muon p_{T};p_{T} (GeV)", 30, 0, 1500); // [GeV]
  h_zmumu_muon2_pt = new TH1F("h_zmumu_muon2_pt", "Muon p_{T};p_{T} (GeV)", 30, 0, 1500); // [GeV]
  h_zmumu_muon1_phi = new TH1F("h_zmumu_muon1_phi", "Muon Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zmumu_muon2_phi = new TH1F("h_zmumu_muon2_phi", "Muon Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zmumu_muon1_eta = new TH1F("h_zmumu_muon1_eta", "Muon eta;#eta", 100, -5, 5);
  h_zmumu_muon2_eta = new TH1F("h_zmumu_muon2_eta", "Muon eta;#eta", 100, -5, 5);
  h_zmumu_mll = new TH1F("h_zmumu_mll", "Invariant mass;m_{ll} (GeV)", 150, 0, 300);
  wk()->addOutput (h_zmumu_muon1_pt);
  wk()->addOutput (h_zmumu_muon2_pt);
  wk()->addOutput (h_zmumu_muon1_phi);
  wk()->addOutput (h_zmumu_muon2_phi);
  wk()->addOutput (h_zmumu_muon1_eta);
  wk()->addOutput (h_zmumu_muon2_eta);
  wk()->addOutput (h_zmumu_mll);


  // Wmunu
  // MET
  h_wmunu_met = new TH1F("h_wmunu_met", "Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  h_wmunu_emulmet_nomu = new TH1F("h_wmunu_emulmet_nomu", "Emulated Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  wk()->addOutput (h_wmunu_met);
  wk()->addOutput (h_wmunu_emulmet_nomu);
  // Jets
  h_wmunu_njet = new TH1F("h_wmunu_njet", "Number of Jets;Numbers", 40, 0, 40);
  h_wmunu_jet1_pt = new TH1F("h_wmunu_jet1_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_wmunu_jet2_pt = new TH1F("h_wmunu_jet2_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_wmunu_jet3_pt = new TH1F("h_wmunu_jet3_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_wmunu_jet1_phi = new TH1F("h_wmunu_jet1_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wmunu_jet2_phi = new TH1F("h_wmunu_jet2_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wmunu_jet3_phi = new TH1F("h_wmunu_jet3_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wmunu_jet1_eta = new TH1F("h_wmunu_jet1_eta", "Jet eta;#eta", 100, -5, 5);
  h_wmunu_jet2_eta = new TH1F("h_wmunu_jet2_eta", "Jet eta;#eta", 100, -5, 5);
  h_wmunu_jet3_eta = new TH1F("h_wmunu_jet3_eta", "Jet eta;#eta", 100, -5, 5);
  h_wmunu_jet1_rap = new TH1F("h_wmunu_jet1_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_wmunu_jet2_rap = new TH1F("h_wmunu_jet2_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_wmunu_jet3_rap = new TH1F("h_wmunu_jet3_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_wmunu_mjj = new TH1F("h_wmunu_mjj", "Invariant mass;m_{jj} (GeV)", 80, 0, 4000);
  h_wmunu_dPhijj = new TH1F("h_wmunu_dPhijj", "Delta Phi between J1 and J2;#Delta#phi(j_{1},j_{2}) (rad)", 32, 0, 3.2);
  h_wmunu_dRjj = new TH1F("h_wmunu_dRjj", "Delta R between J1 and J2;#DeltaR(j_{1},j_{2}) (rad)", 50, 0, 5);
  h_wmunu_dPhimetj1 = new TH1F("h_wmunu_dPhimetj1", "Delta Phi between MET and J1;#Delta#phi(MET,j_{1}) (rad)", 32, 0, 3.2);
  h_wmunu_dPhimetj2 = new TH1F("h_wmunu_dPhimetj2", "Delta Phi between MET and J2;#Delta#phi(MET,j_{2}) (rad)", 32, 0, 3.2);
  h_wmunu_dPhimetj3 = new TH1F("h_wmunu_dPhimetj3", "Delta Phi between MET and J3;#Delta#phi(MET,j_{3}) (rad)", 32, 0, 3.2);
  h_wmunu_dPhiMinmetjet = new TH1F("h_wmunu_dPhiMinmetjet", "Minimum Delta Phi between MET and Jet;#Delta#phi_{min}(MET,jet) (rad)", 32, 0, 3.2);
  wk()->addOutput (h_wmunu_njet);
  wk()->addOutput (h_wmunu_jet1_pt);
  wk()->addOutput (h_wmunu_jet2_pt);
  wk()->addOutput (h_wmunu_jet3_pt);
  wk()->addOutput (h_wmunu_jet1_phi);
  wk()->addOutput (h_wmunu_jet2_phi);
  wk()->addOutput (h_wmunu_jet3_phi);
  wk()->addOutput (h_wmunu_jet1_eta);
  wk()->addOutput (h_wmunu_jet2_eta);
  wk()->addOutput (h_wmunu_jet3_eta);
  wk()->addOutput (h_wmunu_jet1_rap);
  wk()->addOutput (h_wmunu_jet2_rap);
  wk()->addOutput (h_wmunu_jet3_rap);
  wk()->addOutput (h_wmunu_mjj);
  wk()->addOutput (h_wmunu_dPhijj);
  wk()->addOutput (h_wmunu_dRjj);
  wk()->addOutput (h_wmunu_dPhimetj1);
  wk()->addOutput (h_wmunu_dPhimetj2);
  wk()->addOutput (h_wmunu_dPhimetj3);
  wk()->addOutput (h_wmunu_dPhiMinmetjet);
  // Leptons
  h_wmunu_muon_pt = new TH1F("h_wmunu_muon_pt", "Muon p_{T};p_{T} (GeV)", 30, 0, 1500); // [GeV]
  h_wmunu_muon_phi = new TH1F("h_wmunu_muon_phi", "Muon Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wmunu_muon_eta = new TH1F("h_wmunu_muon_eta", "Muon eta;#eta", 100, -5, 5);
  h_wmunu_mT = new TH1F("h_wmunu_mT", "Transverse mass;m_{T} (GeV)", 75, 0, 150);
  wk()->addOutput (h_wmunu_muon_pt);
  wk()->addOutput (h_wmunu_muon_phi);
  wk()->addOutput (h_wmunu_muon_eta);
  wk()->addOutput (h_wmunu_mT);



  // Zee
  // MET
  h_zee_met = new TH1F("h_zee_met", "Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  h_zee_emulmet_noelec = new TH1F("h_zee_emulmet_noelec", "Emulated Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  wk()->addOutput (h_zee_met);
  wk()->addOutput (h_zee_emulmet_noelec);
  // Jets
  h_zee_njet = new TH1F("h_zee_njet", "Number of Jets;Numbers", 40, 0, 40);
  h_zee_jet1_pt = new TH1F("h_zee_jet1_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_zee_jet2_pt = new TH1F("h_zee_jet2_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_zee_jet3_pt = new TH1F("h_zee_jet3_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_zee_jet1_phi = new TH1F("h_zee_jet1_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zee_jet2_phi = new TH1F("h_zee_jet2_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zee_jet3_phi = new TH1F("h_zee_jet3_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zee_jet1_eta = new TH1F("h_zee_jet1_eta", "Jet eta;#eta", 100, -5, 5);
  h_zee_jet2_eta = new TH1F("h_zee_jet2_eta", "Jet eta;#eta", 100, -5, 5);
  h_zee_jet3_eta = new TH1F("h_zee_jet3_eta", "Jet eta;#eta", 100, -5, 5);
  h_zee_jet1_rap = new TH1F("h_zee_jet1_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_zee_jet2_rap = new TH1F("h_zee_jet2_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_zee_jet3_rap = new TH1F("h_zee_jet3_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_zee_mjj = new TH1F("h_zee_mjj", "Invariant mass;m_{jj} (GeV)", 80, 0, 4000);
  h_zee_dPhijj = new TH1F("h_zee_dPhijj", "Delta Phi between J1 and J2;#Delta#phi(j_{1},j_{2}) (rad)", 32, 0, 3.2);
  h_zee_dRjj = new TH1F("h_zee_dRjj", "Delta R between J1 and J2;#DeltaR(j_{1},j_{2}) (rad)", 50, 0, 5);
  h_zee_dPhimetj1 = new TH1F("h_zee_dPhimetj1", "Delta Phi between MET and J1;#Delta#phi(MET,j_{1}) (rad)", 32, 0, 3.2);
  h_zee_dPhimetj2 = new TH1F("h_zee_dPhimetj2", "Delta Phi between MET and J2;#Delta#phi(MET,j_{2}) (rad)", 32, 0, 3.2);
  h_zee_dPhimetj3 = new TH1F("h_zee_dPhimetj3", "Delta Phi between MET and J3;#Delta#phi(MET,j_{3}) (rad)", 32, 0, 3.2);
  h_zee_dPhiMinmetjet = new TH1F("h_zee_dPhiMinmetjet", "Minimum Delta Phi between MET and Jet;#Delta#phi_{min}(MET,jet) (rad)", 32, 0, 3.2);
  wk()->addOutput (h_zee_njet);
  wk()->addOutput (h_zee_jet1_pt);
  wk()->addOutput (h_zee_jet2_pt);
  wk()->addOutput (h_zee_jet3_pt);
  wk()->addOutput (h_zee_jet1_phi);
  wk()->addOutput (h_zee_jet2_phi);
  wk()->addOutput (h_zee_jet3_phi);
  wk()->addOutput (h_zee_jet1_eta);
  wk()->addOutput (h_zee_jet2_eta);
  wk()->addOutput (h_zee_jet3_eta);
  wk()->addOutput (h_zee_jet1_rap);
  wk()->addOutput (h_zee_jet2_rap);
  wk()->addOutput (h_zee_jet3_rap);
  wk()->addOutput (h_zee_mjj);
  wk()->addOutput (h_zee_dPhijj);
  wk()->addOutput (h_zee_dRjj);
  wk()->addOutput (h_zee_dPhimetj1);
  wk()->addOutput (h_zee_dPhimetj2);
  wk()->addOutput (h_zee_dPhimetj3);
  wk()->addOutput (h_zee_dPhiMinmetjet);
  // Leptons
  h_zee_electron1_pt = new TH1F("h_zee_electron1_pt", "Electron p_{T};p_{T} (GeV)", 30, 0, 1500); // [GeV]
  h_zee_electron2_pt = new TH1F("h_zee_electron2_pt", "Electron p_{T};p_{T} (GeV)", 30, 0, 1500); // [GeV]
  h_zee_electron1_phi = new TH1F("h_zee_electron1_phi", "Electron Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zee_electron2_phi = new TH1F("h_zee_electron2_phi", "Electron Phi;#phi (rad)", 64, -3.2, 3.2);
  h_zee_electron1_eta = new TH1F("h_zee_electron1_eta", "Electron eta;#eta", 100, -5, 5);
  h_zee_electron2_eta = new TH1F("h_zee_electron2_eta", "Electron eta;#eta", 100, -5, 5);
  h_zee_mll = new TH1F("h_zee_mll", "Invariant mass;m_{ll} (GeV)", 150, 0, 300);
  wk()->addOutput (h_zee_electron1_pt);
  wk()->addOutput (h_zee_electron2_pt);
  wk()->addOutput (h_zee_electron1_phi);
  wk()->addOutput (h_zee_electron2_phi);
  wk()->addOutput (h_zee_electron1_eta);
  wk()->addOutput (h_zee_electron2_eta);
  wk()->addOutput (h_zee_mll);


  // Wenu
  // MET
  h_wenu_met = new TH1F("h_wenu_met", "Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  h_wenu_emulmet_noelec = new TH1F("h_wenu_emulmet_noelec", "Emulated Offline |Missing E_{T}|;ME_{T} (GeV)", 30, 0, 1500);
  wk()->addOutput (h_wenu_met);
  wk()->addOutput (h_wenu_emulmet_noelec);
  // Jets
  h_wenu_njet = new TH1F("h_wenu_njet", "Number of Jets;Numbers", 40, 0, 40);
  h_wenu_jet1_pt = new TH1F("h_wenu_jet1_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_wenu_jet2_pt = new TH1F("h_wenu_jet2_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_wenu_jet3_pt = new TH1F("h_wenu_jet3_pt", "Jet p_{T};p_{T} (GeV)", 60, 0, 3000); // Jet pt [GeV]
  h_wenu_jet1_phi = new TH1F("h_wenu_jet1_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wenu_jet2_phi = new TH1F("h_wenu_jet2_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wenu_jet3_phi = new TH1F("h_wenu_jet3_phi", "Jet Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wenu_jet1_eta = new TH1F("h_wenu_jet1_eta", "Jet eta;#eta", 100, -5, 5);
  h_wenu_jet2_eta = new TH1F("h_wenu_jet2_eta", "Jet eta;#eta", 100, -5, 5);
  h_wenu_jet3_eta = new TH1F("h_wenu_jet3_eta", "Jet eta;#eta", 100, -5, 5);
  h_wenu_jet1_rap = new TH1F("h_wenu_jet1_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_wenu_jet2_rap = new TH1F("h_wenu_jet2_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_wenu_jet3_rap = new TH1F("h_wenu_jet3_rap", "Jet rapidity;#eta", 100, -5, 5);
  h_wenu_mjj = new TH1F("h_wenu_mjj", "Invariant mass;m_{jj} (GeV)", 80, 0, 4000);
  h_wenu_dPhijj = new TH1F("h_wenu_dPhijj", "Delta Phi between J1 and J2;#Delta#phi(j_{1},j_{2}) (rad)", 32, 0, 3.2);
  h_wenu_dRjj = new TH1F("h_wenu_dRjj", "Delta R between J1 and J2;#DeltaR(j_{1},j_{2}) (rad)", 50, 0, 5);
  h_wenu_dPhimetj1 = new TH1F("h_wenu_dPhimetj1", "Delta Phi between MET and J1;#Delta#phi(MET,j_{1}) (rad)", 32, 0, 3.2);
  h_wenu_dPhimetj2 = new TH1F("h_wenu_dPhimetj2", "Delta Phi between MET and J2;#Delta#phi(MET,j_{2}) (rad)", 32, 0, 3.2);
  h_wenu_dPhimetj3 = new TH1F("h_wenu_dPhimetj3", "Delta Phi between MET and J3;#Delta#phi(MET,j_{3}) (rad)", 32, 0, 3.2);
  h_wenu_dPhiMinmetjet = new TH1F("h_wenu_dPhiMinmetjet", "Minimum Delta Phi between MET and Jet;#Delta#phi_{min}(MET,jet) (rad)", 32, 0, 3.2);
  wk()->addOutput (h_wenu_njet);
  wk()->addOutput (h_wenu_jet1_pt);
  wk()->addOutput (h_wenu_jet2_pt);
  wk()->addOutput (h_wenu_jet3_pt);
  wk()->addOutput (h_wenu_jet1_phi);
  wk()->addOutput (h_wenu_jet2_phi);
  wk()->addOutput (h_wenu_jet3_phi);
  wk()->addOutput (h_wenu_jet1_eta);
  wk()->addOutput (h_wenu_jet2_eta);
  wk()->addOutput (h_wenu_jet3_eta);
  wk()->addOutput (h_wenu_jet1_rap);
  wk()->addOutput (h_wenu_jet2_rap);
  wk()->addOutput (h_wenu_jet3_rap);
  wk()->addOutput (h_wenu_mjj);
  wk()->addOutput (h_wenu_dPhijj);
  wk()->addOutput (h_wenu_dRjj);
  wk()->addOutput (h_wenu_dPhimetj1);
  wk()->addOutput (h_wenu_dPhimetj2);
  wk()->addOutput (h_wenu_dPhimetj3);
  wk()->addOutput (h_wenu_dPhiMinmetjet);
  // Leptons
  h_wenu_electron_pt = new TH1F("h_wenu_electron_pt", "Electron p_{T};p_{T} (GeV)", 30, 0, 1500); // [GeV]
  h_wenu_electron_phi = new TH1F("h_wenu_electron_phi", "Electron Phi;#phi (rad)", 64, -3.2, 3.2);
  h_wenu_electron_eta = new TH1F("h_wenu_electron_eta", "Electron eta;#eta", 100, -5, 5);
  h_wenu_mT = new TH1F("h_wenu_mT", "Transverse mass;m_{T} (GeV)", 75, 0, 150);
  wk()->addOutput (h_wenu_electron_pt);
  wk()->addOutput (h_wenu_electron_phi);
  wk()->addOutput (h_wenu_electron_eta);
  wk()->addOutput (h_wenu_mT);

  */


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  // get TEvent and TStore - must be done here b/c we need to retrieve CutBookkeepers container from TEvent!
  m_event = wk()->xaodEvent();

  // Event Bookkeepers
  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AnalysisMetadata#Luminosity_Bookkeepers

  // get the MetaData tree once a new file is opened, with
  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (!MetaData) {
    Error("fileExecute()", "MetaData not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  MetaData->LoadTree(0);

  //check if file is from a DxAOD
  bool m_isDerivation = !MetaData->GetBranch("StreamAOD");

  if(m_isDerivation){

    // check for corruption
    const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
    if(!m_event->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess()){
      Error("initializeEvent()","Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    if ( incompleteCBC->size() != 0 ) {
      Error("initializeEvent()","Found incomplete Bookkeepers! Check file for corruption.");
      return EL::StatusCode::FAILURE;
    }

    // Now, let's find the actual information
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if(!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()){
      Error("initializeEvent()","Failed to retrieve CutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }

    // First, let's find the smallest cycle number,
    // i.e., the original first processing step/cycle
    int minCycle = 10000;
    for ( auto cbk : *completeCBC ) {
      if ( ! cbk->name().empty()  && minCycle > cbk->cycle() ){ minCycle = cbk->cycle(); }
    }

    // Now, let's actually find the right one that contains all the needed info...
    const xAOD::CutBookkeeper* allEventsCBK=0;
    const xAOD::CutBookkeeper* DxAODEventsCBK=0;
    std::string derivationName = "EXOT5Kernel"; //need to replace by appropriate name
    int maxCycle = -1;
    for (const auto& cbk: *completeCBC) {
      if (cbk->cycle() > maxCycle && cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD") {
        allEventsCBK = cbk;
        maxCycle = cbk->cycle();
      }
      if ( cbk->name() == derivationName){
        DxAODEventsCBK = cbk;
      }
    }
    uint64_t nEventsProcessed  = allEventsCBK->nAcceptedEvents();
    double sumOfWeights        = allEventsCBK->sumOfEventWeights();
    double sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();

    uint64_t nEventsDxAOD           = DxAODEventsCBK->nAcceptedEvents();
    double sumOfWeightsDxAOD        = DxAODEventsCBK->sumOfEventWeights();
    double sumOfWeightsSquaredDxAOD = DxAODEventsCBK->sumOfEventWeightsSquared();



    //Info("execute()", " Event # = %llu, sumOfweights = %f, mcEventWeight = %f", eventInfo->eventNumber(), sumOfWeights, eventInfo->mcEventWeight());
    //Info("execute()", " Event # = %llu, nEventsProcessed = %d, sumOfweights = %f, sumOfWeightsSquared = %f", eventInfo->eventNumber(), nEventsProcessed, sumOfWeights, sumOfWeightsSquared);

    h_sumOfWeights -> Fill(1, sumOfWeightsDxAOD);
    h_sumOfWeights -> Fill(2, sumOfWeightsSquaredDxAOD);
    h_sumOfWeights -> Fill(3, nEventsDxAOD);
    h_sumOfWeights -> Fill(4, sumOfWeights);
    h_sumOfWeights -> Fill(5, sumOfWeightsSquared);
    h_sumOfWeights -> Fill(6, nEventsProcessed);

    //Info("execute()", " Event # = %llu, sumOfWeights/nEventsProcessed = %f", eventInfo->eventNumber(), sumOfWeights/double(nEventsProcessed));
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  //----------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection in initialise. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // check if the event is data or MC
  // (many tools are applied either to data or MC)
  m_isData = true;
  // check if the event is MC
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    m_isData = false; // can do something with this later
  }

  // Retrieve MC channel number
  if (m_isData)
    mcChannelNumber = 1;
  else mcChannelNumber = eventInfo->mcChannelNumber();

  // count number of events
  m_eventCounter = 0;
  m_numCleanEvents = 0;

  // Enable Cutflow plot
  m_useArrayCutflow = false;
  m_useBitsetCutflow = true;
  m_isEmilyCutflow = false;

  // Event Channel
  m_isZnunu = false;
  m_isZmumu = true;
  m_isWmunu = false;
  m_isZee = true;
  m_isWenu = false;

  // Enable Systematics
  m_doSys = true;

  // Enable Overlap Removal tool
  m_doORtool = false;
  m_doORmanual = true;

  // Cut values
  m_muonPtCut = 7.; /// GeV
  m_muonEtaCut = 2.5;
  m_elecPtCut = 7.; /// GeV
  m_elecEtaCut = 2.47;
  m_photPtCut = 20.; /// GeV
  m_photEtaCut = 2.47;
  m_jetPtCut = 20.; /// GeV
  m_jetEtaCut = 4.5;
  m_monoJetPtCut = 120.; /// GeV
  m_monoJetEtaCut = 2.4;
  m_diJet1PtCut = 80.; /// GeV
  m_diJet2PtCut = 50.; /// GeV
  m_diJetRapCut = 4.4;
  m_CJVptCut = 25.; ///GeV
  m_metCut = 200.; ///GeV
  m_mjjCut = 200.; ///GeV
  m_LeadLepPtCut = 80.; ///GeV
  m_SubLeadLepPtCut = 7.; ///GeV
  m_ORJETdeltaR = 0.2;
  m_isoMuonPtMin = 10.; ///GeV
  m_isoMuonPtMax = 500.; ///GeV
  m_recoSF = true;
  m_idSF = true;
  m_isoSF = true;
  m_ttvaSF = true; // for muon
  m_trigSF = true; // for electron
  m_isoPtCut = false; // Initialize Setting (Do not change), This bool is for if you want to set pt range for muon Isolation

  // GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  std::vector<std::string> vecStringGRL;
  // GRL xml file should be put in ZinvAnalysis/share directory
  vecStringGRL.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/ZinvAnalysis/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml"));
  EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
  EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
  EL_RETURN_CHECK("initialize()",m_grl->initialize());

  // Initialize and configure trigger tools
  m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
  EL_RETURN_CHECK( "initialize", m_trigConfigTool->initialize() );
  ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
  m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ) ); // connect the TrigDecisionTool to the ConfigTool
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision" ) );
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->initialize() );

  // initialize the muon calibration and smearing tool
  m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" );
  //m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::DEBUG );
  m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::INFO );
  EL_RETURN_CHECK("initialize()",m_muonCalibrationAndSmearingTool->initialize());

  // initialize the electron and photon calibration and smearing tool
  m_egammaCalibrationAndSmearingTool = new CP::EgammaCalibrationAndSmearingTool( "EgammaCorrectionTool" );
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "ESModel", "es2015PRE" ));  // see below for options
  //EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "decorrelationModel", "FULL_ETACORRELATED_v1" ));  // see below for options
  //EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "decorrelationModel", "FULL_v1" ));  // see below for options
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "decorrelationModel", "1NP_v1" ));  // see below for options
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->initialize());

  // Initialize the MC fudge tool
  m_electronPhotonShowerShapeFudgeTool = new ElectronPhotonShowerShapeFudgeTool( "ElectronPhotonShowerShapeFudgeTool" );
  int FFset = 16; // for MC15 samples, which are based on a geometry derived from GEO-21
  EL_RETURN_CHECK("initialize()",m_electronPhotonShowerShapeFudgeTool->setProperty("Preselection", FFset));
  EL_RETURN_CHECK("initialize()",m_electronPhotonShowerShapeFudgeTool->initialize() );

  // Muon identification (Medium)
  // initialize the muon selection tool
  m_muonSelection = new CP::MuonSelectionTool( "MuonSelection" );
  //EL_RETURN_CHECK("initialize()",m_muonSelection->setProperty( "MaxEta", 2.5 ));
  EL_RETURN_CHECK("initialize()",m_muonSelection->setProperty( "MuQuality", 1)); // 0 tight, 1 medium, 2 loose, 3 very loose
  //m_muonSelection->msg().setLevel( MSG::VERBOSE );
  m_muonSelection->msg().setLevel( MSG::INFO );
  //m_muonSelection->msg().setLevel( MSG::ERROR );
  EL_RETURN_CHECK("initialize()",m_muonSelection->initialize());
  // Muon identification (Loose)
  m_loosemuonSelection = new CP::MuonSelectionTool( "MuonLooseSelection" );
  //m_loosemuonSelection->msg().setLevel( MSG::VERBOSE );
  m_loosemuonSelection->msg().setLevel( MSG::INFO );
  //m_loosemuonSelection->msg().setLevel( MSG::ERROR );
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->setProperty( "MaxEta", 2.5 ));
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->setProperty( "MuQuality", 2));
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->initialize());

  // LH Electron identification
  // initialize the electron selection tool
  m_LHToolTight2015    = new AsgElectronLikelihoodTool ("m_LHToolTight2015");
  m_LHToolMedium2015   = new AsgElectronLikelihoodTool ("m_LHToolMedium2015"); 
  m_LHToolLoose2015    = new AsgElectronLikelihoodTool ("m_LHToolLoose2015");
  // initialize the primary vertex container for the tool to have access to the number of vertices used to adapt cuts based on the pileup
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  // define the config files
  std::string confDir_2015 = "ElectronPhotonSelectorTools/offline/mc15_20150712/";
  std::string confDir_2016 = "ElectronPhotonSelectorTools/offline/mc15_20160113/";
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->setProperty("ConfigFile",confDir_2016+"ElectronLikelihoodTightOfflineConfig2015.conf"));
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->setProperty("ConfigFile",confDir_2015+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
  //EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("ConfigFile",confDir_2015+"ElectronLikelihoodLooseOfflineConfig2015.conf"));
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("ConfigFile",confDir_2015+"ElectronLikelihoodLooseOfflineConfig2015_CutBL.conf"));
  // initialize
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->initialize());
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->initialize());
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->initialize());

  // Photon identification (Tight)
  // create the selector
  m_photonTightIsEMSelector = new AsgPhotonIsEMSelector ( "PhotonTightIsEMSelector" );
  // decide which kind of selection (Loose/Medium/Tight) you want to use
  EL_RETURN_CHECK("initialize()",m_photonTightIsEMSelector->setProperty("isEMMask",egammaPID::PhotonTight));
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  EL_RETURN_CHECK("initialize()",m_photonTightIsEMSelector->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf"));
  // initialise the tool
  EL_RETURN_CHECK("initialize()",m_photonTightIsEMSelector->initialize());
  // Photon identification (Medium)
  // create the selector
  m_photonMediumIsEMSelector = new AsgPhotonIsEMSelector ( "PhotonMediumIsEMSelector" );
  // decide which kind of selection (Loose/Medium/Tight) you want to use
  EL_RETURN_CHECK("initialize()",m_photonMediumIsEMSelector->setProperty("isEMMask",egammaPID::PhotonMedium));
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  EL_RETURN_CHECK("initialize()",m_photonMediumIsEMSelector->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMMediumSelectorCutDefs.conf"));
  // initialise the tool
  EL_RETURN_CHECK("initialize()",m_photonMediumIsEMSelector->initialize());
  // Photon identification (Loose)
  // create the selector
  m_photonLooseIsEMSelector = new AsgPhotonIsEMSelector ( "PhotonLooseIsEMSelector" );
  // decide which kind of selection (Loose/Medium/Tight) you want to use
  EL_RETURN_CHECK("initialize()",m_photonLooseIsEMSelector->setProperty("isEMMask",egammaPID::PhotonLoose));
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  EL_RETURN_CHECK("initialize()",m_photonLooseIsEMSelector->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMLooseSelectorCutDefs.conf"));
  // initialise the tool
  EL_RETURN_CHECK("initialize()",m_photonLooseIsEMSelector->initialize());



  // IsolationSelectionTool
  m_IsolationSelectionTool = new CP::IsolationSelectionTool("IsolationSelectionTool");
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->setProperty("MuonWP","Gradient"));
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->setProperty("ElectronWP","Gradient"));
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->setProperty("PhotonWP","Cone40"));
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->initialize());
  // IsolationSelectionTool for VBF signal
  m_IsoToolVBF = new CP::IsolationSelectionTool("IsoToolVBF");
  //EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("MuonWP","FixedCutLoose"));
  //EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("ElectronWP","FixedCutLoose"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("MuonWP","LooseTrackOnly"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("ElectronWP","LooseTrackOnly"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("PhotonWP","Cone40"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->initialize());


  // Tau identification
  // initialize the tau selection tool
  m_tauSelTool = new TauAnalysisTools::TauSelectionTool( "TauSelectionTool" );
  // define the config files
  std::string confPath = "$ROOTCOREBIN/data/ZinvAnalysis/";
  EL_RETURN_CHECK("initialize()",m_tauSelTool->setProperty( "ConfigPath", confPath+"recommended_selection_mc15.conf"));
  m_tauSelTool->msg().setLevel( MSG::INFO );
  //m_tauSelTool->msg().setLevel( MSG::DEBUG );
  // initialize
  EL_RETURN_CHECK("initialize()",m_tauSelTool->initialize());

  // initialize the tau selection tool for VBF analysis
  m_tauSelToolVBF = new TauAnalysisTools::TauSelectionTool( "TauSelectionToolVBF" );
  // define the config files
  EL_RETURN_CHECK("initialize()",m_tauSelToolVBF->setProperty( "ConfigPath", confPath+"recommended_selection_mc15_VBF.conf"));
  m_tauSelToolVBF->msg().setLevel( MSG::INFO );
  // initialize
  EL_RETURN_CHECK("initialize()",m_tauSelToolVBF->initialize());

  // Initialise tau smearing tool
  m_tauSmearingTool = new TauAnalysisTools::TauSmearingTool( "TauSmaringTool" );
  m_tauSmearingTool->msg().setLevel( MSG::INFO );
  // initialize
  EL_RETURN_CHECK("initialize()",m_tauSmearingTool->initialize());

  // Initialise TauOverlappingElectronLLHDecorator
  m_tauOverlappingElectronLLHDecorator = new TauAnalysisTools::TauOverlappingElectronLLHDecorator("TauOverlappingElectronLLHDecorator"); 
  EL_RETURN_CHECK("initialize()",m_tauOverlappingElectronLLHDecorator->initialize());


  // Jet
  // JES Calibration (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JES_calibration_AN1)
  const std::string name = "ZinvxAODAnalysis"; //string describing the current thread, for logging
  TString jetAlgo = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo
  //TString config = "JES_MC15Prerecommendation_April2015.config"; //Path to global config used to initialize the tool
  TString config = "JES_2015dataset_recommendation_Feb2016.config"; //Path to global config used to initialize the tool
  TString calibSeq = "JetArea_Residual_Origin_EtaJES_GSC"; //String describing the calibration sequence to apply
  if (m_isData) calibSeq += "_Insitu";
  //Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
  m_jetCalibration = new JetCalibrationTool(name, jetAlgo, config, calibSeq, m_isData);
  //Initialize the tool
  EL_RETURN_CHECK("initialize()",m_jetCalibration->initializeTool(name));

  // JES uncertainty (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JES_uncertainty)
  m_jetUncertaintiesTool = new JetUncertaintiesTool("JetUncertaintiesTool");
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("JetDefinition", "AntiKt4EMTopo"));
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("MCType", "MC15"));
  //EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("ConfigFile", "JES_2015/Prerec/PrerecJES2015_AllNuisanceParameters_25ns.config"));
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("ConfigFile", "JES_2015/Moriond2016/JES2015_SR_Scenario1.config"));
  // Initialise jet uncertainty tool
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->initialize());

  // JER uncertainty  (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JER_uncertainty)
  // Jet Resolution (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetResolution2015Prerecom)
  // Configure the JERTool
  m_jerTool = new JERTool("JERTool");
  EL_RETURN_CHECK("initialize()",m_jerTool->setProperty("PlotFileName", "JetResolution/Prerec2015_xCalib_2012JER_ReducedTo9NP_Plots_v2.root"));
  EL_RETURN_CHECK("initialize()",m_jerTool->setProperty("CollectionName", "AntiKt4EMTopoJets"));
  EL_RETURN_CHECK("initialize()",m_jerTool->initialize());
  // Configure the JERSmearingTool
  m_jerHandle = ToolHandle<IJERTool>(m_jerTool->name());
  m_jerSmearingTool = new JERSmearingTool("JERSmearingTool");
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("ApplyNominalSmearing", false));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("JERTool", m_jerHandle));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("isMC", !m_isData));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("SystematicMode", "Simple")); //"Simple" provides one NP (smearing only in MC), "Full" provides 10NPs (smearing both on data and MC)
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->initialize());

  // Configure the JVT tool.
  m_jvtag = 0;
  m_jvtag = new JetVertexTaggerTool("jvtag");
  m_jvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
  EL_RETURN_CHECK("initialize()",m_jvtag->setProperty("JVTFileName","JetMomentTools/JVTlikelihood_20140805.root"));
  EL_RETURN_CHECK("initialize()",m_jvtag->initialize());

  // Initialize and configure the jet cleaning tool
  m_jetCleaningTight = new JetCleaningTool("JetCleaningTight");
  m_jetCleaningLoose = new JetCleaningTool("JetCleaningLoose");
  m_jetCleaningTight->msg().setLevel( MSG::DEBUG ); 
  m_jetCleaningLoose->msg().setLevel( MSG::DEBUG ); 
  EL_RETURN_CHECK("initialize()",m_jetCleaningTight->setProperty( "CutLevel", "TightBad"));
  EL_RETURN_CHECK("initialize()",m_jetCleaningLoose->setProperty( "CutLevel", "LooseBad"));
  //EL_RETURN_CHECK("initialize()",m_jetCleaningTight->setProperty("DoUgly", false));
  //EL_RETURN_CHECK("initialize()",m_jetCleaningLoose->setProperty("DoUgly", false));
  EL_RETURN_CHECK("initialize()",m_jetCleaningTight->initialize());
  EL_RETURN_CHECK("initialize()",m_jetCleaningLoose->initialize());

  // Initialise MET tools
  m_metMaker = new met::METMaker("METMakerTool");
  EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoRemoveMuonJets", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoSetMuonJetEMScale", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoMuonEloss", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoIsolMuonEloss", true));
  EL_RETURN_CHECK("initialize()",m_metMaker->setProperty("JetMinWeightedPt", 20000.));
  EL_RETURN_CHECK("initialize()",m_metMaker->setProperty("JetMinEFrac", 0.0));
  //m_metMaker->msg().setLevel( MSG::VERBOSE ); // or DEBUG or VERBOSE
  EL_RETURN_CHECK("initialize()",m_metMaker->initialize());

  // Initialize the harmonization reccommendation tools
  const bool doTaus = true, doPhotons = true;
  const bool boostedLeptons = false;
  EL_RETURN_CHECK("initialize()",ORUtils::recommendedTools(m_toolBox, "OverlapRemovalTool", 
                                                          inputLabel, outputLabel, bJetLabel, 
                                                          boostedLeptons, outputPassValue, 
                                                          doTaus, doPhotons));
  //EL_RETURN_CHECK("initialize()",ORUtils::harmonizedTools(m_toolBox, "OverlapRemovalTool", 
  //      inputLabel, outputLabel,
  //      outputPassValue, doTaus, doPhotons));
  // Set message level for all tools
  //m_toolBox->setMsgLevel(MSG::DEBUG);
  // Initialize all tools
  m_orTool = static_cast<ORUtils::OverlapRemovalTool*>(m_toolBox.getMasterTool());
  m_orTool->setName("ORTool");
  EL_RETURN_CHECK("initialize()",m_toolBox.initialize());


  // Initialise Muon Efficiency Tool
  m_muonEfficiencySFTool = new CP::MuonEfficiencyScaleFactors( "MuonEfficiencySFTool" );
  EL_RETURN_CHECK("initialize()",m_muonEfficiencySFTool->setProperty("WorkingPoint", "Loose") );
  EL_RETURN_CHECK("initialize()",m_muonEfficiencySFTool->setProperty("CalibrationRelease", "Data15_allPeriods_260116"));
  EL_RETURN_CHECK("initialize()",m_muonEfficiencySFTool->initialize() );
  // Initialise Muon Isolation Tool
  m_muonIsolationSFTool = new CP::MuonEfficiencyScaleFactors( "MuonIsolationSFTool" );
  EL_RETURN_CHECK("initialize()",m_muonIsolationSFTool->setProperty("WorkingPoint", "LooseTrackOnlyIso") );
  EL_RETURN_CHECK("initialize()",m_muonIsolationSFTool->setProperty("CalibrationRelease", "Data15_allPeriods_260116"));
  EL_RETURN_CHECK("initialize()",m_muonIsolationSFTool->initialize() );
  // Initialise Muon TTVA Efficiency Tool
  m_muonTTVAEfficiencySFTool = new CP::MuonEfficiencyScaleFactors( "MuonTTVAEfficiencySFTool" );
  EL_RETURN_CHECK("initialize()",m_muonTTVAEfficiencySFTool->setProperty("WorkingPoint", "TTVA") );
  EL_RETURN_CHECK("initialize()",m_muonTTVAEfficiencySFTool->setProperty("CalibrationRelease", "Data15_allPeriods_260116"));
  EL_RETURN_CHECK("initialize()",m_muonTTVAEfficiencySFTool->initialize() );

  // Initialise Muon Trigger Scale Factor Tool
  m_muonTriggerSFTool = new CP::MuonTriggerScaleFactors( "MuonTriggerSFTool" );
  EL_RETURN_CHECK("initialize()",m_muonTriggerSFTool->setProperty("MuonQuality", "Loose"));
  EL_RETURN_CHECK("initialize()",m_muonTriggerSFTool->setProperty("Isolation", "LooseTrackOnly"));
  EL_RETURN_CHECK("initialize()",m_muonTriggerSFTool->initialize() );

  // Initialise Electron Efficiency Tool
  m_elecEfficiencySFTool_reco = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_reco");
  std::vector< std::string > corrFileNameList_reco;
  corrFileNameList_reco.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.RecoTrk.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("CorrectionFileNameList", corrFileNameList_reco) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->initialize() );

  m_elecEfficiencySFTool_id_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Loose");
  std::vector< std::string > corrFileNameList_id_Loose;
  corrFileNameList_id_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.LooseAndBLayerLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("CorrectionFileNameList", corrFileNameList_id_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->initialize() );

  m_elecEfficiencySFTool_id_Medium = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Medium");
  std::vector< std::string > corrFileNameList_id_Medium;
  corrFileNameList_id_Medium.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.MediumLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("CorrectionFileNameList", corrFileNameList_id_Medium) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->initialize() );

  m_elecEfficiencySFTool_id_Tight = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Tight");
  std::vector< std::string > corrFileNameList_id_Tight;
  corrFileNameList_id_Tight.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.TightLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->setProperty("CorrectionFileNameList", corrFileNameList_id_Tight) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->initialize() );

  m_elecEfficiencySFTool_iso_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Loose");
  std::vector< std::string > corrFileNameList_iso_Loose;
  corrFileNameList_iso_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("CorrectionFileNameList", corrFileNameList_iso_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->initialize() );

  m_elecEfficiencySFTool_iso_Medium = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Medium");
  std::vector< std::string > corrFileNameList_iso_Medium;
  corrFileNameList_iso_Medium.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.MediumLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("CorrectionFileNameList", corrFileNameList_iso_Medium) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->initialize() );

  m_elecEfficiencySFTool_iso_Tight = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Tight");
  std::vector< std::string > corrFileNameList_iso_Tight;
  corrFileNameList_iso_Tight.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.TightLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->setProperty("CorrectionFileNameList", corrFileNameList_iso_Tight) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->initialize() );

  m_elecEfficiencySFTool_trigEff = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_trigEff");
  std::vector< std::string > corrFileNameList_trigEff;
  corrFileNameList_trigEff.push_back("ElectronEfficiencyCorrection/efficiency.e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->setProperty("CorrectionFileNameList", corrFileNameList_trigEff) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->initialize() );

  m_elecEfficiencySFTool_trigSF_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_trigSF_Loose");
  std::vector< std::string > corrFileNameList_trigSF_Loose;
  corrFileNameList_trigSF_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->setProperty("CorrectionFileNameList", corrFileNameList_trigSF_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->initialize() );

  // Initialise Jet JVT Efficiency Tool
  m_jvtefficiencyTool = new CP::JetJvtEfficiency("JvtEfficiencyTool");
  //EL_RETURN_CHECK("initialize()",m_jvtefficiencyTool->setProperty("WorkingPoint",) );
  EL_RETURN_CHECK("initialize()",m_jvtefficiencyTool->initialize() );

  // Initialise Tau Efficiency Tool
  m_tauEffTool = new TauAnalysisTools::TauEfficiencyCorrectionsTool("TauEffTool");
  EL_RETURN_CHECK("initialize()",m_tauEffTool->initialize() );

  // Initialise MET Tools
  m_metSystTool = new met::METSystematicsTool("METSystTool");
  EL_RETURN_CHECK("initialize()",m_metSystTool->setProperty("JetColl", "AntiKt4EMTopoJets") );
  EL_RETURN_CHECK("initialize()",m_metSystTool->setProperty("ConfigSoftTrkFile", "TrackSoftTerms.config") );
  EL_RETURN_CHECK("initialize()",m_metSystTool->initialize() );

  // Initialise Isolation Correction Tool
  m_isoCorrTool = new CP::IsolationCorrectionTool( "IsoCorrTool" );
  EL_RETURN_CHECK("initialize()",m_isoCorrTool->setProperty( "IsMC", !m_isData) );
  //EL_RETURN_CHECK("initialize()",m_isoCorrTool->setProperty( "AFII_corr", m_isAtlfast) );
  EL_RETURN_CHECK("initialize()",m_isoCorrTool->initialize() );

  // Initialise PileupReweighting Tool
  m_prwTool = new CP::PileupReweightingTool("PrwTool");
  std::vector<std::string> file_conf;
  // xml file should be put in ZinvAnalysis/share directory
  file_conf.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/ZinvAnalysis/PRW.root"));
  std::vector<std::string> file_ilumi;
  file_ilumi.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/ZinvAnalysis/ilumicalc_histograms_None_276262-284484.root"));
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("ConfigFiles", file_conf) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("LumiCalcFiles", file_ilumi) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactor",     1. / 1.16) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactorUP",   1.) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactorDOWN", 1. / 1.23) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("UnrepresentedDataAction", 2));
  if ( !(mcChannelNumber == 363121 || mcChannelNumber == 363351) ) {
    EL_RETURN_CHECK("initialize()",m_prwTool->initialize() );
  }    


  // Initialize PMGTools (PMGSherpa22VJetsWeightTool)
  // See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/CentralMC15ProductionList#NEW_Sherpa_v2_2_V_jets_NJet_rewe 
  // Function of njettruth which is the number of truth jets with 
  // pt>20 and |eta|<4.5
  m_PMGSherpa22VJetsWeightTool = new PMGSherpa22VJetsWeightTool("PMGSherpa22VJetsWeightTool");
  EL_RETURN_CHECK("initialize()",m_PMGSherpa22VJetsWeightTool->setProperty("TruthJetContainer","AntiKt4TruthJets") );
  EL_RETURN_CHECK("initialize()",m_PMGSherpa22VJetsWeightTool->setProperty("TruthParticleContainer","TruthParticles") );
  EL_RETURN_CHECK("initialize()",m_PMGSherpa22VJetsWeightTool->initialize() );


  // Get the systematics registry and add the recommended systematics into our list of systematics to run over (+/-1 sigma):
  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  const CP::SystematicSet& recommendedSystematics = registry.recommendedSystematics(); // get list of recommended systematics
  m_sysList = CP::make_systematics_vector(recommendedSystematics); 


  // Initialize Cutflow
  if (m_useBitsetCutflow)
    m_BitsetCutflow = new BitsetCutflow(wk());


  // Initialize Cutflow count array
  for (int i=0; i<40; i++) {
    m_eventCutflow[i]=0;
  }


  ///////////////////////////////
  // Create Histograms //////////
  ///////////////////////////////

  // Publication bins
  // Monojet and VBF MET
  Float_t binsMET[] = {200.,250.,300.,350.,500.,700.,1000.,1400.};
  Int_t nbinMET = sizeof(binsMET)/sizeof(Float_t) - 1;
  // VBF mjj
  Float_t binsMjj[] = {200.,400.,600.,1000.,2000.,3000.,4000.};
  Int_t nbinMjj = sizeof(binsMjj)/sizeof(Float_t) - 1;
  // VBF dPhi(j1,j2)
  Float_t pi = TMath::Pi();
  Float_t binsDPhi[] = {0., pi/(7-1), 2*pi/(7-1), 3*pi/(7-1), 4*pi/(7-1), 5*pi/(7-1), pi}; // where 7 is number of bins
  Int_t nbinDPhi = sizeof(binsDPhi)/sizeof(Float_t) - 1;


  TH1::SetDefaultSumw2(kTRUE);

  for (const auto &sysList : m_sysList){
    if ((!m_doSys || m_isData) && (sysList).name() != "") continue;
    std::string sysName = (sysList).name();


    if (m_doSys && (sysName.find("TAUS_")!=std::string::npos || sysName.find("PH_")!=std::string::npos )) continue;
    if (m_isZmumu && m_doSys && ((sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos))) continue;
    if (m_isZee && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos))) continue;
    if (m_isZnunu && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos || sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos)) ) continue;


    if (m_isZmumu) {
      h_channel = "h_zmumu_";

      // Publication plots
      addHist(hMap1D, "ZmumuMET_mono"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "ZmumuMET_search"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "ZmumuMjj_search"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, "ZmumuDeltaPhiAll"+sysName, nbinDPhi, binsDPhi);

      ////////////////////////
      // Monojet phasespace //
      ////////////////////////
      // MET
      addHist(hMap1D, h_channel+"monojet_met"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"monojet_met_emulmet"+sysName, 30, 0., 1500.);
      // Jets
      addHist(hMap1D, h_channel+"monojet_njet"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"monojet_jet_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"monojet_jet_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"monojet_jet_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_jet_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_dPhimetjet"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"monojet_dPhiMinmetjet"+sysName, 16, 0., 3.2);
      // Leptons
      addHist(hMap1D, h_channel+"monojet_lepton1_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"monojet_lepton2_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"monojet_lepton1_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"monojet_lepton2_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"monojet_lepton1_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_lepton2_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_mll"+sysName, 150, 0., 300.);

      ////////////////////
      // VBF phasespace //
      ////////////////////
      // MET
      addHist(hMap1D, h_channel+"vbf_met"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"vbf_met_emulmet"+sysName, 30, 0., 1500.);
      // Jets
      addHist(hMap1D, h_channel+"vbf_mjj"+sysName, 80, 0., 4000.);
      addHist(hMap1D, h_channel+"vbf_dPhijj"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_njet"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"vbf_jet1_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"vbf_jet2_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"vbf_jet3_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"vbf_jet1_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_jet2_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_jet3_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_jet1_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet2_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet3_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet1_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet2_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet3_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_dRjj"+sysName, 25, 0., 5.);
      addHist(hMap1D, h_channel+"vbf_dPhimetj1"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_dPhimetj2"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_dPhimetj3"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_dPhiMinmetjet"+sysName, 16, 0., 3.2);
      // Leptons
      addHist(hMap1D, h_channel+"vbf_lepton1_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"vbf_lepton2_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"vbf_lepton1_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_lepton2_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_lepton1_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_lepton2_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_mll"+sysName, 150, 0., 300.);
    }

    if (m_isZee) {
      h_channel = "h_zee_";

      // Publication plots
      addHist(hMap1D, "ZeeMET_mono"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "ZeeMET_search"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "ZeeMjj_search"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, "ZeeDeltaPhiAll"+sysName, nbinDPhi, binsDPhi);

      ////////////////////////
      // Monojet phasespace //
      ////////////////////////
      // MET
      addHist(hMap1D, h_channel+"monojet_met"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"monojet_met_emulmet"+sysName, 30, 0., 1500.);
      // Jets
      addHist(hMap1D, h_channel+"monojet_njet"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"monojet_jet_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"monojet_jet_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"monojet_jet_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_jet_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_dPhimetjet"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"monojet_dPhiMinmetjet"+sysName, 16, 0., 3.2);
      // Leptons
      addHist(hMap1D, h_channel+"monojet_lepton1_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"monojet_lepton2_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"monojet_lepton1_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"monojet_lepton2_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"monojet_lepton1_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_lepton2_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"monojet_mll"+sysName, 150, 0., 300.);

      ////////////////////
      // VBF phasespace //
      ////////////////////
      // MET
      addHist(hMap1D, h_channel+"vbf_met"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"vbf_met_emulmet"+sysName, 30, 0., 1500.);

      // Jets
      addHist(hMap1D, h_channel+"vbf_mjj"+sysName, 80, 0., 4000.);
      addHist(hMap1D, h_channel+"vbf_dPhijj"+sysName, 16, 0., 3.2);

      addHist(hMap1D, h_channel+"vbf_njet"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"vbf_jet1_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"vbf_jet2_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"vbf_jet3_pt"+sysName, 60, 0., 3000.);
      addHist(hMap1D, h_channel+"vbf_jet1_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_jet2_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_jet3_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_jet1_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet2_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet3_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet1_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet2_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_jet3_rap"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_dRjj"+sysName, 25, 0., 5.);
      addHist(hMap1D, h_channel+"vbf_dPhimetj1"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_dPhimetj2"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_dPhimetj3"+sysName, 16, 0., 3.2);
      addHist(hMap1D, h_channel+"vbf_dPhiMinmetjet"+sysName, 16, 0., 3.2);

      // Leptons
      addHist(hMap1D, h_channel+"vbf_lepton1_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"vbf_lepton2_pt"+sysName, 30, 0., 1500.);
      addHist(hMap1D, h_channel+"vbf_lepton1_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_lepton2_phi"+sysName, 32, -3.2, 3.2);
      addHist(hMap1D, h_channel+"vbf_lepton1_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_lepton2_eta"+sysName, 25, -5., 5.);
      addHist(hMap1D, h_channel+"vbf_mll"+sysName, 150, 0., 300.);
    }




  } // Systematic loop end




  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.


  // push cutflow bitset to cutflow hist
  if (m_useBitsetCutflow)
    m_BitsetCutflow->PushBitSet();

  // print every 100 events, so we know where we are:
  if( (m_eventCounter % 100) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;
  if (m_useArrayCutflow) m_eventCutflow[0]+=1;

  //----------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection in execute. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Calculate EventWeight
  mcEventWeight = 1.;
  mcEventWeight_Zmumu = 1.;
  mcEventWeight_Wmunu = 1.;
  mcEventWeight_Zee = 1.;
  mcEventWeight_Wenu = 1.;

  float mcWeight = 1.;
  if (!m_isData) {
    mcWeight = eventInfo->mcEventWeight();
  }
  if (((mcChannelNumber > 363100 && mcChannelNumber < 363500) || (mcChannelNumber > 304010 && mcChannelNumber < 304025)) && std::abs(mcWeight) > 100.) mcWeight = 1.;

  // Sherpa v2.2 V+jets NJet reweight 
  if (!m_isData) {
    // See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/CentralMC15ProductionList#NEW_Sherpa_v2_2_V_jets_NJet_rewe
    // only reweight W and Z samples
    if (mcChannelNumber > 363100 &&  mcChannelNumber < 363500) {
      double sherpaReweightValue = m_PMGSherpa22VJetsWeightTool->getWeight(); 
      mcWeight*=sherpaReweightValue;
    }
  }




  /*
  // BCID Information
  if (m_isData){
    int Bcid = eventInfo->bcid();
    if (Bcid > 324 && Bcid < 417){
      Info("execute()", "  BCID = %d", Bcid);
    }
  }
  */


  // if data check if event passes GRL
  if(m_isData){ // it's data!
    if(!m_grl->passRunLB(*eventInfo)){
      return EL::StatusCode::SUCCESS; // go to next event
    }
  } // end if not MC
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("GRL");
  if (m_useArrayCutflow) m_eventCutflow[1]+=1;


  //------------------------------------------------------------
  // Apply event cleaning to remove events due to 
  // problematic regions of the detector, and incomplete events.
  // Apply to data.
  //------------------------------------------------------------
  // reject event if:
  if(m_isData){
    if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) 
        || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) 
        || (eventInfo->errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error) 
        || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) )  )
    {
      return EL::StatusCode::SUCCESS; // go to the next event
    } // end if event flags check
  } // end if the event is data
  m_numCleanEvents++;
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("LAr_Tile_Core");
  if (m_useArrayCutflow) m_eventCutflow[2]+=1;


  //---------------------
  // EVENT SELECTION
  //---------------------

  //---------------------
  // Retrive vertex object and select events with at least one good primary vertex with at least 2 tracks
  //---------------------
  const xAOD::VertexContainer* vertices(0);
  /// retrieve arguments: container type, container key
  if ( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){ 
    Error("execute()","Failed to retrieve PrimaryVertices container. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  xAOD::Vertex* primVertex = 0;
  int nGoodVtx = 0;
  for( ; vtx_itr != vtx_end; ++vtx_itr ) {
    if ((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx){
      primVertex = (*vtx_itr);
      nGoodVtx++;
    }
  }


  //--------------
  // Preseletion
  //--------------
  //if (nGoodVtx==1)
  //   Info("execute()", "  Found one prim.vertex: nGoodVtx = %d", nGoodVtx); // just to print out something
  if (nGoodVtx>1)
    Info("execute()", "  WARNING!!!! Found more than one prim.vertex: nGoodVtx = %d", nGoodVtx); // just to print out something
  if (nGoodVtx==0)
    //   Info("execute()", "  %s", "No one prim.vertex found"); // just to print out something
    return EL::StatusCode::SUCCESS;
  if (primVertex->nTrackParticles() < 2) return EL::StatusCode::SUCCESS;

  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Primary vertex");
  if (m_useArrayCutflow) m_eventCutflow[3]+=1;




  //------------
  // MUONS
  //------------

  /// full copy 
  // get muon container of interest
  const xAOD::MuonContainer* m_muons(0);
  if ( !m_event->retrieve( m_muons, "Muons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //------------
  // ELECTRONS
  //------------

  /// full copy 
  // get electron container of interest
  const xAOD::ElectronContainer* m_electrons(0);
  if ( !m_event->retrieve( m_electrons, "Electrons" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Electron container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //------------
  // PHOTONS
  //------------

  /// full copy 
  // get photon container of interest
  const xAOD::PhotonContainer* m_photons(0);
  if ( !m_event->retrieve( m_photons, "Photons" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Photon container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //------------
  // TAUS
  //------------

  /// full copy 
  // get tau container of interest
  const xAOD::TauJetContainer* m_taus(0);
  if ( !m_event->retrieve( m_taus, "TauJets" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Tau container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //------------
  // JETS
  //------------

  /// full copy 
  // get jet container of interest
  const xAOD::JetContainer* m_jets(0);
  if ( !m_event->retrieve( m_jets, jetType ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Jet container. Exiting." );
    return EL::StatusCode::FAILURE;
  }


  //-----------------------
  // Systematics Start
  //-----------------------
  // loop over recommended systematics
  for (const auto &sysList : m_sysList){
    std::string sysName = (sysList).name();
    if ((!m_doSys || m_isData) && sysName != "") continue;

    if (m_doSys && (sysName.find("TAUS_")!=std::string::npos || sysName.find("PH_")!=std::string::npos )) continue;
    if (m_isZmumu && m_doSys && ((sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos))) continue;
    if (m_isZee && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos))) continue;
    if (m_isZnunu && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos || sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos)) ) continue;

    //if(sysName=="") std::cout << "Nominal (no syst) "  << std::endl;
    //else std::cout << "Systematic: " << sysName << std::endl;


    if (!m_isData) {

      if (m_jerSmearingTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure JERSmearingTool for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_jetUncertaintiesTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure JetUncertaintiesTool for systematics");
        return EL::StatusCode::FAILURE; 
      }
      if (m_muonEfficiencySFTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure MuonEfficiencyScaleFactorsToolSF for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_muonIsolationSFTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure MuonIsolationEfficiencyScaleFactorsToolSF for systematics");
        return EL::StatusCode::FAILURE; 
      } 
      if (m_muonTTVAEfficiencySFTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure MuonTTVAEfficiencyScaleFactorsToolSF for systematics");
        return EL::StatusCode::FAILURE; 
      } 
      if (m_muonCalibrationAndSmearingTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure MuonCalibrationAndSmearingTool for systematics");
        return EL::StatusCode::FAILURE;
      } 
      if (m_elecEfficiencySFTool_reco->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure electronEfficiencyCorrectionToolRecoSF for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_elecEfficiencySFTool_id_Loose->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure electronEfficiencyCorrectionToolIdLooseSF for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_elecEfficiencySFTool_id_Tight->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure electronEfficiencyCorrectionToolIdTightSF for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_elecEfficiencySFTool_iso_Loose->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure electronEfficiencyCorrectionToolIsoSFlooseID for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_elecEfficiencySFTool_iso_Tight->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure electronEfficiencyCorrectionToolIsoSFtightID for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_elecEfficiencySFTool_trigSF_Loose->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure electronEfficiencyCorrectionToolTriggerSFloose for systematics");
        return EL::StatusCode::FAILURE;
      }
      if (m_egammaCalibrationAndSmearingTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure EgammaCalibrationAndSmearingTool for systematics");
        return EL::StatusCode::FAILURE;
      } 
      if (m_isoCorrTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure IsolationCorrectionTool for systematics");
        return EL::StatusCode::FAILURE;
      } 
      if (m_tauEffTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure TauEfficiencyCorrectionsTool for systematics");
        return EL::StatusCode::FAILURE;
      } 
      if (m_tauSmearingTool->applySystematicVariation(sysList) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure TauSmearingTool for systematics");
        return EL::StatusCode::FAILURE;
      } 
      if (m_metSystTool->applySystematicVariation( sysList ) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure METSystematicsTool for systematics");
        return EL::StatusCode::FAILURE;
      }   
      if (m_prwTool->applySystematicVariation( sysList ) != CP::SystematicCode::Ok) {
        Error("execute()", "Cannot configure PileupReweightingTool for systematics");
        return EL::StatusCode::FAILURE;
      }

    }





    //---------------------
    // Pile-up reweighting
    //---------------------
    if (!m_isData) {
      if (mcChannelNumber == 363121 || mcChannelNumber == 363351) { // These samples have missing mu values and the pileup reweighting tool doesn't like that and crashes.
        mcEventWeight = mcWeight;
      }
      else {
        float pu_weight = m_prwTool->getCombinedWeight(*eventInfo); // Get Pile-up weight
        mcEventWeight = mcWeight * pu_weight;
      }
    }




    ///////////////////
    // For VBF study //
    ///////////////////

    //------------
    // MUONS
    //------------
    /// shallow copy for muon calibration and smearing tool
    // create a shallow copy of the muons container for MET building
    std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > muons_shallowCopy = xAOD::shallowCopyContainer( *m_muons );
    xAOD::MuonContainer* muonSC = muons_shallowCopy.first;

    // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
    // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
    // The method is defined in the header file xAODBase/IParticleHelpers.h
    bool setLinksMuon = xAOD::setOriginalObjectLink(*m_muons,*muonSC);
    if(!setLinksMuon) {
      Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
      return StatusCode::FAILURE;
    }
    // iterate over our shallow copy
    for (const auto& muon : *muonSC) { // C++11 shortcut
      // VBF Muon Selection
      passMuonVBF(*muon, eventInfo, primVertex);
      //Info("execute()", "  VBF muon pt = %.2f GeV", (muon->pt() * 0.001));
    } // end for loop over shallow copied muons


    //------------
    // ELECTRONS
    //------------
    /// shallow copy for electron calibration tool
    // create a shallow copy of the electrons container for MET building
    std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* > elec_shallowCopy = xAOD::shallowCopyContainer( *m_electrons );
    xAOD::ElectronContainer* elecSC = elec_shallowCopy.first;

    // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
    // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
    // The method is defined in the header file xAODBase/IParticleHelpers.h
    bool setLinksElec = xAOD::setOriginalObjectLink(*m_electrons,*elecSC);
    if(!setLinksElec) {
      Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
      return StatusCode::FAILURE;
    }
    // iterate over our shallow copy
    for (const auto& electron : *elecSC) { // C++11 shortcut
      // VBF Electron Selection
      passElectronVBF(*electron, eventInfo, primVertex);
      //Info("execute()", "  VBF electron pt = %.2f GeV", (electron->pt() * 0.001));
    } // end for loop over shallow copied electrons





    //------------
    // PHOTONS
    //------------
    /// shallow copy for photon calibration tool
    // create a shallow copy of the photons container for MET building
    std::pair< xAOD::PhotonContainer*, xAOD::ShallowAuxContainer* > phot_shallowCopy = xAOD::shallowCopyContainer( *m_photons );
    xAOD::PhotonContainer* photSC = phot_shallowCopy.first;

    // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
    // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
    // The method is defined in the header file xAODBase/IParticleHelpers.h
    bool setLinksPhoton = xAOD::setOriginalObjectLink(*m_photons,*photSC);
    if(!setLinksPhoton) {
      Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
      return StatusCode::FAILURE;
    }
    // iterate over our shallow copy
    for (const auto& photon : *photSC) { // C++11 shortcut
      // VBF Tau Selection
      passPhotonVBF(*photon, eventInfo); 
    } // end for loop over shallow copied photons



    //------------
    // TAUS
    //------------
    m_tauOverlappingElectronLLHDecorator->initializeEvent();
    /// shallow copy for tau calibration tool
    // create a shallow copy of the taus container for MET building
    std::pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* > tau_shallowCopy = xAOD::shallowCopyContainer( *m_taus );
    xAOD::TauJetContainer* tauSC = tau_shallowCopy.first;

    // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
    // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
    // The method is defined in the header file xAODBase/IParticleHelpers.h
    bool setLinksTau = xAOD::setOriginalObjectLink(*m_taus,*tauSC);
    if(!setLinksTau) {
      Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
      return StatusCode::FAILURE;
    }
    // iterate over our shallow copy
    for (const auto& taujet : *tauSC) { // C++11 shortcut
      // VBF Tau Selection
      passTauVBF(*taujet, eventInfo);
    } // end for loop over shallow copied taus



    //------------
    // JETS
    //------------
    /// shallow copy for jet calibration tool
    // create a shallow copy of the jets container for MET building
    std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > jet_shallowCopy = xAOD::shallowCopyContainer( *m_jets );
    xAOD::JetContainer* jetSC = jet_shallowCopy.first;

    // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
    // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
    // The method is defined in the header file xAODBase/IParticleHelpers.h
    bool setLinksJet = xAOD::setOriginalObjectLink(*m_jets,*jetSC);
    if(!setLinksJet) {
      Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
      return StatusCode::FAILURE;
    }

    // iterate over our shallow copy
    for (const auto& jets : *jetSC) { // C++11 shortcut
      //Info("execute()", "  original jet pt = %.2f GeV", jets->pt() * 0.001);

      // According to https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15

      // JES calibration
      if ( !m_jetCalibration->applyCalibration(*jets).isSuccess() ){
        Error("execute()", "Failed to apply calibration to Jet objects. Exiting." );
        return EL::StatusCode::FAILURE;
      }

      float jetPt = jets->pt() * 0.001; /// GeV
      float jetEta = jets->eta();

      // JES correction
      if (!m_isData){
        if ( m_jetUncertaintiesTool->applyCorrection(*jets) != CP::CorrectionCode::Ok){ // apply correction and check return code
          Error("execute()", "Failed to apply JES correction to Jet objects. Exiting." );
          return EL::StatusCode::FAILURE;
        }
      }

      // JER smearing
      if (!m_isData){
        if ( m_jerSmearingTool->applyCorrection(*jets) != CP::CorrectionCode::Ok){ // apply correction and check return code
          Error("execute()", "Failed to apply JER smearing. Exiting. ");
          return EL::StatusCode::FAILURE;
        }
      }

      // JVT Tool
      float newjvt = m_jvtagup->updateJvt(*jets);
      acc_jvt(*jets) = newjvt;

      //Info("execute()", "  corrected jet pt = %.2f GeV", jets->pt() * 0.001);
      //Info("execute()", "  updated jet jvt = %.2f ", newjvt);

      dec_baseline(*jets) = false;
      selectDec(*jets) = false; // To select objects for Overlap removal

      // pT cut
      if (jetPt > m_jetPtCut && fabs(jetEta) < m_jetEtaCut) {
        dec_baseline(*jets) = true;
        selectDec(*jets) = true; // To select objects for Overlap removal
      }

    } // end for loop over shallow copied jets


    //----------------
    // Overlap Removal
    //----------------

    if (m_doORtool){
      //auto m_orTool = m_toolBox->getMasterHandle();
      if ( !m_orTool->removeOverlaps(elecSC, muonSC, jetSC, tauSC, photSC).isSuccess() ){
        Error("execute()", "Failed to apply the overlap removal to all objects. Exiting." );
        return EL::StatusCode::FAILURE;
      }


      // Now, dump all of the results

      // electrons
      for(auto electron : *elecSC){
        if(overlapAcc(*electron)) {
          nOverlapElectrons++;
          //Info("execute()", "  EventNumber : %i |  Overlap electron pt = %.2f GeV", EventNumber, (electron->pt() * 0.001));
        }
        nInputElectrons++;
      }
      // muons
      for(auto muon : *muonSC){
        if(overlapAcc(*muon)) nOverlapMuons++;
        nInputMuons++;
      }
      // jets
      for (auto jet : *jetSC) {
        if(overlapAcc(*jet)){
          nOverlapJets++;
          //Info("execute()", "  EventNumber : %i |  Overlap jet pt = %.2f GeV", EventNumber, (jet->pt() * 0.001));
        }
        nInputJets++;
      }
      // taus
      for(auto tau : *tauSC){
        if(overlapAcc(*tau)) nOverlapTaus++;
        nInputTaus++;
      }
      // photons
      for(auto photon : *photSC){
        if(overlapAcc(*photon)) nOverlapPhotons++;
        nInputPhotons++;
      }
    }


    // --------------------------------------------------------
    // Select Signal Jet and Bad Jet applying overlap removal
    // --------------------------------------------------------
    /// Creating New Hard Object Containers
    // [For jet identification] filter the Jet container m_jets, placing selected jets into m_goodJet
    xAOD::JetContainer* m_goodJet = new xAOD::JetContainer(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Jet>

    bool isBadJet = false;

    // iterate over our shallow copy
    for (const auto& jets : *jetSC) { // C++11 shortcut

      // Veto Jet (cleaning Jet)
      if (IsBadJet(*jets)) isBadJet = true;

      // Jet Signal Selection
      if (IsSignalJet(*jets)) {
        //double jetPt = (jets->pt()) * 0.001; /// GeV

        m_goodJet->push_back( jets );
      }
    } // end for loop over shallow copied jets



    //------------------------------------
    // Event Cleaning (Jet cleaning tool)
    //------------------------------------
    //if (isBadJet) return EL::StatusCode::SUCCESS;
    if (isBadJet){

      delete muons_shallowCopy.first;
      delete muons_shallowCopy.second;

      delete elec_shallowCopy.first;
      delete elec_shallowCopy.second;

      delete phot_shallowCopy.first;
      delete phot_shallowCopy.second;

      delete tau_shallowCopy.first;
      delete tau_shallowCopy.second;

      delete jet_shallowCopy.first;
      delete jet_shallowCopy.second;

      delete m_goodJet;

      continue;
    }
    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Jet Cleaning");
    if (m_useArrayCutflow) m_eventCutflow[4]+=1;



    //-----------------
    // Rebuild the MET
    //-----------------

    // Create a MissingETContainer with its aux store for each systematic
    xAOD::MissingETContainer* m_met = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* m_metAux = new xAOD::MissingETAuxContainer();
    m_met->setStore(m_metAux);

    //retrieve the original containers
    const xAOD::MissingETContainer* m_metCore(0);
    std::string coreMetKey = "MET_Core_" + jetType;
    coreMetKey.erase(coreMetKey.length() - 4); //this removes the Jets from the end of the jetType
    if ( !m_event->retrieve( m_metCore, coreMetKey ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Unable to retrieve MET core container: " );
      return EL::StatusCode::FAILURE;
    }

    //retrieve the MET association map
    const xAOD::MissingETAssociationMap* m_metMap(0);
    std::string metAssocKey = "METAssoc_" + jetType;
    metAssocKey.erase(metAssocKey.length() - 4 );//this removes the Jets from the end of the jetType
    if ( !m_event->retrieve( m_metMap, metAssocKey ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Unable to retrieve MissingETAssociationMap: " );
      return EL::StatusCode::FAILURE;
    }

    // do CST or TST
    //std::string softTerm = "SoftClus";
    std::string softTerm = "PVSoftTrk";

    // It is necessary to reset the selected objects before every MET calculation
    m_metMap->resetObjSelectionFlags();



    //==========================================================
    // For rebuild the emulated MET (by marking Muon invisible)
    //==========================================================


    // Create a MissingETContainer with its aux store for each systematic
    xAOD::MissingETContainer* m_emulmet_nomu = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* m_emulmetAux_nomu = new xAOD::MissingETAuxContainer();
    m_emulmet_nomu->setStore(m_emulmetAux_nomu);

    //retrieve the original containers
    const xAOD::MissingETContainer* m_emulmetCore_nomu(0);
    if ( !m_event->retrieve( m_emulmetCore_nomu, coreMetKey ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Unable to retrieve MET core container: " );
      return EL::StatusCode::FAILURE;
    }

    //retrieve the MET association map
    const xAOD::MissingETAssociationMap* m_emulmetMap_nomu(0);
    if ( !m_event->retrieve( m_emulmetMap_nomu, metAssocKey ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Unable to retrieve MissingETAssociationMap: " );
      return EL::StatusCode::FAILURE;
    }

    // It is necessary to reset the selected objects before every MET calculation
    m_emulmetMap_nomu->resetObjSelectionFlags();


    //==============================================================
    // For rebuild the emulated MET (by marking Electron invisible)
    //==============================================================


    // Create a MissingETContainer with its aux store for each systematic
    xAOD::MissingETContainer* m_emulmet_noelec = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* m_emulmetAux_noelec = new xAOD::MissingETAuxContainer();
    m_emulmet_noelec->setStore(m_emulmetAux_noelec);

    //retrieve the original containers
    const xAOD::MissingETContainer* m_emulmetCore_noelec(0);
    if ( !m_event->retrieve( m_emulmetCore_noelec, coreMetKey ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Unable to retrieve MET core container: " );
      return EL::StatusCode::FAILURE;
    }

    //retrieve the MET association map
    const xAOD::MissingETAssociationMap* m_emulmetMap_noelec(0);
    if ( !m_event->retrieve( m_emulmetMap_noelec, metAssocKey ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Unable to retrieve MissingETAssociationMap: " );
      return EL::StatusCode::FAILURE;
    }

    // It is necessary to reset the selected objects before every MET calculation
    m_emulmetMap_noelec->resetObjSelectionFlags();





    //here we apply some basic cuts and rebuild the met at each step

    // Electron
    //-----------------

    /// Creat New Hard Object Containers
    // [For MET building] filter the Electron container m_electrons, placing selected electrons into m_MetElectrons
    ConstDataVector<xAOD::ElectronContainer> m_MetElectrons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Electron>

    // iterate over our shallow copy
    for (const auto& electron : *elecSC) { // C++11 shortcut
      // For MET rebuilding
      if (dec_baseline(*electron)) {
        m_MetElectrons.push_back( electron );
      }
    } // end for loop over shallow copied electrons
    //const xAOD::ElectronContainer* p_MetElectrons = m_MetElectrons.asDataVector();
    // For real MET
    m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
        xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
        m_met,                                      //filling this met container
        m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
        m_metMap);                                  //and this association map

    // For emulated MET marking electrons invisible
    m_metMaker->markInvisible(m_MetElectrons.asDataVector(), m_emulmetMap_noelec);

    // For emulated MET marking muons invisible
    m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
        xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
        m_emulmet_nomu,                             //filling this met container
        m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
        m_emulmetMap_nomu);                         //and this association map


    // Photon
    //-----------------

    /// Creat New Hard Object Containers
    // [For MET building] filter the Photon container m_photons, placing selected photons into m_MetPhotons
    ConstDataVector<xAOD::PhotonContainer> m_MetPhotons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Photon>

    // iterate over our shallow copy
    for (const auto& photon : *photSC) { // C++11 shortcut
      // For MET rebuilding
      if (dec_baseline(*photon)) {
        m_MetPhotons.push_back( photon );
      }
    } // end for loop over shallow copied photons
    // For real MET
    m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
        xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
        m_met,                                    //filling this met container
        m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
        m_metMap);                                //and this association map

    // For emulated MET marking electrons invisible
    m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
        xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
        m_emulmet_noelec,                         //filling this met container
        m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
        m_emulmetMap_noelec);                     //and this association map

    // For emulated MET marking muons invisible
    m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
        xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
        m_emulmet_nomu,                           //filling this met container
        m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
        m_emulmetMap_nomu);                       //and this association map


    // TAUS
    //-----------------

    /// Creat New Hard Object Containers
    // [For MET building] filter the TauJet container m_taus, placing selected taus into m_MetTaus
    ConstDataVector<xAOD::TauJetContainer> m_MetTaus(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::TauJet>

    // iterate over our shallow copy
    for (const auto& taujet : *tauSC) { // C++11 shortcut
      // For MET rebuilding
      if (dec_baseline(*taujet)) {
        m_MetTaus.push_back( taujet );
      }
    } // end for loop over shallow copied taus
    // For real MET
    m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
        xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
        m_met,                                 //filling this met container
        m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
        m_metMap);                             //and this association map

    // For emulated MET marking electrons invisible
    m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
        xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
        m_emulmet_noelec,                      //filling this met container
        m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
        m_emulmetMap_noelec);                  //and this association map

    // For emulated MET marking muons invisible
    m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
        xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
        m_emulmet_nomu,                        //filling this met container
        m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
        m_emulmetMap_nomu);                    //and this association map


    // Muon
    //-----------------

    /// Creat New Hard Object Containers
    // [For MET building] filter the Muon container m_muons, placing selected muons into m_MetMuons
    ConstDataVector<xAOD::MuonContainer> m_MetMuons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Muon>

    // iterate over our shallow copy
    for (const auto& muon : *muonSC) { // C++11 shortcut
      // For MET rebuilding
      if (dec_baseline(*muon)) {
        m_MetMuons.push_back( muon );
      }
    } // end for loop over shallow copied muons
    // For real MET
    m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
        xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
        m_met,                                  //filling this met container
        m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
        m_metMap);                              //and this association map

    // For emulated MET marking electrons invisible
    m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
        xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
        m_emulmet_noelec,                       //filling this met container
        m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
        m_emulmetMap_noelec);                   //and this association map

    // For emulated MET marking muons invisible
    m_metMaker->markInvisible(m_MetMuons.asDataVector(), m_emulmetMap_nomu);

    met::addGhostMuonsToJets(*m_muons, *jetSC);

    // JET
    //-----------------

    //Now time to rebuild jetMet and get the soft term
    //This adds the necessary soft term for both CST and TST
    //these functions create an xAODMissingET object with the given names inside the container
    // For real MET
    m_metMaker->rebuildJetMET("RefJet",          //name of jet met
        "SoftClus",        //name of soft cluster term met
        "PVSoftTrk",       //name of soft track term met
        m_met,             //adding to this new met container
        jetSC,             //using this jet collection to calculate jet met
        m_metCore,         //core met container
        m_metMap,          //with this association map
        true);             //apply jet jvt cut

    // For emulated MET marking electrons invisible
    m_metMaker->rebuildJetMET("RefJet",          //name of jet met
        "SoftClus",           //name of soft cluster term met
        "PVSoftTrk",          //name of soft track term met
        m_emulmet_noelec,     //adding to this new met container
        jetSC,                //using this jet collection to calculate jet met
        m_emulmetCore_noelec, //core met container
        m_emulmetMap_noelec,  //with this association map
        true);                //apply jet jvt cut

    // For emulated MET marking muons invisible
    m_metMaker->rebuildJetMET("RefJet",          //name of jet met
        "SoftClus",           //name of soft cluster term met
        "PVSoftTrk",          //name of soft track term met
        m_emulmet_nomu,       //adding to this new met container
        jetSC,                //using this jet collection to calculate jet met
        m_emulmetCore_nomu,   //core met container
        m_emulmetMap_nomu,    //with this association map
        true);                //apply jet jvt cut



    // Soft term uncertainties
    //-------------------------
    if (!m_isData) {
      // Get the track soft term (For real MET)
      xAOD::MissingET* softTrkmet = (*m_met)[softTerm];
      if (m_metSystTool->applyCorrection(*softTrkmet) != CP::CorrectionCode::Ok) {
        Error("execute()", "METSystematicsTool returns Error CorrectionCode");
      }
      // Get the track soft term (For emulated MET marking electrons invisible)
      xAOD::MissingET* softTrkmet_noelec = (*m_emulmet_noelec)[softTerm];
      if (m_metSystTool->applyCorrection(*softTrkmet_noelec) != CP::CorrectionCode::Ok) {
        Error("execute()", "METSystematicsTool returns Error CorrectionCode");
      }
      // Get the track soft term (For emulated MET marking muons invisible)
      xAOD::MissingET* softTrkmet_nomu = (*m_emulmet_nomu)[softTerm];
      if (m_metSystTool->applyCorrection(*softTrkmet_nomu) != CP::CorrectionCode::Ok) {
        Error("execute()", "METSystematicsTool returns Error CorrectionCode");
      }
    }



    // MET Build
    //-----------------

    //m_metMaker->rebuildTrackMET("RefJetTrk", softTerm, m_met, jetSC, m_metCore, m_metMap, true);

    //this builds the final track or cluster met sums, using systematic varied container
    //In the future, you will be able to run both of these on the same container to easily output CST and TST
    // For real MET
    m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());

    // For emulated MET marking electrons invisible
    m_metMaker->buildMETSum("Final", m_emulmet_noelec, (*m_emulmet_noelec)[softTerm]->source());

    // For emulated MET marking muons invisible
    m_metMaker->buildMETSum("Final", m_emulmet_nomu, (*m_emulmet_nomu)[softTerm]->source());



    // Fill real MET
    float MET_ex = -9e9;
    float MET_ey = -9e9;
    float MET = -9e9;
    float SumET = -9e9;
    float MET_phi = -9e9;

    MET_ex = ((*m_met)["Final"]->mpx()) * 0.001;
    MET_ey = ((*m_met)["Final"]->mpy()) * 0.001;
    //MET = sqrt(MET_ex*MET_ex+MET_ey*MET_ey);
    MET = ((*m_met)["Final"]->met()) * 0.001;
    SumET = ((*m_met)["Final"]->sumet()) * 0.001;
    //MET_phi = atan2(MET_ey, MET_ex);
    MET_phi = ((*m_met)["Final"]->phi());

    /*
    h_met_ex->Fill( MET_ex ); // GeV
    h_met_ey->Fill( MET_ey ); // GeV
    h_met->Fill( MET ); // GeV
    h_sumet->Fill( SumET ); // GeV
    h_met_phi->Fill( MET_phi ); // GeV
    */


    // Fill emulated MET (by marking electrons invisible)
    float emulMET_noelec_ex = -9e9;
    float emulMET_noelec_ey = -9e9;
    float emulMET_noelec = -9e9;
    float emulSumET_noelec = -9e9;
    float emulMET_noelec_phi = -9e9;

    emulMET_noelec_ex = ((*m_emulmet_noelec)["Final"]->mpx()) * 0.001;
    emulMET_noelec_ey = ((*m_emulmet_noelec)["Final"]->mpy()) * 0.001;
    emulMET_noelec = ((*m_emulmet_noelec)["Final"]->met()) * 0.001;
    emulSumET_noelec = ((*m_emulmet_noelec)["Final"]->sumet()) * 0.001;
    emulMET_noelec_phi = ((*m_emulmet_noelec)["Final"]->phi());

    /*
    h_emulmet_noelec_ex->Fill( emulMET_noelec_ex ); // GeV
    h_emulmet_noelec_ey->Fill( emulMET_noelec_ey ); // GeV
    h_emulmet_noelec->Fill( emulMET_noelec ); // GeV
    h_emulsumet_noelec->Fill( emulSumET_noelec ); // GeV
    h_emulmet_noelec_phi->Fill( emulMET_noelec_phi ); // GeV
    */


    // Fill emulated MET (by marking muons invisible)
    float emulMET_nomu_ex = -9e9;
    float emulMET_nomu_ey = -9e9;
    float emulMET_nomu = -9e9;
    float emulSumET_nomu = -9e9;
    float emulMET_nomu_phi = -9e9;

    emulMET_nomu_ex = ((*m_emulmet_nomu)["Final"]->mpx()) * 0.001;
    emulMET_nomu_ey = ((*m_emulmet_nomu)["Final"]->mpy()) * 0.001;
    emulMET_nomu = ((*m_emulmet_nomu)["Final"]->met()) * 0.001;
    emulSumET_nomu = ((*m_emulmet_nomu)["Final"]->sumet()) * 0.001;
    emulMET_nomu_phi = ((*m_emulmet_nomu)["Final"]->phi());

    /*
    h_emulmet_nomu_ex->Fill( emulMET_nomu_ex ); // GeV
    h_emulmet_nomu_ey->Fill( emulMET_nomu_ey ); // GeV
    h_emulmet_nomu->Fill( emulMET_nomu ); // GeV
    h_emulsumet_nomu->Fill( emulSumET_nomu ); // GeV
    h_emulmet_nomu_phi->Fill( emulMET_nomu_phi ); // GeV
    */


    //-----------------------------------------------
    // Define Good Leptons and Calculate Scale Factor
    //-----------------------------------------------

    ///////////////
    // Good Muon //
    ///////////////
    xAOD::MuonContainer* m_goodMuon = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);
    // iterate over our shallow copy
    for (const auto& muon : *muonSC) { // C++11 shortcut
      // Muon Selection for VBF study
      if (dec_baseline(*muon)) {
        if(m_doORtool){
          if (!overlapAcc(*muon)){
            m_goodMuon->push_back( muon );
          }
        }
        else m_goodMuon->push_back( muon );
        //Info("execute()", "  Good muon pt = %.2f GeV", (muon->pt() * 0.001));
      }
    } // end for loop over shallow copied muons

    ///////////////////
    // Good Electron //
    ///////////////////
    xAOD::ElectronContainer* m_goodElectron = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
    // iterate over our shallow copy
    for (const auto& electron : *elecSC) { // C++11 shortcut
      // Electron Selection for VBF study
      if (dec_baseline(*electron)) {
        if(m_doORtool){
          if(!overlapAcc(*electron)){
            m_goodElectron->push_back( electron );
          }
        }
        else m_goodElectron->push_back( electron );
        //Info("execute()", "  Good electron pt = %.2f GeV", (electron->pt() * 0.001));
      }
    } // end for loop over shallow copied electrons

    //////////////
    // Good Tau //
    //////////////
    xAOD::TauJetContainer* m_goodTau = new xAOD::TauJetContainer(SG::VIEW_ELEMENTS);
    // iterate over our shallow copy
    for (const auto& taujet : *tauSC) { // C++11 shortcut
      // Muon Overlap removal
      bool closeToMuonOrEle = false;
      for (const auto& muon : *m_goodMuon) {
       if (DeltaR(taujet->eta(), muon->eta(), taujet->phi(), muon->phi()) < 0.3) closeToMuonOrEle = true;
      }
      for (const auto& electron : *m_goodElectron) {
        if (DeltaR(taujet->eta(), electron->eta(), taujet->phi(), electron->phi()) < 0.3) closeToMuonOrEle = true;
      }
      // Tau Selection for VBF study
      if (dec_baseline(*taujet)) {
        if(m_doORtool){
          if(!overlapAcc(*taujet)){
            m_goodTau->push_back( taujet );
          }
        }
        else {
          if (!closeToMuonOrEle) m_goodTau->push_back( taujet );
        }
        //Info("execute()", "  Good tau pt = %.2f GeV", (taujet->pt() * 0.001));
      }
    } // end for loop over shallow copied electrons


    /////////////////////////
    // Overlap removed Jet //
    /////////////////////////
    xAOD::JetContainer* m_signalJet = new xAOD::JetContainer(SG::VIEW_ELEMENTS);
    // iterate over our shallow copy
    for (const auto& jet : *m_goodJet) { // C++11 shortcut

      if (m_doORmanual) { // If jets are selected jet for overlap removal and if you want overlap removal manually
        bool isORjet = false;

        for (const auto& muon : *m_goodMuon) {
          if (DeltaR(jet->eta(), muon->eta(), jet->phi(), muon->phi()) < m_ORJETdeltaR) {
            int ntrks = 0;
            std::vector<int> ntrks_vec = jet->auxdata<std::vector<int> >("NumTrkPt1000");          
            if (ntrks_vec.size() > 0) ntrks = ntrks_vec[primVertex->index()];
            if (ntrks < 5) isORjet = true;
          }
        }

        for (const auto& electron : *m_goodElectron) {
          if (DeltaR(jet->eta(), electron->eta(), jet->phi(), electron->phi()) < m_ORJETdeltaR) isORjet = true;
        }

        for (const auto& tau : *m_goodTau) {
          if (DeltaR(jet->eta(), tau->eta(), jet->phi(), tau->phi()) < m_ORJETdeltaR) isORjet = true;
        }

        if ( !isORjet ) m_signalJet->push_back( jet );
      }
      else m_signalJet->push_back( jet );

    } // end for loop over shallow copied jets



    //-------------------------------------
    // Define Monojet and DiJet Properties
    //-------------------------------------

    float monojet_pt = 0;
    float monojet_phi = 0;
    float monojet_eta = 0;
    float monojet_rapidity = 0;
    float dPhiMonojetMet = 0;
    float dPhiMonojetMet_nomu = 0;
    float dPhiMonojetMet_noelec = 0;
    TLorentzVector jet1;
    TLorentzVector jet2;
    float jet1_pt = 0;
    float jet2_pt = 0;
    float jet3_pt = 0;
    float jet1_phi = 0;
    float jet2_phi = 0;
    float jet3_phi = 0;
    float jet1_eta = 0;
    float jet2_eta = 0;
    float jet3_eta = 0;
    float jet1_rapidity = 0;
    float jet2_rapidity = 0;
    float jet3_rapidity = 0;

    float signalJet_ht = 0;
    float dPhiJet1Met = 0;
    float dPhiJet2Met = 0;
    float dPhiJet3Met = 0;
    float dPhiJet1Met_nomu = 0;
    float dPhiJet2Met_nomu = 0;
    float dPhiJet3Met_nomu = 0;
    float dPhiJet1Met_noelec = 0;
    float dPhiJet2Met_noelec = 0;
    float dPhiJet3Met_noelec = 0;

    float mjj = 0;
    bool pass_monoJet = false; // Select monoJet
    bool pass_diJet = false; // Select DiJet
    bool pass_CJV = true; // Central Jet Veto (CJV)
    bool pass_dPhijetmet = true; // deltaPhi(Jet_i,MET)
    bool pass_dPhijetmet_nomu = true; // deltaPhi(Jet_i,MET_nomu)
    bool pass_dPhijetmet_noelec = true; // deltaPhi(Jet,MET_noelec)
    float dPhiMinjetmet = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)
    float dPhiMinjetmet_nomu = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)
    float dPhiMinjetmet_noelec = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)


    // Sort signal Jet
    std::sort(m_signalJet->begin(), m_signalJet->end(), DescendingPt());

    // Monojet Selection
    if (m_signalJet->size() > 0) {

      monojet_pt = m_signalJet->at(0)->pt() * 0.001;
      monojet_phi = m_signalJet->at(0)->phi();
      monojet_eta = m_signalJet->at(0)->eta();
      monojet_rapidity = m_signalJet->at(0)->rapidity();


      // Define Monojet
      if ( monojet_pt >  m_monoJetPtCut ){
        if ( fabs(monojet_eta) < m_monoJetEtaCut){
          if ( m_jetCleaningTight->accept( *m_signalJet->at(0) ) ){ //Tight Leading Jet 
            pass_monoJet = true;
            //Info("execute()", "  Leading jet pt = %.2f GeV", monojet_pt);
          }
        }
      }

      // deltaPhi(monojet,MET) decision
      // For Znunu
      if (m_isZnunu){
        dPhiMonojetMet = DeltaPhi(monojet_phi, MET_phi);
      }
      // For Muon channel
      if (m_isZmumu || m_isWmunu){
        dPhiMonojetMet_nomu = DeltaPhi(monojet_phi, emulMET_nomu_phi);
      }
      // For Electron channel
      if (m_isZee || m_isWenu){
        dPhiMonojetMet_noelec = DeltaPhi(monojet_phi, emulMET_noelec_phi);
      }

    } // MonoJet selection 



    // DiJet Selection
    if (m_signalJet->size() > 1) {

      jet1 = m_signalJet->at(0)->p4();
      jet2 = m_signalJet->at(1)->p4();
      jet1_pt = m_signalJet->at(0)->pt() * 0.001;
      jet2_pt = m_signalJet->at(1)->pt() * 0.001;
      jet1_phi = m_signalJet->at(0)->phi();
      jet2_phi = m_signalJet->at(1)->phi();
      jet1_eta = m_signalJet->at(0)->eta();
      jet2_eta = m_signalJet->at(1)->eta();
      jet1_rapidity = m_signalJet->at(0)->rapidity();
      jet2_rapidity = m_signalJet->at(1)->rapidity();
      auto dijet = jet1 + jet2;
      mjj = dijet.M() * 0.001;

      //Info("execute()", "  jet1 = %.2f GeV, jet2 = %.2f GeV", jet1_pt, jet2_pt);
      //Info("execute()", "  mjj = %.2f GeV", mjj);

      // Define Dijet
      if ( jet1_pt > m_diJet1PtCut && jet2_pt > m_diJet2PtCut ){
        if ( fabs(jet1_rapidity) < m_diJetRapCut && fabs(jet2_rapidity) < m_diJetRapCut ){
          if ( m_jetCleaningTight->accept( *m_signalJet->at(0) ) ){ //Tight Leading Jet 
            pass_diJet = true;
          }
        }
      }

      // deltaPhi(Jet1,MET) or deltaPhi(Jet2,MET) decision
      // For Znunu
      if (m_isZnunu){
        dPhiJet1Met = DeltaPhi(jet1_phi, MET_phi);
        dPhiJet2Met = DeltaPhi(jet2_phi, MET_phi);
      }
      // For Muon channel
      if (m_isZmumu || m_isWmunu){
        dPhiJet1Met_nomu = DeltaPhi(jet1_phi, emulMET_nomu_phi);
        dPhiJet2Met_nomu = DeltaPhi(jet2_phi, emulMET_nomu_phi);
      }
      // For Electron channel
      if (m_isZee || m_isWenu){
        dPhiJet1Met_noelec = DeltaPhi(jet1_phi, emulMET_noelec_phi);
        dPhiJet2Met_noelec = DeltaPhi(jet2_phi, emulMET_noelec_phi);
      }

    } // DiJet selection 




    // For jet3
    if (m_signalJet->size() > 2) {
      jet3_pt = m_signalJet->at(2)->pt() * 0.001;
      jet3_phi = m_signalJet->at(2)->phi();
      jet3_eta = m_signalJet->at(2)->eta();
      jet3_rapidity = m_signalJet->at(2)->rapidity();
      // deltaPhi(Jet3,MET)
      dPhiJet3Met = DeltaPhi(jet3_phi, MET_phi);
      dPhiJet3Met_nomu = DeltaPhi(jet3_phi, emulMET_nomu_phi);
      dPhiJet3Met_noelec = DeltaPhi(jet3_phi, emulMET_noelec_phi);
    }


    // Define DeltaPhi(Jet_i,MET) cut and Central Jet Veto (CJV)
    if (m_signalJet->size() > 0) {

      // loop over the jets in the Good Jets Container
      for (const auto& jet : *m_signalJet) {
        float signal_jet_pt = jet->pt() * 0.001;
        float signal_jet_rapidity = jet->rapidity();
        float signal_jet_phi = jet->phi();

        // Calculate dPhi(Jet_i,MET) and dPhi_min(Jet_i,MET)
        // For Znunu
        if (m_isZnunu){
          float dPhijetmet = DeltaPhi(signal_jet_phi,MET_phi);
          //Info("execute()", " [Znunu] Event # = %llu", eventInfo->eventNumber());
          //Info("execute()", " [Znunu] dPhi = %.2f", dPhijetmet);
          if ( signal_jet_pt > 30. && fabs(signal_jet_rapidity) < 4.4 && dPhijetmet < 0.4 ) pass_dPhijetmet = false;
          dPhiMinjetmet = std::min(dPhiMinjetmet, dPhijetmet);
          //Info("execute()", " [Znunu] dPhi_min = %.2f", dPhiMinjetmet);
        }
        // For muon channel
        if (m_isZmumu || m_isWmunu){
          float dPhijetmet_nomu = DeltaPhi(signal_jet_phi,emulMET_nomu_phi);
          //Info("execute()", " [Zmumu] Event # = %llu", eventInfo->eventNumber());
          //Info("execute()", " [Zmumu] dPhi = %.2f", dPhijetmet_nomu);
          if ( signal_jet_pt > 30. && fabs(signal_jet_rapidity) < 4.4 && dPhijetmet_nomu < 0.4 ) pass_dPhijetmet_nomu = false;
          dPhiMinjetmet_nomu = std::min(dPhiMinjetmet_nomu, dPhijetmet_nomu);
          //Info("execute()", " [Zmumu] dPhi_min = %.2f", dPhiMinjetmet_nomu);
        }
        // For electron channel
        if (m_isZee || m_isWenu){
          float dPhijetmet_noelec = DeltaPhi(signal_jet_phi,emulMET_noelec_phi);
          //Info("execute()", " [Zee] Event # = %llu", eventInfo->eventNumber());
          //Info("execute()", " [Zee] dPhi = %.2f", dPhijetmet_noelec);
          if ( signal_jet_pt > 30. && fabs(signal_jet_rapidity) < 4.4 && dPhijetmet_noelec < 0.4 ) pass_dPhijetmet_noelec = false;
          dPhiMinjetmet_noelec = std::min(dPhiMinjetmet_noelec, dPhijetmet_noelec);
          //Info("execute()", " [Zee] dPhi_min = %.2f", dPhiMinjetmet_noelec);
        }

        // Central Jet Veto (CJV)
        if ( m_signalJet->size() > 2 && pass_diJet ){
          if (m_signalJet->at(0) != jet && m_signalJet->at(1) != jet){
            //cout << "m_signalJet->at(0) = " << m_signalJet->at(0) << " jet = " << jet << endl;
            if (signal_jet_pt > m_CJVptCut && fabs(signal_jet_rapidity) < m_diJetRapCut) {
              if ( (jet1_rapidity > jet2_rapidity) && (signal_jet_rapidity < jet1_rapidity && signal_jet_rapidity > jet2_rapidity)){
                pass_CJV = false;
              }
              if ( (jet1_rapidity < jet2_rapidity) && (signal_jet_rapidity > jet1_rapidity && signal_jet_rapidity < jet2_rapidity)){
                pass_CJV = false;
              }
            }
          }
        }

        //Info("execute()", "  Znunu Signal Jet pt = %.2f GeV, eta = %.2f", signal_pt_jet, signal_eta_jet);
        signalJet_ht += signal_jet_pt;
      } // Jet loop

    } // End DeltaPhi(Jet_i,MET) cut and Central Jet Veto (CJV)




    //----------------------------------
    // Define Zmumu and Wmunu Selection
    //----------------------------------

    // Wmunu
    bool pass_Wmunu = false; // Select Wmunu channel
    // Zmumu
    float mll_muon = 0.; // For Zmumu channel
    float mT_muon = 0.;// For Wmunu channel
    float muon1_pt = 0.;
    float muon2_pt = 0.;
    float muon1_charge = 0.;
    float muon2_charge = 0.;
    bool pass_dimuonPtCut = false; // dimuon pT cut
    bool pass_OSmuon = false; // Opposite sign change muon
    bool pass_SSmuon = false; // Same sign change muon

    if (m_isZmumu || m_isWmunu){

      // Zmumu Selection
      if (m_goodMuon->size() > 1) {
        std::partial_sort(m_goodMuon->begin(), m_goodMuon->begin()+2, m_goodMuon->end(), DescendingPt());

        TLorentzVector muon1 = m_goodMuon->at(0)->p4();
        TLorentzVector muon2 = m_goodMuon->at(1)->p4();
        muon1_pt = m_goodMuon->at(0)->pt() * 0.001;
        muon2_pt = m_goodMuon->at(1)->pt() * 0.001;
        muon1_charge = m_goodMuon->at(0)->charge();
        muon2_charge = m_goodMuon->at(1)->charge();
        auto Zmass_muon = muon1 + muon2;
        mll_muon = Zmass_muon.M() * 0.001;

        //Info("execute()", "  muon1 = %.2f GeV, muon2 = %.2f GeV", muon1_pt, muon2_pt);
        //Info("execute()", "  mll (Zmumu) = %.2f GeV", mll_muon);

        if ( muon1_pt >  m_LeadLepPtCut && muon2_pt > m_SubLeadLepPtCut ) pass_dimuonPtCut = true;
        if ( muon1_charge * muon2_charge < 0 ) pass_OSmuon = true;
        if ( muon1_charge * muon2_charge > 0 ) pass_SSmuon = true;
        //Info("execute()", "  muon1 charge = %f, muon2 charge = %f, pass_OSmuon = %d, pass_SSmuon = %d", muon1_charge, muon2_charge, pass_OSmuon, pass_SSmuon);

      } // Zmumu selection loop

      // Wmunu Selection
      if (m_goodMuon->size() == 1) {
        float muon_pt = m_goodMuon->at(0)->pt() * 0.001;
        float muon_phi = m_goodMuon->at(0)->phi();
        mT_muon = TMath::Sqrt( 2. * muon_pt * MET * ( 1. - TMath::Cos(muon_phi - MET_phi) ) );

        if ( muon_pt > 25. ){
          pass_Wmunu = true;
        }
      } //Wmunu Selection loop

    }


    //-------------------------------
    // Define Zee and Wenu Selection
    // ------------------------------

    // Wenu
    bool pass_Wenu = false; // Select Wenu channel
    // Zee
    float mll_electron = 0.; // Select Zee channel
    float mT_electron = 0.;// For Wenu channel
    float electron1_pt = 0.;
    float electron2_pt = 0.;
    float electron1_charge = 0.;
    float electron2_charge = 0.;
    bool pass_dielectronPtCut = false; // dielectron pT cut
    bool pass_OSelectron = false; // Opposite sign change electron
    bool pass_SSelectron = false; // Same sign change electron

    if (m_isZee || m_isWenu){

      // Zee Selection
      if (m_goodElectron->size() > 1) {
        std::partial_sort(m_goodElectron->begin(), m_goodElectron->begin()+2, m_goodElectron->end(), DescendingPt());

        TLorentzVector electron1 = m_goodElectron->at(0)->p4();
        TLorentzVector electron2 = m_goodElectron->at(1)->p4();
        electron1_pt = m_goodElectron->at(0)->pt() * 0.001;
        electron2_pt = m_goodElectron->at(1)->pt() * 0.001;
        electron1_charge = m_goodElectron->at(0)->charge();
        electron2_charge = m_goodElectron->at(1)->charge();
        auto Zmass_electron = electron1 + electron2;
        mll_electron = Zmass_electron.M() * 0.001;

        //Info("execute()", "  electron1 = %.2f GeV, electron2 = %.2f GeV", electron1_pt, electron2_pt);
        //Info("execute()", "  mll (Zee) = %.2f GeV", mll_electron);


        if ( electron1_pt >  m_LeadLepPtCut && electron2_pt > m_SubLeadLepPtCut ) pass_dielectronPtCut = true;
        if ( electron1_charge * electron2_charge < 0 ) pass_OSelectron = true;
        if ( electron1_charge * electron2_charge > 0 ) pass_SSelectron = true;

      } // Zee selection loop

      // Wenu Selection
      if (m_goodElectron->size() == 1) {
        float electron_pt = m_goodElectron->at(0)->pt() * 0.001;
        float electron_phi = m_goodElectron->at(0)->phi();
        mT_electron = TMath::Sqrt( 2. * electron_pt * MET * ( 1. - TMath::Cos(electron_phi - MET_phi) ) );

        if ( electron_pt > 25. ){
          pass_Wenu = true;
        }
      } //Wenu Selection loop

    }




    // ------------------
    // Get isolated track
    // ------------------
    /*
    // Retrieve main TrackParticle collection
    const xAOD::TrackParticleContainer* inTracks(0);
    if ( !m_event->retrieve( inTracks, "InDetTrackParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve TrackParticle container. Exiting." );
    return EL::StatusCode::FAILURE;
    }

    int Nisotrk = NumIsoTracks(inTracks, primVertex, 3., 10.) - NumMuonIsoTrack(muonSC, inTracks, primVertex, 3., 10.) - NumElecIsoTrack(elecSC, inTracks, primVertex, 3., 10.);

    bool passIsoTrk = true;
    if (Nisotrk > 0) {
    passIsoTrk = false;
    //Info("execute()", "  The number of Isolated track counted = %i (N_SignalMuon = %lu, N_SignalElec = %lu)", Nisotrk, m_signalMuon->size(), m_signalElectron->size() );
    }
    */



    //-----------
    // VBF study 
    //-----------

    //-------------------------------
    // Z -> vv + JET EVENT SELECTION
    //-------------------------------

    if (m_isZnunu){
      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]MET Trigger");
        m_eventCutflow[5]+=1;
        if ( MET > m_metCut ) {
          if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]MET cut");
          m_eventCutflow[6]+=1;
          /*
             Info("execute()", "=====================================");
             Info("execute()", " Event # = %llu", eventInfo->eventNumber());
             Info("execute()", " Good Event number = %i", m_eventCutflow[6]);
             Info("execute()", " MET = %.3f GeV", MET);
             Info("execute()", " RefElectron = %.3f GeV", ((*m_met)["RefElectron"]->met()) * 0.001);
             Info("execute()", " RefPhoton = %.3f GeV", ((*m_met)["RefPhoton"]->met()) * 0.001);
             Info("execute()", " RefTau = %.3f GeV", ((*m_met)["RefTau"]->met()) * 0.001);
             Info("execute()", " RefMuon = %.3f GeV", ((*m_met)["RefMuon"]->met()) * 0.001);
             Info("execute()", " RefJet = %.3f GeV", ((*m_met)["RefJet"]->met()) * 0.001);
             Info("execute()", " SoftClus = %.3f GeV", ((*m_met)["SoftClus"]->met()) * 0.001);
             Info("execute()", " PVSoftTrk = %.3f GeV", ((*m_met)["PVSoftTrk"]->met()) * 0.001);
             Info("execute()", " # of good jets = %lu", m_signalJet->size());
             if (m_signalJet->size() > 0){
             int jetCount = 0;
             for (const auto& jet : *m_signalJet) {
             jetCount++;
             Info("execute()", " jet # : %i", jetCount);
             Info("execute()", " jet pt = %.3f GeV", jet->pt() * 0.001);
             Info("execute()", " jet eta = %.3f GeV", jet->eta());
             Info("execute()", " jet phi = %.3f GeV", jet->phi());
             }
             }
             if (m_goodElectron->size() > 0){
             int eleCount = 0;
             for (const auto& electron : *m_goodElectron) {
             eleCount++;
             Info("execute()", " electron # : %i", eleCount);
             Info("execute()", " electron pt = %.3f GeV", electron->pt() * 0.001);
             Info("execute()", " electron eta = %.3f GeV", electron->eta());
             Info("execute()", " electron phi = %.3f GeV", electron->phi());
             }
             }
             if (m_goodMuon->size() > 0){
             int muCount = 0;
             for (const auto& muon : *m_goodMuon) {
             muCount++;
             Info("execute()", " muon # : %i", muCount);
             Info("execute()", " muon pt = %.3f GeV", muon->pt() * 0.001);
             Info("execute()", " muon eta = %.3f GeV", muon->eta());
             Info("execute()", " muon phi = %.3f GeV", muon->phi());
             }
             }
             if (m_goodTau->size() > 0){
             int tauCount = 0;
             for (const auto& tau : *m_goodTau) {
             tauCount++;
             Info("execute()", " tau # : %i", tauCount);
             Info("execute()", " tau pt = %.3f GeV", tau->pt() * 0.001);
             Info("execute()", " tau eta = %.3f GeV", tau->eta());
             Info("execute()", " tau phi = %.3f GeV", tau->phi());
             }
             }
             */
          if (m_goodElectron->size() == 0) {
            if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]Electron Veto");
            m_eventCutflow[7]+=1;
            if ( m_goodMuon->size() == 0) {
              if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]Muon Veto");
              m_eventCutflow[8]+=1;
              if (m_goodTau->size() == 0) {
                if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]Tau Veto");
                m_eventCutflow[9]+=1;
                if ( m_signalJet->size() > 0 ) {
                  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]At least One Jets");
                  m_eventCutflow[10]+=1;
                  if ( pass_diJet ) {
                    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]DiJet");
                    m_eventCutflow[11]+=1;
                    if ( mjj > m_mjjCut ) {
                      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]mjj cut");
                      m_eventCutflow[12]+=1;
                      if ( pass_CJV ) {
                        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]CJV");
                        m_eventCutflow[13]+=1;
                        if ( pass_dPhijetmet ) {
                          if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]dPhi(jet_i,MET) cut");
                          m_eventCutflow[14]+=1;
                          /*
                          // Fill histogram
                          // MET
                          h_znunu_met->Fill(MET);
                          // Jets
                          h_znunu_njet->Fill(m_signalJet->size());
                          h_znunu_jet1_pt->Fill(jet1_pt);
                          h_znunu_jet2_pt->Fill(jet2_pt);
                          h_znunu_jet1_phi->Fill(jet1_phi);
                          h_znunu_jet2_phi->Fill(jet2_phi);
                          h_znunu_jet1_eta->Fill(jet1_eta);
                          h_znunu_jet2_eta->Fill(jet2_eta);
                          h_znunu_jet1_rap->Fill(jet1_rapidity);
                          h_znunu_jet2_rap->Fill(jet2_rapidity);
                          h_znunu_mjj->Fill(mjj);
                          h_znunu_dPhijj->Fill(DeltaPhi(jet1_phi, jet2_phi));
                          h_znunu_dRjj->Fill(DeltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi));
                          h_znunu_dPhimetj1->Fill(dPhiJet1Met);
                          h_znunu_dPhimetj2->Fill(dPhiJet2Met);
                          h_znunu_dPhiMinmetjet->Fill(dPhiMinjetmet);
                          // For jet3
                          if (m_signalJet->size() > 2){
                            h_znunu_jet3_pt->Fill(jet3_pt);
                            h_znunu_jet3_phi->Fill(jet3_phi);
                            h_znunu_jet3_eta->Fill(jet3_eta);
                            h_znunu_jet3_rap->Fill(jet3_rapidity);
                            h_znunu_dPhimetj3->Fill(dPhiJet3Met);
                          }
                          */


                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    //---------------------------------
    // Z -> mumu + JET EVENT SELECTION
    //---------------------------------

    if (m_isZmumu){
      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]MET Trigger");
        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[16]+=1;
        if ( emulMET_nomu > m_metCut ) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]MET cut");
          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[17]+=1;
          if (m_goodElectron->size() == 0) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]Electron Veto");
            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[18]+=1;
            if ( m_goodMuon->size() > 1) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]At least Two Electrons");
              if (sysName == "" && m_useArrayCutflow) m_eventCutflow[19]+=1;
              if (m_goodTau->size() == 0) {
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]Tau Veto");
                if (sysName == "" && m_useArrayCutflow) m_eventCutflow[20]+=1;
                if ( pass_dimuonPtCut && pass_OSmuon && m_goodMuon->size() == 2 && mll_muon > 66. && mll_muon < 116. ){
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]mll cut");
                  if (sysName == "" && m_useArrayCutflow) m_eventCutflow[21]+=1;
                  if ( m_signalJet->size() > 0 ) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]At least One Jets");
                    if (sysName == "" && m_useArrayCutflow) m_eventCutflow[22]+=1;

                    // Calculate muon SF
                    if (!m_isData) {
                      double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoSF, m_ttvaSF);
                      //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                      mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                    }

                    ////////////////////////
                    // MonoJet phasespace //
                    ////////////////////////
                    if ( pass_monoJet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, monojet]MonoJet");
                      if ( pass_dPhijetmet_nomu ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, monojet]dPhi(jet_i,MET) cut");

                        // Fill histogram
                        h_channel = "h_zmumu_";
                        // Publication plot
                        hMap1D["ZmumuMET_mono"+sysName]->Fill(emulMET_nomu, mcEventWeight_Zmumu);
                        // MET
                        hMap1D[h_channel+"monojet_met"+sysName]->Fill(MET, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_met_emulmet"+sysName]->Fill(emulMET_nomu, mcEventWeight_Zmumu);
                        // Jets
                        hMap1D[h_channel+"monojet_njet"+sysName]->Fill(m_signalJet->size(), mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_jet_pt"+sysName]->Fill(monojet_pt, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_jet_phi"+sysName]->Fill(monojet_phi, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_jet_eta"+sysName]->Fill(monojet_eta, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_jet_rap"+sysName]->Fill(monojet_rapidity, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_dPhimetjet"+sysName]->Fill(dPhiMonojetMet_nomu, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_nomu, mcEventWeight_Zmumu);
                        // Leptons
                        hMap1D[h_channel+"monojet_lepton1_pt"+sysName]->Fill(m_goodMuon->at(0)->pt() * 0.001, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_lepton2_pt"+sysName]->Fill(m_goodMuon->at(1)->pt() * 0.001, mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_lepton1_phi"+sysName]->Fill(m_goodMuon->at(0)->phi(), mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_lepton2_phi"+sysName]->Fill(m_goodMuon->at(1)->phi(), mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_lepton1_eta"+sysName]->Fill(m_goodMuon->at(0)->eta(), mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_lepton2_eta"+sysName]->Fill(m_goodMuon->at(1)->eta(), mcEventWeight_Zmumu);
                        hMap1D[h_channel+"monojet_mll"+sysName]->Fill(mll_muon, mcEventWeight_Zmumu);

                      } // pass dPhijetmet_noelec
                    } // pass monojet

                    ////////////////////
                    // VBF phasespace //
                    ////////////////////
                    if ( pass_diJet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, VBF]DiJet");
                      if (sysName == "" && m_useArrayCutflow) m_eventCutflow[23]+=1;
                      if ( mjj > m_mjjCut ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, VBF]mjj cut");
                        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[24]+=1;
                        if ( pass_CJV ) {
                          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, VBF]CJV cut");
                          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[25]+=1;
                          if ( pass_dPhijetmet_nomu ) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, VBF]dPhi(jet_i,MET) cut");
                            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[26]+=1;
                            // Fill histogram
                            h_channel = "h_zmumu_";
                            // Publication plot
                            hMap1D["ZmumuMET_search"+sysName]->Fill(emulMET_nomu, mcEventWeight_Zmumu);
                            hMap1D["ZmumuMjj_search"+sysName]->Fill(mjj, mcEventWeight_Zmumu);
                            hMap1D["ZmumuDeltaPhiAll"+sysName]->Fill(DeltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                            // MET
                            hMap1D[h_channel+"vbf_met"+sysName]->Fill(MET, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_met_emulmet"+sysName]->Fill(emulMET_nomu, mcEventWeight_Zmumu);
                            // Jets
                            hMap1D[h_channel+"vbf_mjj"+sysName]->Fill(mjj, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_dPhijj"+sysName]->Fill(DeltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_njet"+sysName]->Fill(m_signalJet->size(), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet1_pt"+sysName]->Fill(jet1_pt, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet2_pt"+sysName]->Fill(jet2_pt, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet1_phi"+sysName]->Fill(jet1_phi, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet2_phi"+sysName]->Fill(jet2_phi, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet1_eta"+sysName]->Fill(jet1_eta, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet2_eta"+sysName]->Fill(jet2_eta, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet1_rap"+sysName]->Fill(jet1_rapidity, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_jet2_rap"+sysName]->Fill(jet2_rapidity, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_dRjj"+sysName]->Fill(DeltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_dPhimetj1"+sysName]->Fill(dPhiJet1Met_nomu, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_dPhimetj2"+sysName]->Fill(dPhiJet2Met_nomu, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_nomu, mcEventWeight_Zmumu);
                            // For jet3
                            if (m_signalJet->size() > 2){
                              hMap1D[h_channel+"vbf_jet3_pt"+sysName]->Fill(jet3_pt, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet3_phi"+sysName]->Fill(jet3_phi, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet3_eta"+sysName]->Fill(jet3_eta, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet3_rap"+sysName]->Fill(jet3_rapidity, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_dPhimetj3"+sysName]->Fill(dPhiJet3Met_nomu, mcEventWeight_Zmumu);
                            }
                            // Leptons
                            hMap1D[h_channel+"vbf_lepton1_pt"+sysName]->Fill(m_goodMuon->at(0)->pt() * 0.001, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_lepton2_pt"+sysName]->Fill(m_goodMuon->at(1)->pt() * 0.001, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_lepton1_phi"+sysName]->Fill(m_goodMuon->at(0)->phi(), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_lepton2_phi"+sysName]->Fill(m_goodMuon->at(1)->phi(), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_lepton1_eta"+sysName]->Fill(m_goodMuon->at(0)->eta(), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_lepton2_eta"+sysName]->Fill(m_goodMuon->at(1)->eta(), mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_mll"+sysName]->Fill(mll_muon, mcEventWeight_Zmumu);


                          } // pass diPhijetMET
                        } // pass CJV
                      } // mjj Cut
                    } // pass diJet

                  } // at least 1 jet
                } // pass Zmumu 
              } // Tau veto
            } // at least 1 muon
          } // Electron veto
        } // MET cut
      } // HLT_xe70
    } // m_isZmumu


    //---------------------------------
    // W -> munu + JET EVENT SELECTION
    //---------------------------------

    if (m_isWmunu){
      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
        if ( emulMET_nomu > m_metCut ) {
          if (m_goodElectron->size() == 0) {
            if ( m_goodMuon->size() > 0 ) {
              if (m_goodTau->size() == 0) {
                if ( pass_Wmunu && m_goodMuon->size() == 1 && mT_muon > 30. && mT_muon < 100. ){
                  if ( m_signalJet->size() > 1 ) {
                    if ( pass_diJet ) {
                      if ( mjj > m_mjjCut ) {
                        if ( pass_CJV ) {
                          // Calculate muon SF
                          if (!m_isData) {
                            double totalMuonSF_Wmunu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoSF, m_ttvaSF);
                            //Info("execute()", " Wmunu Total Muon SF = %.3f ", totalMuonSF_Wmunu);
                            mcEventWeight_Wmunu = mcEventWeight * totalMuonSF_Wmunu;
                          }
                          /*
                          // Fill histogram
                          // MET
                          h_wmunu_met->Fill(MET, mcEventWeight_Wmunu);
                          h_wmunu_emulmet_nomu->Fill(emulMET_nomu, mcEventWeight_Wmunu);
                          // Jets
                          h_wmunu_njet->Fill(m_signalJet->size(), mcEventWeight_Wmunu);
                          h_wmunu_jet1_pt->Fill(jet1_pt, mcEventWeight_Wmunu);
                          h_wmunu_jet2_pt->Fill(jet2_pt, mcEventWeight_Wmunu);
                          h_wmunu_jet1_phi->Fill(jet1_phi, mcEventWeight_Wmunu);
                          h_wmunu_jet2_phi->Fill(jet2_phi, mcEventWeight_Wmunu);
                          h_wmunu_jet1_eta->Fill(jet1_eta, mcEventWeight_Wmunu);
                          h_wmunu_jet2_eta->Fill(jet2_eta, mcEventWeight_Wmunu);
                          h_wmunu_jet1_rap->Fill(jet1_rapidity, mcEventWeight_Wmunu);
                          h_wmunu_jet2_rap->Fill(jet2_rapidity, mcEventWeight_Wmunu);
                          h_wmunu_mjj->Fill(mjj, mcEventWeight_Wmunu);
                          h_wmunu_dPhijj->Fill(DeltaPhi(jet1_phi, jet2_phi), mcEventWeight_Wmunu);
                          h_wmunu_dRjj->Fill(DeltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Wmunu);
                          h_wmunu_dPhimetj1->Fill(dPhiJet1Met_nomu, mcEventWeight_Wmunu);
                          h_wmunu_dPhimetj2->Fill(dPhiJet2Met_nomu, mcEventWeight_Wmunu);
                          h_wmunu_dPhiMinmetjet->Fill(dPhiMinjetmet_nomu, mcEventWeight_Wmunu);
                          // For jet3
                          if (m_signalJet->size() > 2){
                          h_wmunu_jet3_pt->Fill(jet3_pt, mcEventWeight_Wmunu);
                          h_wmunu_jet3_phi->Fill(jet3_phi, mcEventWeight_Wmunu);
                          h_wmunu_jet3_eta->Fill(jet3_eta, mcEventWeight_Wmunu);
                          h_wmunu_jet3_rap->Fill(jet3_rapidity, mcEventWeight_Wmunu);
                          h_wmunu_dPhimetj3->Fill(dPhiJet3Met_nomu, mcEventWeight_Wmunu);
                          }
                          // Leptons
                          h_wmunu_muon_pt->Fill(m_goodMuon->at(0)->pt() * 0.001, mcEventWeight_Wmunu);
                          h_wmunu_muon_phi->Fill(m_goodMuon->at(0)->phi(), mcEventWeight_Wmunu);
                          h_wmunu_muon_eta->Fill(m_goodMuon->at(0)->eta(), mcEventWeight_Wmunu);
                          h_wmunu_mT->Fill(mT_muon, mcEventWeight_Wmunu);
                          */

                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }




    //-------------------------------
    // Z -> ee + JET EVENT SELECTION
    //-------------------------------

    if (m_isZee){
      if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]Electron Trigger");
        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[28]+=1;
        if ( emulMET_noelec > m_metCut ) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]MET cut");
          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[29]+=1;
          if (m_goodElectron->size() > 1) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]At least Two Electron");
            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[30]+=1;
            if ( m_goodMuon->size() == 0) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]Muon Veto");
              if (sysName == "" && m_useArrayCutflow) m_eventCutflow[31]+=1;
              if (m_goodTau->size() == 0) {
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]Tau Veto");
                if (sysName == "" && m_useArrayCutflow) m_eventCutflow[32]+=1;
                if ( pass_dielectronPtCut && pass_OSelectron && m_goodElectron->size() == 2 && mll_electron > 66. && mll_electron < 116. ) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]mll cut");
                  if (sysName == "" && m_useArrayCutflow) m_eventCutflow[33]+=1;
                  if ( m_signalJet->size() > 0 ) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]At least One Jets");
                    if (sysName == "" && m_useArrayCutflow) m_eventCutflow[34]+=1;

                    // Calculate electron SF
                    if (!m_isData) {
                      float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoSF, m_trigSF);
                      //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                      mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                    }

                    ////////////////////////
                    // MonoJet phasespace //
                    ////////////////////////
                    if ( pass_monoJet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, monojet]MonoJet");
                      if ( pass_dPhijetmet_noelec ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, monojet]dPhi(jet_i,MET) cut");

                        // Fill histogram
                        h_channel = "h_zee_";
                        // Publication plot
                        hMap1D["ZeeMET_mono"+sysName]->Fill(emulMET_noelec, mcEventWeight_Zee);
                        // MET
                        hMap1D[h_channel+"monojet_met"+sysName]->Fill(MET, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_met_emulmet"+sysName]->Fill(emulMET_noelec, mcEventWeight_Zee);
                        // Jets
                        hMap1D[h_channel+"monojet_njet"+sysName]->Fill(m_signalJet->size(), mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_jet_pt"+sysName]->Fill(monojet_pt, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_jet_phi"+sysName]->Fill(monojet_phi, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_jet_eta"+sysName]->Fill(monojet_eta, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_jet_rap"+sysName]->Fill(monojet_rapidity, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_dPhimetjet"+sysName]->Fill(dPhiMonojetMet_noelec, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_noelec, mcEventWeight_Zee);
                        // Leptons
                        hMap1D[h_channel+"monojet_lepton1_pt"+sysName]->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_lepton2_pt"+sysName]->Fill(m_goodElectron->at(1)->pt() * 0.001, mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_lepton1_phi"+sysName]->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_lepton2_phi"+sysName]->Fill(m_goodElectron->at(1)->phi(), mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_lepton1_eta"+sysName]->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_lepton2_eta"+sysName]->Fill(m_goodElectron->at(1)->eta(), mcEventWeight_Zee);
                        hMap1D[h_channel+"monojet_mll"+sysName]->Fill(mll_electron, mcEventWeight_Zee);

                      } // pass dPhijetmet_noelec
                    } // pass monojet

                    ////////////////////
                    // VBF phasespace //
                    ////////////////////
                    if ( pass_diJet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, VBF]DiJet");
                      if (sysName == "" && m_useArrayCutflow) m_eventCutflow[35]+=1;
                      if ( mjj > m_mjjCut ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, VBF]mjj cut");
                        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[36]+=1;
                        if ( pass_CJV ) {
                          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, VBF]CJV cut");
                          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[37]+=1;
                          if ( pass_dPhijetmet_noelec ) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, VBF]dPhi(jet_i,MET) cut");
                            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[38]+=1;

                            // Fill histogram
                            h_channel = "h_zee_";
                            // Publication plot
                            hMap1D["ZeeMET_search"+sysName]->Fill(emulMET_noelec, mcEventWeight_Zee);
                            hMap1D["ZeeMjj_search"+sysName]->Fill(mjj, mcEventWeight_Zee);
                            hMap1D["ZeeDeltaPhiAll"+sysName]->Fill(DeltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                            // MET
                            hMap1D[h_channel+"vbf_met"+sysName]->Fill(MET, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_met_emulmet"+sysName]->Fill(emulMET_noelec, mcEventWeight_Zee);
                            // Jets
                            hMap1D[h_channel+"vbf_mjj"+sysName]->Fill(mjj, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_dPhijj"+sysName]->Fill(DeltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_njet"+sysName]->Fill(m_signalJet->size(), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet1_pt"+sysName]->Fill(jet1_pt, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet2_pt"+sysName]->Fill(jet2_pt, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet1_phi"+sysName]->Fill(jet1_phi, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet2_phi"+sysName]->Fill(jet2_phi, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet1_eta"+sysName]->Fill(jet1_eta, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet2_eta"+sysName]->Fill(jet2_eta, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet1_rap"+sysName]->Fill(jet1_rapidity, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_jet2_rap"+sysName]->Fill(jet2_rapidity, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_dRjj"+sysName]->Fill(DeltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_dPhimetj1"+sysName]->Fill(dPhiJet1Met_noelec, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_dPhimetj2"+sysName]->Fill(dPhiJet2Met_noelec, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_noelec, mcEventWeight_Zee);
                            // For jet3
                            if (m_signalJet->size() > 2){
                              hMap1D[h_channel+"vbf_jet3_pt"+sysName]->Fill(jet3_pt, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet3_phi"+sysName]->Fill(jet3_phi, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet3_eta"+sysName]->Fill(jet3_eta, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet3_rap"+sysName]->Fill(jet3_rapidity, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_dPhimetj3"+sysName]->Fill(dPhiJet3Met_noelec, mcEventWeight_Zee);
                            }
                            // Leptons
                            hMap1D[h_channel+"vbf_lepton1_pt"+sysName]->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_lepton2_pt"+sysName]->Fill(m_goodElectron->at(1)->pt() * 0.001, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_lepton1_phi"+sysName]->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_lepton2_phi"+sysName]->Fill(m_goodElectron->at(1)->phi(), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_lepton1_eta"+sysName]->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_lepton2_eta"+sysName]->Fill(m_goodElectron->at(1)->eta(), mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_mll"+sysName]->Fill(mll_electron, mcEventWeight_Zee);

                          } // pass dPhijetmet_noelec
                        } // pass CJV
                      } // mjj cut
                    } // pass diJet

                  } // at least 1 jet
                } // pass Zee
              } // Tau Veto
            } // Muon Veto
          } // at least 1 electron
        } // MET cut
      } // sigle electron trigger
    } // m_isZee



    //---------------------------------
    // W -> enu + JET EVENT SELECTION
    //---------------------------------

    if (m_isWenu){
      if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){
        if ( emulMET_noelec > m_metCut ) {
          if (m_goodElectron->size() > 0) {
            if ( m_goodMuon->size() == 0 ) {
              if (m_goodTau->size() == 0) {
                if ( pass_Wenu && m_goodElectron->size() == 1 && mT_electron > 30. && mT_electron < 100. ){
                  if ( m_signalJet->size() > 1 ) {
                    if ( pass_diJet ) {
                      if ( mjj > m_mjjCut ) {
                        if ( pass_CJV ) {
                          // Calculate electron SF
                          if (!m_isData) {
                            double totalElectronSF_Wenu = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoSF, m_trigSF);
                            //Info("execute()", " Wenu Total Electron SF = %.3f ", totalElectronSF_Wenu);
                            mcEventWeight_Wenu = mcEventWeight * totalElectronSF_Wenu;
                          }
                          /*
                          // Fill histogram
                          // MET
                          h_wenu_met->Fill(MET, mcEventWeight_Wenu);
                          h_wenu_emulmet_noelec->Fill(emulMET_noelec, mcEventWeight_Wenu);
                          // Jets
                          h_wenu_njet->Fill(m_signalJet->size(), mcEventWeight_Wenu);
                          h_wenu_jet1_pt->Fill(jet1_pt, mcEventWeight_Wenu);
                          h_wenu_jet2_pt->Fill(jet2_pt, mcEventWeight_Wenu);
                          h_wenu_jet1_phi->Fill(jet1_phi, mcEventWeight_Wenu);
                          h_wenu_jet2_phi->Fill(jet2_phi, mcEventWeight_Wenu);
                          h_wenu_jet1_eta->Fill(jet1_eta, mcEventWeight_Wenu);
                          h_wenu_jet2_eta->Fill(jet2_eta, mcEventWeight_Wenu);
                          h_wenu_jet1_rap->Fill(jet1_rapidity, mcEventWeight_Wenu);
                          h_wenu_jet2_rap->Fill(jet2_rapidity, mcEventWeight_Wenu);
                          h_wenu_mjj->Fill(mjj, mcEventWeight_Wenu);
                          h_wenu_dPhijj->Fill(DeltaPhi(jet1_phi, jet2_phi), mcEventWeight_Wenu);
                          h_wenu_dRjj->Fill(DeltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Wenu);
                          h_wenu_dPhimetj1->Fill(dPhiJet1Met_noelec, mcEventWeight_Wenu);
                          h_wenu_dPhimetj2->Fill(dPhiJet2Met_noelec, mcEventWeight_Wenu);
                          h_wenu_dPhiMinmetjet->Fill(dPhiMinjetmet_noelec, mcEventWeight_Wenu);
                          // For jet3
                          if (m_signalJet->size() > 2){
                            h_wenu_jet3_pt->Fill(jet3_pt, mcEventWeight_Wenu);
                            h_wenu_jet3_phi->Fill(jet3_phi, mcEventWeight_Wenu);
                            h_wenu_jet3_eta->Fill(jet3_eta, mcEventWeight_Wenu);
                            h_wenu_jet3_rap->Fill(jet3_rapidity, mcEventWeight_Wenu);
                            h_wenu_dPhimetj3->Fill(dPhiJet3Met_noelec, mcEventWeight_Wenu);
                          }
                          // Leptons
                          h_wenu_electron_pt->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Wenu);
                          h_wenu_electron_phi->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Wenu);
                          h_wenu_electron_eta->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Wenu);
                          h_wenu_mT->Fill(mT_electron, mcEventWeight_Wenu);
                          */

                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }




    //-------------------------------
    // Z -> mumu + JET EVENT Cutflow
    //-------------------------------

    if (m_isZmumu && m_isEmilyCutflow){

      if ( (m_signalJet->size() > 0 && monojet_pt > 100.) || (m_signalJet->size() > 1 && jet1_pt > 55. && jet2_pt > 45.) ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Skim cuts");
        if (m_goodMuon->size() > 1) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]At least Two Muon");
          if (pass_OSmuon) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Opposite sign charge");
            if (pass_dimuonPtCut) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Dimuon pT cut");
              if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]MET Trigger");
                if (emulMET_nomu > m_metCut) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]MET cut");
                  if (mll_muon > 66. && mll_muon < 116.) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Zmass window");
                    if (m_goodElectron->size() == 0 && m_goodTau->size() == 0 && m_goodMuon->size() == 2) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Lepton Veto");
                      ////////////////////////
                      // MonoJet phasespace //
                      ////////////////////////
                      if (pass_monoJet && pass_dPhijetmet_nomu) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Monojet cut");
                      }
                      ////////////////////
                      // VBF phasespace //
                      ////////////////////
                      if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_nomu) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]VBF cut");
                      }
                    } // Lepton Veto
                  } // Zmass window
                } // MET cut
              } // MET trigger
            } // Dimuon pT cut
          } // Opposite sign change cut
        } // At least 2 muons
      } // Skim cuts

    } // end cutflow





    //-------------------------------
    // Z -> ee + JET EVENT Cutflow
    //-------------------------------

    if (m_isZee && m_isEmilyCutflow){

      if ( (m_signalJet->size() > 0 && monojet_pt > 100.) || (m_signalJet->size() > 1 && jet1_pt > 55. && jet2_pt > 45.) ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Skim cuts");
        if (m_goodElectron->size() > 1) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]At least Two Electron");
          if (pass_OSelectron) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Opposite sign charge");
            if (pass_dielectronPtCut) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Dielectron pT cut");
              if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Electron Trigger");
                if (emulMET_noelec > m_metCut) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]MET cut");
                  if (mll_electron > 66. && mll_electron < 116.) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Zmass window");
                    if (m_goodMuon->size() == 0 && m_goodTau->size() == 0 && m_goodElectron->size() == 2) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Lepton Veto");
                      ////////////////////////
                      // MonoJet phasespace //
                      ////////////////////////
                      if (pass_monoJet && pass_dPhijetmet_noelec) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Monojet cut");
                      }
                      ////////////////////
                      // VBF phasespace //
                      ////////////////////
                      if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_noelec) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]VBF cut");
                      }
                    } // Lepton Veto
                  } // Zmass window
                } // MET cut
              } // Electron trigger
            } // Dielectron pT cut
          } // Opposite sign change cut
        } // At least 2 muons
      } // Skim cuts

    } // end Cutflow





    //////////////////////////
    // Delete copy containers
    //////////////////////////

    // Deep copies. Clearing containers deletes contents including AuxStore.
    delete m_goodJet;

    // MET study
    //delete m_signalMuon;
    //delete m_signalElectron;

    // VBF study
    delete m_goodMuon;
    delete m_goodElectron;
    delete m_goodTau;
    delete m_signalJet;


    //////////////////////////////////
    // Delete shallow copy containers
    //////////////////////////////////

    // The containers created by the shallow copy are owned by you. Remember to delete them
    delete m_met;
    delete m_metAux;

    delete m_emulmet_nomu;
    delete m_emulmetAux_nomu;

    delete m_emulmet_noelec;
    delete m_emulmetAux_noelec;

    delete muons_shallowCopy.first;
    delete muons_shallowCopy.second;

    delete elec_shallowCopy.first;
    delete elec_shallowCopy.second;

    delete phot_shallowCopy.first;
    delete phot_shallowCopy.second;

    delete tau_shallowCopy.first;
    delete tau_shallowCopy.second;

    delete jet_shallowCopy.first;
    delete jet_shallowCopy.second;

  } // end for loop over systematics

  return EL::StatusCode::SUCCESS;
}



  EL::StatusCode ZinvxAODAnalysis :: postExecute ()
  {
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ZinvxAODAnalysis :: finalize ()
  {
    // This method is the mirror image of initialize(), meaning it gets
    // called after the last event has been processed on the worker node
    // and allows you to finish up any objects you created in
    // initialize() before they are written to disk.  This is actually
    // fairly rare, since this happens separately for each worker node.
    // Most of the time you want to do your post-processing on the
    // submission node after all your histogram outputs have been
    // merged.  This is different from histFinalize() in that it only
    // gets called on worker nodes that processed input events.

    // cutflow
    if (m_useBitsetCutflow) m_BitsetCutflow->PushBitSet();;

    //*************************
    // deleting of all tools
    // ************************

    // GRL
    if (m_grl) {
      delete m_grl;
      m_grl = 0;
    }

    // cleaning up trigger tools
    if( m_trigConfigTool ) {
      delete m_trigConfigTool;
      m_trigConfigTool = 0;
    }
    if( m_trigDecisionTool ) {
      delete m_trigDecisionTool;
      m_trigDecisionTool = 0;
    }


    // in finalize, delete tool:
    /// Muon Calibration and Smearing Tool
    if(m_muonCalibrationAndSmearingTool){
      delete m_muonCalibrationAndSmearingTool;
      m_muonCalibrationAndSmearingTool = 0;
    }

    /// Electron Calibration and Smearing Tool
    if(m_egammaCalibrationAndSmearingTool){
      delete m_egammaCalibrationAndSmearingTool;
      m_egammaCalibrationAndSmearingTool = 0;
    }

    /// Muon selector tool
    if(m_muonSelection){
      delete m_muonSelection;
      m_muonSelection = 0;
    }

    /// Loose Muon selector tool
    if(m_loosemuonSelection){
    delete m_loosemuonSelection;
    m_loosemuonSelection = 0;
    }

    /// Electron selector tool
    if(m_LHToolTight2015){
      delete m_LHToolTight2015;
      m_LHToolTight2015 = 0;
    }
    if(m_LHToolMedium2015){
      delete m_LHToolMedium2015;
      m_LHToolMedium2015 = 0;
    }
    if(m_LHToolLoose2015){
      delete m_LHToolLoose2015;
      m_LHToolLoose2015 = 0;
    }


    /// Recomputing the photon ID flags
    if(m_photonTightIsEMSelector){
      delete m_photonTightIsEMSelector;
      m_photonTightIsEMSelector = 0;
    }
    if(m_photonMediumIsEMSelector){
      delete m_photonMediumIsEMSelector;
      m_photonMediumIsEMSelector = 0;
    }
    if(m_photonLooseIsEMSelector){
      delete m_photonLooseIsEMSelector;
      m_photonLooseIsEMSelector = 0;
    }

    /// IsolationSelectionTool
    if(m_IsolationSelectionTool){
      delete m_IsolationSelectionTool;
      m_IsolationSelectionTool = 0;
    }

    /// IsolationSelectionTool for VBF signal
    if(m_IsoToolVBF){
      delete m_IsoToolVBF;
      m_IsoToolVBF = 0;
    }

    /// Tau Smearing Tool
    if(m_tauSmearingTool){
      delete m_tauSmearingTool;
      m_tauSmearingTool = 0;
    }

    /// Tau Selection Tool
    if(m_tauSelTool){
      delete m_tauSelTool;
      m_tauSelTool = 0;
    }

    /// Tau Selection Tool for VBF signal
    if(m_tauSelToolVBF){
      delete m_tauSelToolVBF;
      m_tauSelToolVBF = 0;
    }

    /// JES Calibration
    if(m_jetCalibration){
      delete m_jetCalibration;
      m_jetCalibration = 0;
    }

    /// JES uncertainty
    if(m_jetUncertaintiesTool){
      delete m_jetUncertaintiesTool;
      m_jetUncertaintiesTool = 0;
    }

    /// JER Tool
    if(m_jerTool){
      delete m_jerTool;
      m_jerTool = 0;
    }

    ///  JER Smearing Tool
    if(m_jerSmearingTool){
      delete m_jerSmearingTool;
      m_jerSmearingTool = 0;
    }

    /// JVT Tool
    if(m_jvtag){
      delete m_jvtag;
      m_jvtag = 0;
    }

    /// Jet Cleaning Tool
    if(m_jetCleaningTight) {
      delete m_jetCleaningTight;
      m_jetCleaningTight = 0;
    }
    if(m_jetCleaningLoose) {
      delete m_jetCleaningLoose;
      m_jetCleaningLoose = 0;
    }

    /// MET Tool
    if(m_metMaker){
      delete m_metMaker;
      m_metMaker = 0;
    }

    /// Muon Efficiency Tool
    if(m_muonEfficiencySFTool){
      delete m_muonEfficiencySFTool;
      m_muonEfficiencySFTool = 0;
    }

    /// Muon Isolation Tool
    if(m_muonIsolationSFTool){
      delete m_muonIsolationSFTool;
      m_muonIsolationSFTool = 0;
    }

    /// Muon TTVA Efficiency Tool
    if(m_muonTTVAEfficiencySFTool){
      delete m_muonTTVAEfficiencySFTool;
      m_muonTTVAEfficiencySFTool = 0;
    }

    /// Muon Trigger Scale Factor Tool
    if(m_muonTriggerSFTool){
      delete m_muonTriggerSFTool;
      m_muonTriggerSFTool = 0;
    }

    /// Electron Efficiency Tool
    if(m_elecEfficiencySFTool_reco){
      delete m_elecEfficiencySFTool_reco;
      m_elecEfficiencySFTool_reco = 0;
    }

    if(m_elecEfficiencySFTool_id_Loose){
      delete m_elecEfficiencySFTool_id_Loose;
      m_elecEfficiencySFTool_id_Loose = 0;
    }

    if(m_elecEfficiencySFTool_id_Medium){
      delete m_elecEfficiencySFTool_id_Medium;
      m_elecEfficiencySFTool_id_Medium = 0;
    }

    if(m_elecEfficiencySFTool_id_Tight){
      delete m_elecEfficiencySFTool_id_Tight;
      m_elecEfficiencySFTool_id_Tight = 0;
    }

    if(m_elecEfficiencySFTool_iso_Loose){
      delete m_elecEfficiencySFTool_iso_Loose;
      m_elecEfficiencySFTool_iso_Loose = 0;
    }

    if(m_elecEfficiencySFTool_iso_Medium){
      delete m_elecEfficiencySFTool_iso_Medium;
      m_elecEfficiencySFTool_iso_Medium = 0;
    }

    if(m_elecEfficiencySFTool_iso_Tight){
      delete m_elecEfficiencySFTool_iso_Tight;
      m_elecEfficiencySFTool_iso_Tight = 0;
    }

    if(m_elecEfficiencySFTool_trigEff){
      delete m_elecEfficiencySFTool_trigEff;
      m_elecEfficiencySFTool_trigEff = 0;
    }

    if(m_elecEfficiencySFTool_trigSF_Loose){
      delete m_elecEfficiencySFTool_trigSF_Loose;
      m_elecEfficiencySFTool_trigSF_Loose = 0;
    }

    /// Jet JVT Efficiency Tool
    if(m_jvtefficiencyTool){
      delete m_jvtefficiencyTool;
      m_jvtefficiencyTool = 0;
    }

    /// Tau Efficiency Tool
    if(m_tauEffTool){
      delete m_tauEffTool;
      m_tauEffTool = 0;
    }

    /// MET Tools
    if(m_metSystTool){
      delete m_metSystTool;
      m_metSystTool = 0;
    }

    /// Isolation Correction Tool
    if(m_isoCorrTool){
      delete m_isoCorrTool;
      m_isoCorrTool = 0;
    }

    /// PileupReweighting Tool
    if(m_prwTool){
      delete m_prwTool;
      m_prwTool = 0;
    }

    /// Cutflow
    if(m_useBitsetCutflow && m_BitsetCutflow){
      delete m_BitsetCutflow;
      m_BitsetCutflow = 0;
    }



    // print out the number of Overlap removal
    Info("finalize()", "======================================================");
    Info("finalize()", "Number overlap electrons:    %i / %i", nOverlapElectrons, nInputElectrons);
    Info("finalize()", "Number overlap muons:    %i / %i", nOverlapMuons, nInputMuons);
    Info("finalize()", "Number overlap jets:    %i / %i", nOverlapJets, nInputJets);
    Info("finalize()", "Number overlap taus:    %i / %i", nOverlapTaus, nInputTaus);
    Info("finalize()", "Number overlap photons:    %i / %i", nOverlapPhotons, nInputPhotons);
    Info("finalize()", "======================================================");

    // print out the final number of clean events
    Info("finalize()", "Number of clean events = %i", m_numCleanEvents);

    // print out Cutflow
    Info("finalize()", "================================================");
    Info("finalize()", "===============  Base Cutflow  =================");
    for(int i=0; i<5; ++i) {
      int j = i+1;
      Info("finalize()", "Event cutflow (%i) = %i", j, m_eventCutflow[i]);
    }
    if (m_isZnunu){
      Info("finalize()", "===============  Znunu Cutflow  ==================");
      for(int i=5; i<15 ; ++i) {
        int j = i+1;
        Info("finalize()", "Znunu Event cutflow (%i) = %i", j, m_eventCutflow[i]);
      }
    }
    if (m_isZmumu){
      Info("finalize()", "===============  Zmumu Cutflow  =================");
      for(int i=16; i<27 ; ++i) {
        int j = i-10;
        Info("finalize()", "Zmumu Event cutflow (%i) = %i", j, m_eventCutflow[i]);
      }
    }
    if (m_isZee){
      Info("finalize()", "===============  Zee Cutflow  ==================");
      for(int i=28; i<39 ; ++i) {
        int j = i-22;
        Info("finalize()", "Zee Event cutflow (%i) = %i", j, m_eventCutflow[i]);
      }
    }


    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ZinvxAODAnalysis :: histFinalize ()
  {
    // This method is the mirror image of histInitialize(), meaning it
    // gets called after the last event has been processed on the worker
    // node and allows you to finish up any objects you created in
    // histInitialize() before they are written to disk.  This is
    // actually fairly rare, since this happens separately for each
    // worker node.  Most of the time you want to do your
    // post-processing on the submission node after all your histogram
    // outputs have been merged.  This is different from finalize() in
    // that it gets called on all worker nodes regardless of whether
    // they processed input events.
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ZinvxAODAnalysis :: addHist(std::map<std::string, TH1*> &hMap, std::string tag,
      int bins, double min, double max) {

    std::string label = tag;
    TH1* h_temp = new TH1F(label.c_str(), "", bins, min, max);
    wk()->addOutput (h_temp);
    hMap[label] = h_temp;

    return EL::StatusCode::SUCCESS;
  }

  EL::StatusCode ZinvxAODAnalysis :: addHist(std::map<std::string, TH1*> &hMap, std::string tag,
      int bins, Float_t binArray[]) {

    std::string label = tag;
    TH1* h_temp = new TH1F(label.c_str(), "", bins, binArray);
    wk()->addOutput (h_temp);
    hMap[label] = h_temp;

    return EL::StatusCode::SUCCESS;
  }


  EL::StatusCode ZinvxAODAnalysis :: passMuonSelection(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(mu) = false;
    selectDec(mu) = false; // To select objects for Overlap removal
    dec_signal(mu) = false;

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double muPt = (mu.pt()) * 0.001; /// GeV
    //if ( muPt < 4. ) return EL::StatusCode::SUCCESS;

    // Muon Calibration
    if (!m_isData){
      if(m_muonCalibrationAndSmearingTool->applyCorrection(mu) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
        Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // MuonSelectionTool(Medium)
    if(!m_muonSelection->accept(mu)) return EL::StatusCode::SUCCESS;
    // MuonSelectionTool (Loose)
    //if(!m_loosemuonSelection->accept(mu)) return EL::StatusCode::SUCCESS;

    // Muon tranverse momentum cut
    if (muPt <= 10. ) return EL::StatusCode::SUCCESS; /// veto muon

    // Muon eta cut
    double muEta = mu.eta();
    if (fabs(muEta) >= 2.47) return EL::StatusCode::SUCCESS;

    //  if (mu.muonType()=xAOD::Muon_v1::Combined) return EL::StatusCode::SUCCESS;
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return EL::StatusCode::SUCCESS;


    // Baseline Muon
    dec_baseline(mu) = true;
    selectDec(mu) = true; // To select objects for Overlap removal


    // Muon pt cut
    if (muPt <= 25. ) return EL::StatusCode::SUCCESS;
    // Muon eta cut
    if (fabs(muEta) >= 2.4) return EL::StatusCode::SUCCESS;

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
      tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 3.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(mu)) return EL::StatusCode::SUCCESS;

    // Signal Muon
    dec_signal(mu) = true;

    return EL::StatusCode::SUCCESS;

  }





  EL::StatusCode ZinvxAODAnalysis :: passMuonVBF(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(mu) = false;
    selectDec(mu) = false; // To select objects for Overlap removal

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSignal. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double muPt = (mu.pt()) * 0.001; /// GeV
    //if ( muPt < 4. ) return EL::StatusCode::SUCCESS;

    // Muon Calibration
    if (!m_isData){
      if(m_muonCalibrationAndSmearingTool->applyCorrection(mu) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
        Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // MuonSelectionTool (Loose)
    if(!m_loosemuonSelection->accept(mu)) return EL::StatusCode::SUCCESS;

    // Muon tranverse momentum
    if (muPt < m_muonPtCut ) return EL::StatusCode::SUCCESS;

    // Muon eta cut
    double muEta = mu.eta();
    if (fabs(muEta) > m_muonEtaCut) return EL::StatusCode::SUCCESS;

    // Combined (CB) or Segment-tagged (ST) muons (excluding Stand-alone (SA), Calorimeter-tagged (CaloTag) muons etc..)
    //if (!(mu.muonType() == xAOD::Muon::Combined || mu.muonType() == xAOD::Muon::SegmentTagged)) return EL::StatusCode::SUCCESS;
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return EL::StatusCode::SUCCESS;

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
      tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 3.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsoToolVBF->accept(mu)) return EL::StatusCode::SUCCESS;
    // Isolation for specific muon pT range
    //m_isoPtCut = true; if (muPt > m_isoMuonPtMin && muPt < m_isoMuonPtMax && !m_IsoToolVBF->accept(mu)) return EL::StatusCode::SUCCESS;


    dec_baseline(mu) = true;
    selectDec(mu) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }




  EL::StatusCode ZinvxAODAnalysis :: passElectronSelection(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(elec) = false;
    selectDec(elec) = false; // To select objects for Overlap removal
    dec_signal(elec) = false;

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_identification

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double elecPt = (elec.pt()) * 0.001; /// GeV
    //if ( elecPt < 4. ) return EL::StatusCode::SUCCESS;

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return EL::StatusCode::SUCCESS;

    // "Please apply the identification to uncalibrated electron object. ID scale factors are to be applied to calibrated objects."
    // LH Electron (Medium)
    bool LHmediumSel = false;
    LHmediumSel = m_LHToolMedium2015->accept(elec);
    if (!LHmediumSel) return EL::StatusCode::SUCCESS;
    // LH Electron (Loose)
    //bool LHlooseSel = false;
    //LHlooseSel = m_LHToolLoose2015->accept(elec);
    //if (!LHlooseSel) return EL::StatusCode::SUCCESS;

    //Info("execute()", "  Selected electron pt from new Electron Container = %.2f GeV", (elec.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(elec) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original electron values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = elec.caloCluster()->eta();
    double Eta = elec.caloCluster()->etaBE(2);
    if ( fabs(Eta) >= 2.47 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

    // pT cut
    double elecPtCut = 10.0; /// GeV
    if (elecPt <= elecPtCut) return EL::StatusCode::SUCCESS; /// veto electron

    // Baseline Electron
    dec_baseline(elec) = true;
    selectDec(elec) = true; // To select objects for Overlap removal


    // pT cut
    if (elecPt <= 25.) return EL::StatusCode::SUCCESS; /// veto electron

    // d0 / z0 cuts applied
    // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle *tp = elec.trackParticle() ; //your input track particle from the electron
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 5.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( elec.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(elec)) return EL::StatusCode::SUCCESS;

    // Signal Electron
    dec_signal(elec) = true;

    return EL::StatusCode::SUCCESS;

  }




  EL::StatusCode ZinvxAODAnalysis :: passElectronVBF(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo,
      xAOD::Vertex* primVertex){

    dec_baseline(elec) = false;
    selectDec(elec) = false; // To select objects for Overlap removal

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double elecPt = (elec.pt()) * 0.001; /// GeV
    //if ( elecPt < 4. ) return EL::StatusCode::SUCCESS;

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return EL::StatusCode::SUCCESS;

    // "Please apply the identification to uncalibrated electron object. ID scale factors are to be applied to calibrated objects."
    // LH Electron identification
    //
    // LH Electron (Loose)
    bool LHlooseSel = false;
    LHlooseSel = m_LHToolLoose2015->accept(elec);
    if (!LHlooseSel) return EL::StatusCode::SUCCESS;
    /*
    // LH Electron (Medium)
    bool LHmediumSel = false;
    LHmediumSel = m_LHToolMedium2015->accept(elec);
    if (!LHmediumSel) return EL::StatusCode::SUCCESS;
    // LH Electron (Tight)
    bool LHtightSel = false;
    LHtightSel = m_LHToolTight2015->accept(elec);
    if (!LHtightSel) return EL::StatusCode::SUCCESS;
    */

    //Info("execute()", "  Selected electron pt from new Electron Container = %.2f GeV", (elec.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(elec) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original electron values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = elec.caloCluster()->eta();
    double Eta = elec.caloCluster()->etaBE(2);
    if ( fabs(Eta) > m_elecEtaCut || (fabs(Eta) > 1.37 && fabs(Eta) < 1.52)) return EL::StatusCode::SUCCESS;
    //if ( fabs(Eta) > m_elecEtaCut ) return EL::StatusCode::SUCCESS;

    /// pT cut
    if (elecPt < m_elecPtCut ) return EL::StatusCode::SUCCESS; /// veto electron

    // d0 / z0 cuts applied
    // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle *tp = elec.trackParticle() ; //your input track particle from the electron
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 5.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( elec.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation Correction Tool
    if(m_isoCorrTool->applyCorrection(elec) == CP::CorrectionCode::Error) { 
      Error("execute()", "IsolationCorrectionTool returns Error CorrectionCode");
    }

    // Isolation requirement
    if (!m_IsoToolVBF->accept(elec)) return EL::StatusCode::SUCCESS;

    dec_baseline(elec) = true;
    selectDec(elec) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }





  EL::StatusCode ZinvxAODAnalysis :: passPhotonSelection(xAOD::Photon& phot,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(phot) = false;
    selectDec(phot) = false; // To select objects for Overlap removal

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passPhotonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Photon author cuts
    if ( !(phot.author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) )
      return EL::StatusCode::SUCCESS;


    //Info("execute()", "  Selected photon pt from new Photon Container = %.2f GeV", (phot.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = phot.caloCluster()->eta();
    double Eta = phot.caloCluster()->etaBE(2);
    if ( fabs(Eta) >= 2.37 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

    // pT cut
    double photPt = (phot.pt()) * 0.001; /// GeV
    double photPtCut = 10.0; /// GeV
    if (photPt <= photPtCut) return EL::StatusCode::SUCCESS; /// veto photon

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return EL::StatusCode::SUCCESS;

    // MC fudge tool
    if (!m_isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    //if (!m_IsolationSelectionTool->accept(phot)) return EL::StatusCode::SUCCESS;


    dec_baseline(phot) = true;
    selectDec(phot) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }




  EL::StatusCode ZinvxAODAnalysis :: passPhotonVBF(xAOD::Photon& phot,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(phot) = false;
    selectDec(phot) = false; // To select objects for Overlap removal

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passPhotonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Photon author cuts
    if ( !(phot.author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) )
      return EL::StatusCode::SUCCESS;


    //Info("execute()", "  Selected photon pt from new Photon Container = %.2f GeV", (phot.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = phot.caloCluster()->eta();
    double Eta = phot.caloCluster()->etaBE(2);
    //if ( fabs(Eta) >= m_photEtaCut || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;
    if ( fabs(Eta) > m_photEtaCut ) return EL::StatusCode::SUCCESS;

    // pT cut
    double photPt = (phot.pt()) * 0.001; /// GeV
    if (photPt < m_photPtCut) return EL::StatusCode::SUCCESS; /// veto photon

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return EL::StatusCode::SUCCESS;

    // MC fudge tool
    if (!m_isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;
    //if (!m_photonLooseIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    //if (!m_IsoToolVBF->accept(phot)) return EL::StatusCode::SUCCESS;


    dec_baseline(phot) = true;
    selectDec(phot) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }



  EL::StatusCode ZinvxAODAnalysis :: passTauSelection(xAOD::TauJet& tau,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(tau) = false;
    selectDec(tau) = false; // To select objects for Overlap removal

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Tau Smearing (for MC)
    if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !m_isData){ // it's MC!
      if(m_tauSmearingTool->applyCorrection(tau) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original tau values are taken.
        Error("execute()", "TauSmearingTool returns Error CorrectionCode");
      }
    }

    //Info("execute()", "  original tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    // TauSelectionTool (Medium)
    if(!m_tauSelTool->accept(tau)) return EL::StatusCode::SUCCESS;
    // TauSelectionTool (Loose for VBF)
    //if(!m_tauSelToolVBF->accept(tau)) return EL::StatusCode::SUCCESS;

    //Info("execute()", "  Selected tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    dec_baseline(tau) = true;
    selectDec(tau) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }


  EL::StatusCode ZinvxAODAnalysis :: passTauVBF(xAOD::TauJet& tau,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(tau) = false;
    selectDec(tau) = false; // To select objects for Overlap removal

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // TauOverlappingElectronLLHDecorator
    m_tauOverlappingElectronLLHDecorator->decorate(tau);

    // Tau Smearing (for MC)
    if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !m_isData){ // it's MC!
      if(m_tauSmearingTool->applyCorrection(tau) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original tau values are taken.
        Error("execute()", "TauSmearingTool returns Error CorrectionCode");
      }
    }

    //Info("execute()", "  original tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    // TauSelectionTool (Loose for VBF)
    if(!m_tauSelToolVBF->accept(tau)) return EL::StatusCode::SUCCESS;

    //Info("execute()", "  Selected tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    dec_baseline(tau) = true;
    selectDec(tau) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }




  bool ZinvxAODAnalysis :: IsBadJet(xAOD::Jet& jet) {

    if (m_doORtool && overlapAcc(jet)) return false;

    double jetPt = (jet.pt()) * 0.001; /// GeV

    //Info("execute()", "  corrected jet pt in IsBadJet function = %.2f GeV", jetPt );
    //Info("execute()", "  updated jet jvt in IsBadJet function = %.2f ", cacc_jvt(jet) );

    // Pile-up
    if ( cacc_jvt(jet) < 0.59 && fabs(jet.eta()) < 2.4 && jetPt < 50.0 ) return false;

    // pT cut
    if ( jetPt < m_jetPtCut || fabs(jet.eta()) > m_jetEtaCut) return false; 

    // Jet Cleaning Tool
    dec_bad(jet) = !m_jetCleaningLoose->accept( jet );

    return dec_bad(jet);

  }


  bool ZinvxAODAnalysis :: IsSignalJet(xAOD::Jet& jet) {

    if (m_doORtool && overlapAcc(jet)) return false;
    if (!dec_baseline(jet)) return false;

    double jetPt = (jet.pt()) * 0.001; /// GeV

    // pT, eta cut
    if ( jetPt < m_jetPtCut || fabs(jet.eta()) > m_jetEtaCut ) return false;

    bool isgoodjet = !dec_bad(jet) && (cacc_jvt(jet) > 0.59 || fabs(jet.eta()) > 2.4 || jetPt > 50.0);

    dec_signal(jet) = isgoodjet;

    return isgoodjet;

  }


  float ZinvxAODAnalysis :: GetGoodMuonSF(xAOD::Muon& mu,
      const bool recoSF, const bool isoSF, const bool ttvaSF) {

    float sf(1.);

    if (recoSF) {
      float sf_reco(1.);
      if (m_muonEfficiencySFTool->getEfficiencyScaleFactor( mu, sf_reco ) == CP::CorrectionCode::OutOfValidityRange) {
        Error("execute()", " GetGoodMuonSF: Reco getEfficiencyScaleFactor out of validity range");
      }
      //Info("execute()", "  GetGoodMuonSF: sf_reco = %.5f ", sf_reco );
      sf *= sf_reco;
    }

    if (isoSF) {
      float sf_iso(1.);
      if (m_muonIsolationSFTool->getEfficiencyScaleFactor( mu, sf_iso ) == CP::CorrectionCode::OutOfValidityRange) {
        Error("execute()", " GetGoodMuonSF: Iso getEfficiencyScaleFactor out of validity range");
      }
      //Info("execute()", "  GetGoodMuonSF: sf_iso = %.5f ", sf_iso );
      sf *= sf_iso;
    }

    if (ttvaSF) {
      float sf_TTVA(1.);
      if (m_muonTTVAEfficiencySFTool->getEfficiencyScaleFactor( mu, sf_TTVA ) == CP::CorrectionCode::OutOfValidityRange) {
        Error("execute()", " GetGoodMuonSF: TTVA getEfficiencyScaleFactor out of validity range");
      }
      //Info("execute()", "  GetGoodMuonSF: sf_TTVA = %.5f ", sf_TTVA );
      sf *= sf_TTVA;
    }

    //Info("execute()", "  GetGoodMuonSF: Good Muon SF = %.5f ", sf );
    dec_scalefactor(mu) = sf;
    return sf;

  }

  double ZinvxAODAnalysis :: GetTotalMuonSF(xAOD::MuonContainer& muons,
      bool recoSF, bool isoSF, bool ttvaSF) {

    double sf(1.);
    bool iso_SF = isoSF;

    for (const auto& muon : muons) {
      iso_SF = true;
      double muPt = muon->pt() * 0.001; /// GeV
      if (m_isoPtCut && (muPt < m_isoMuonPtMin || muPt > m_isoMuonPtMax)) {
        iso_SF = false;
      }
      //Info("execute()", "  GetTotalMuonSF: Muon pT = %.2f GeV ", muPt );
      sf *= GetGoodMuonSF(*muon, recoSF, iso_SF, ttvaSF);
    }

    //Info("execute()", "  GetTotalMuonSF: Total Muon SF = %.5f ", sf );
    return sf;

  }


  float ZinvxAODAnalysis :: GetGoodElectronSF(xAOD::Electron& elec,
      const bool recoSF, const bool idSF, const bool isoSF, const bool trigSF) {

    float sf(1.);

    if (recoSF) {
      double sf_reco(1.);
      if (m_elecEfficiencySFTool_reco->getEfficiencyScaleFactor( elec, sf_reco ) == CP::CorrectionCode::Ok) {
        sf *= sf_reco;
        //Info("execute()", "  GetGoodElectronSF: sf_reco = %.5f ", sf_reco );
      }
      else {
        Error("execute()", " GetGoodElectronSF: Reco getEfficiencyScaleFactor returns Error CorrectionCode");
      }
    }

    if (idSF) {
      double sf_id(1.);
      if (m_elecEfficiencySFTool_id_Loose->getEfficiencyScaleFactor( elec, sf_id ) == CP::CorrectionCode::Ok) {
        sf *= sf_id;
        //Info("execute()", "  GetGoodElectronSF: sf_id = %.5f ", sf_id );
      }
      else {
        Error("execute()", " GetGoodElectronSF: Id (Loose) getEfficiencyScaleFactor returns Error CorrectionCode");
      }
    }

    if (isoSF) {
      double sf_iso(1.);
      if (m_elecEfficiencySFTool_iso_Loose->getEfficiencyScaleFactor( elec, sf_iso ) == CP::CorrectionCode::Ok) {
        sf *= sf_iso;
        //Info("execute()", "  GetGoodElectronSF: sf_iso = %.5f ", sf_iso );
      }
      else {
        Error("execute()", " GetGoodElectronSF: Iso (Loose) getEfficiencyScaleFactor returns Error CorrectionCode");
      }
    }

    if (trigSF) {
      double sf_trig(1.);
      if (m_elecEfficiencySFTool_trigSF_Loose->getEfficiencyScaleFactor( elec, sf_trig ) == CP::CorrectionCode::Ok) {
        sf *= sf_trig;
        //Info("execute()", "  GetGoodElectronSF: sf_trig = %.5f ", sf_trig );
      }
      else {
        Error("execute()", " GetGoodElectronSF: Trigger (Loose) getEfficiencyScaleFactor returns Error CorrectionCode");
      }
    }


    //Info("execute()", "  GetGoodElectronSF: Good Electron SF = %.5f ", sf );
    dec_scalefactor(elec) = sf;
    return sf;

  }


  float ZinvxAODAnalysis :: GetTotalElectronSF(xAOD::ElectronContainer& electrons,
      bool recoSF, bool idSF, bool isoSF, bool trigSF) {

    float sf(1.);
    bool trig_SF = trigSF;

    for (const auto& electron : electrons) {
      trig_SF = true;
      if (electrons.at(0) != electron) trig_SF = false; // Apply trigSF to only leading electron
      sf *= GetGoodElectronSF(*electron, recoSF, idSF, isoSF, trig_SF);
      //Info("execute()", "  Good electron pt = %.2f GeV, trigSF = %d, Electron SF = %.2f", (electron->pt() * 0.001), trigSF, sf);
    }

    //Info("execute()", "  GetTotalElectronSF: Total Electron SF = %.5f ", sf );
    return sf;

  }



  int ZinvxAODAnalysis :: NumIsoTracks(const xAOD::TrackParticleContainer* inTracks,
      xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
    //
    //  Fill track objects with information about isolated tracks. For being isolated, there should be no track
    //  above 3 GeV satisfying quality requirement that are within a cone of 0.4 around the probed track.
    //
    //============================================================================================================

    // Integer to return with this function
    // ------------------------------------

    int NisoTrack = 0;


    // Loop over tracks in the event
    // -----------------------------

    for( auto trk_itr : *inTracks ){

      int NCloseby = 0;

      // Apply quality cuts on the track with Pt>10 GeV
      // ----------------------------------------------
      float pt   = (trk_itr->pt()) * 0.001; /// GeV
      float eta  = trk_itr->eta();
      //float phi  = trk_itr->phi();
      float d0   = trk_itr->d0();
      float z0   = (trk_itr->z0() + trk_itr->vz() - primVertex->z());
      int Ndof    = trk_itr->numberDoF();
      if (Ndof==0) Ndof = 1;

      float Chi2 = trk_itr->chiSquared()/Ndof;
      uint8_t nSCT      = -1;
      uint8_t nPix      = -1;
      if(!trk_itr->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
      if(!trk_itr->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
      uint8_t NHits   = nSCT + nPix;

      if(pt < Pt_High || fabs(eta) > 2.5 || fabs(z0) >= 2.0 || fabs(d0) >= 1.0 || Chi2 >= 3.0 || NHits < 5) continue;

      // Loop over *other* tracks and apply quality cuts too
      // ---------------------------------------------------

      for( auto trk2_itr : *inTracks ){

        if ( trk2_itr == trk_itr ) continue;

        float pt_2   = (trk2_itr->pt()) * 0.001; /// GeV
        float eta_2  = trk2_itr->eta();
        //float phi_2  = trk2_itr->phi();
        float d0_2   = trk2_itr->d0();
        float z0_2   = (trk2_itr->z0() + trk2_itr->vz() - primVertex->z());
        int Ndof_2    = trk2_itr->numberDoF();
        if (Ndof_2==0) Ndof_2 = 1;

        float Chi2_2 = trk2_itr->chiSquared()/Ndof_2;
        uint8_t nSCT_2      = -1;
        uint8_t nPix_2      = -1;
        if(!trk2_itr->summaryValue(nPix_2,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
        if(!trk2_itr->summaryValue(nSCT_2,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
        uint8_t NHits_2   = nSCT_2 + nPix_2;

        if(pt_2 < Pt_Low || fabs(eta_2) > 2.5 || fabs(z0_2) >= 2.0 || fabs(d0_2) >= 1.0 || Chi2_2 >= 3.0 || NHits_2 < 5) continue;

        float dR = trk2_itr->p4().DeltaR(trk_itr->p4());
        //Info("execute()", "  dR(track1 with pt=%f, track2 with pt=%f) = %f ", pt, pt_2, dR);

        // Count the number of track in a 0.4 cone around the probed track
        // ---------------------------------------------------------------

        if(dR < 0.4) NCloseby++;
      } // end 2nd loop

      // Fill the vertex collection
      // --------------------------

      // Note: The object properties are: E, pT, eta, phi, d0, d0_err, d0_vx, d0_err_vx, charge, n_SCTHits, n_PixelHits

      if (NCloseby < 1) NisoTrack++;

    } // end 1st loop

    return NisoTrack;
  }


  int ZinvxAODAnalysis :: NumMuonIsoTrack(xAOD::MuonContainer* muons, const xAOD::TrackParticleContainer* inTracks,
      xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
    //
    //  Apply the same criteria as in the SetIsoTracks function to determine how many muons are isolated
    //  according to this definition.
    //
    //============================================================================================================

    // Integer to return with this function
    // ------------------------------------

    int NisoMuon = 0;


    // Loop over muon objects obtained in the event
    // --------------------------------------------

    for (const auto& muon_itr : *muons) { // C++11 shortcut

      int NCloseby = 0;

      // Apply muon cuts on the baseline muon with pT>10GeV
      float muon_pt = (muon_itr->pt()) * 0.001; /// GeV

      //if (!dec_baseline(muon_itr) || !m_IsolationSelectionTool->accept(muon_itr)) continue; 
      if (!dec_baseline(*muon_itr) || muon_pt < Pt_High) continue; 

      // Loop over *other* tracks and apply quality cuts too
      // ---------------------------------------------------

      for( auto trk_itr : *inTracks ){

        // Apply quality cuts on the track with Pt>3 GeV
        // ----------------------------------------------
        float pt   = (trk_itr->pt()) * 0.001; /// GeV
        float eta  = trk_itr->eta();
        //float phi  = trk_itr->phi();
        float d0   = trk_itr->d0();
        float z0   = (trk_itr->z0() + trk_itr->vz() - primVertex->z());
        int Ndof    = trk_itr->numberDoF();
        if (Ndof==0) Ndof = 1;

        float Chi2 = trk_itr->chiSquared()/Ndof;
        uint8_t nSCT      = -1;
        uint8_t nPix      = -1;
        if(!trk_itr->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
        if(!trk_itr->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
        uint8_t NHits   = nSCT + nPix;

        if(pt < Pt_Low || fabs(eta) > 2.5 || fabs(z0) >= 2.0 || fabs(d0) >= 1.0 || Chi2 >= 3.0 || NHits < 5) continue;


        float dR = trk_itr->p4().DeltaR(muon_itr->p4());
        //Info("execute()", "  dR(Muon with pt=%f, track with pt=%f) = %f ", muon_pt, pt, dR);

        // Count the number of track in a 0.4 cone around the probed track
        // ---------------------------------------------------------------

        if(dR < 0.4) NCloseby++;
      } // end 2nd loop

      // Count the number of isolated muon in the event following track iso criteria
      // ---------------------------------------------------------------------------
      //
      // Note 1: We assume here that the track in the collection which is the closest to the muon is actually the
      //         one that come from the muon. This will happen for each muon, so the muon will be isolated if there
      //         other tracks except this clostest one which is the muon itself (hence NCloseby<=1). It is also
      //         possible that no tracks are found in the cone because of the Chi2 cut on the muon might not exactly
      //         match the one required for the closeby tracks.
      //
      // Note 2: Typically, there will be one Iso muon per W events and two for Z.

      if (NCloseby <= 1) NisoMuon++;

    } // end 1st loop

    return NisoMuon;
  }


  int ZinvxAODAnalysis :: NumElecIsoTrack(xAOD::ElectronContainer* electrons, const xAOD::TrackParticleContainer* inTracks,
      xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
    //
    //  Apply the same criteria as in the SetIsoTracks function to determine how many electrons are isolated
    //  according to this definition.
    //
    //============================================================================================================

    // Integer to return with this function
    // ------------------------------------

    int NisoElec = 0;


    // Loop over electron objects obtained in the event
    // --------------------------------------------

    for (const auto& elec_itr : *electrons) { // C++11 shortcut

      int NCloseby = 0;

      // Apply electron cuts on the baseline electron with pT>10GeV
      float elec_pt = (elec_itr->pt()) * 0.001; /// GeV

      //if (!dec_baseline(elec_itr) || !m_IsolationSelectionTool->accept(elec_itr)) continue; 
      if (!dec_baseline(*elec_itr) || elec_pt < Pt_High) continue; 

      // Loop over *other* tracks and apply quality cuts too
      // ---------------------------------------------------

      for( auto trk_itr : *inTracks ){

        // Apply quality cuts on the track with Pt>3 GeV
        // ----------------------------------------------
        float pt   = (trk_itr->pt()) * 0.001; /// GeV
        float eta  = trk_itr->eta();
        //float phi  = trk_itr->phi();
        float d0   = trk_itr->d0();
        float z0   = (trk_itr->z0() + trk_itr->vz() - primVertex->z());
        int Ndof    = trk_itr->numberDoF();
        if (Ndof==0) Ndof = 1;

        float Chi2 = trk_itr->chiSquared()/Ndof;
        uint8_t nSCT      = -1;
        uint8_t nPix      = -1;
        if(!trk_itr->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
        if(!trk_itr->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
        uint8_t NHits   = nSCT + nPix;

        if(pt < Pt_Low || fabs(eta) > 2.5 || fabs(z0) >= 2.0 || fabs(d0) >= 1.0 || Chi2 >= 3.0 || NHits < 5) continue;


        float dR = trk_itr->p4().DeltaR(elec_itr->p4());
        //Info("execute()", "  dR(Electron with pt=%f, track with pt=%f) = %f ", elec_pt, pt, dR);

        // Count the number of track in a 0.4 cone around the probed track
        // ---------------------------------------------------------------

        if(dR < 0.4) NCloseby++;
      } // end 2nd loop

      // Count the number of isolated elec in the event following track iso criteria
      // ---------------------------------------------------------------------------

      if (NCloseby <= 1) NisoElec++;

    } // end 1st loop

    return NisoElec;
  }



  float ZinvxAODAnalysis :: DeltaPhi(float phi1, float phi2) {

    float dPhi = TMath::Abs(phi1 - phi2);

    if(dPhi > TMath::Pi())
      dPhi = TMath::TwoPi() - dPhi;

    return dPhi;

  }



  float ZinvxAODAnalysis :: DeltaR(float eta1, float eta2, float phi1, float phi2) {

    float dEta = eta1 - eta2;
    float dPhi = DeltaPhi(phi1,phi2);

    return TMath::Sqrt(dEta*dEta + dPhi*dPhi);

  }


