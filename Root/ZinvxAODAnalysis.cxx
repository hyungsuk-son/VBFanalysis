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
static SG::AuxElement::Decorator<char> dec_signal_forZ("signal_forZ");
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



  //TH1::SetDefaultSumw2(kTRUE);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  // get TEvent and TStore - must be done here b/c we need to retrieve CutBookkeepers container from TEvent!
  m_event = wk()->xaodEvent();

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

  if(!m_isData && m_isDerivation){

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
  m_isZnunu = true;
  m_isZmumu = true;
  m_isWmunu = true;
  m_isZee = true;
  m_isWenu = false;

  // Enable Systematics
  m_doSys = true;

  // Cut values
  m_muonPtCut = 7000.; /// MeV
  m_muonEtaCut = 2.5;
  m_elecPtCut = 7000.; /// MeV
  m_elecEtaCut = 2.47;
  m_photPtCut = 20000.; /// MeV
  m_photEtaCut = 2.47;
  m_jetPtCut = 20000.; /// MeV
  m_jetEtaCut = 4.5;
  m_monoJetPtCut = 120000.; /// MeV
  m_monoJetEtaCut = 2.4;
  m_sm1JetPtCut = 220000.; /// MeV
  m_sm1JetEtaCut = 2.4;
  m_diJet1PtCut = 80000.; /// MeV
  m_diJet2PtCut = 50000.; /// MeV
  m_diJetRapCut = 4.4;
  m_CJVptCut = 25000.; ///MeV
  m_metCut = 200000.; ///MeV
  m_mjjCut = 200000.; ///MeV
  m_LeadLepPtCut = 80000.; ///MeV
  m_SubLeadLepPtCut = 7000.; ///MeV
  m_ORJETdeltaR = 0.2;
  m_isoMuonPtMin = 10000.; ///MeV
  m_isoMuonPtMax = 500000.; ///MeV
  m_mllMin = 66000.; ///MeV
  m_mllMax = 116000.; ///MeV
  m_recoSF = true;
  m_idSF = true;
  m_ttvaSF = true; // for muon
  m_isoMuonSF = true;
  m_isoMuonSFforZ = false; // Muons in Z->mumu are not isolated
  m_trigSF = true; // for electron
  m_isoElectronSF = true;
  m_isoPtCut = false; // Initialize Setting (Do not change), This bool is for if you want to set pt range for muon Isolation
  // Blind cut
  m_METblindcut = 500000.; /// MeV
  m_Mjjblindcut = 750000.; /// MeV

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

  //////////
  // Muon //
  //////////
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


  //////////////
  // Electron //
  //////////////
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

  // Initialise Electron Efficiency Tool
  m_elecEfficiencySFTool_reco = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_reco");
  std::vector< std::string > corrFileNameList_reco;
  corrFileNameList_reco.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.RecoTrk.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("CorrectionFileNameList", corrFileNameList_reco) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->initialize() );

  m_elecEfficiencySFTool_id_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Loose");
  std::vector< std::string > corrFileNameList_id_Loose;
  corrFileNameList_id_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.LooseAndBLayerLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("CorrectionFileNameList", corrFileNameList_id_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->initialize() );

  m_elecEfficiencySFTool_id_Medium = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Medium");
  std::vector< std::string > corrFileNameList_id_Medium;
  corrFileNameList_id_Medium.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.MediumLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("CorrectionFileNameList", corrFileNameList_id_Medium) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->initialize() );

  m_elecEfficiencySFTool_id_Tight = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Tight");
  std::vector< std::string > corrFileNameList_id_Tight;
  corrFileNameList_id_Tight.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.TightLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->setProperty("CorrectionFileNameList", corrFileNameList_id_Tight) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->initialize() );

  m_elecEfficiencySFTool_iso_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Loose");
  std::vector< std::string > corrFileNameList_iso_Loose;
  corrFileNameList_iso_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("CorrectionFileNameList", corrFileNameList_iso_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->initialize() );

  m_elecEfficiencySFTool_iso_Medium = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Medium");
  std::vector< std::string > corrFileNameList_iso_Medium;
  corrFileNameList_iso_Medium.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.MediumLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("CorrectionFileNameList", corrFileNameList_iso_Medium) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->initialize() );

  m_elecEfficiencySFTool_iso_Tight = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Tight");
  std::vector< std::string > corrFileNameList_iso_Tight;
  corrFileNameList_iso_Tight.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.TightLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->setProperty("CorrectionFileNameList", corrFileNameList_iso_Tight) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->initialize() );

  m_elecEfficiencySFTool_trigEff = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_trigEff");
  std::vector< std::string > corrFileNameList_trigEff;
  corrFileNameList_trigEff.push_back("ElectronEfficiencyCorrection/efficiency.e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->setProperty("CorrectionFileNameList", corrFileNameList_trigEff) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->initialize() );

  m_elecEfficiencySFTool_trigSF_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_trigSF_Loose");
  std::vector< std::string > corrFileNameList_trigSF_Loose;
  corrFileNameList_trigSF_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->setProperty("CorrectionFileNameList", corrFileNameList_trigSF_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->setProperty("ForceDataType", 1) );
  //EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->setProperty("CorrelationModel", "FULL") );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF_Loose->initialize() );



  ////////////
  // Photon //
  ////////////
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


  ///////////////
  // Isolation //
  ///////////////
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
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("PhotonWP","FixedCutTight"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->initialize());

  /////////
  // Tau //
  /////////
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


  /////////
  // Jet //
  /////////
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
  //m_jvtag = 0;
  m_jvtag = new JetVertexTaggerTool("jvtag");
  //m_jvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
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

  /////////
  // MET //
  /////////
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
  auto t = m_toolBox.getTool("MuJetORT");
  EL_RETURN_CHECK("initialize()",t->setProperty("NumJetTrk", 5) );
  EL_RETURN_CHECK("initialize()",m_toolBox.initialize());


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
  file_ilumi.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/ZinvAnalysis/ilumicalc_histograms_None_276262-284484_OflLumi-13TeV-004.root"));
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("ConfigFiles", file_conf) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("LumiCalcFiles", file_ilumi) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactor",     1. / 1.16) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactorUP",   1.) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactorDOWN", 1. / 1.23) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("UnrepresentedDataAction", 2));
  if ( !(mcChannelNumber == 363121 || mcChannelNumber == 363351 || // One of the Ztautau or Wtaunu (from Valentinos)
        (mcChannelNumber >= 361063 && mcChannelNumber <= 361068) || mcChannelNumber == 361088 || mcChannelNumber == 361089) ) // Diboson samples
  { // These samples have missing mu values and the pileup reweighting tool doesn't like that and crashes.
    EL_RETURN_CHECK("initialize()",m_prwTool->initialize() );
  }    


  // Initialize PMGTools (PMGSherpa22VJetsWeightTool)
  // See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/CentralMC15ProductionList#Sherpa_v2_2_0_V_jets_NJet_reweig
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
/*
  std::vector<std::string> variations = {"JET_GroupedNP_1", "JET_GroupedNP_2", "JET_GroupedNP_3", "JET_EtaIntercalibration_NonClosure"};
  for (auto &s : variations) {
    CP::SystematicSet sysSet_up = {CP::SystematicVariation(s, 5)};
    CP::SystematicSet sysSet_down = {CP::SystematicVariation(s, -5)};
    m_sysList.push_back(sysSet_up);
    m_sysList.push_back(sysSet_down);
  }
*/

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
    if (m_isZmumu && !m_isZee && !m_isZnunu && m_doSys && ((sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos))) continue;
    if (m_isZee && !m_isZmumu && !m_isZnunu && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos))) continue;
    if (m_isZnunu && !m_isZmumu && !m_isZee && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos || sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos)) ) continue;

    //if (m_isZee && m_doSys && sysName.find("CorrUncertaintyNP")!=std::string::npos) continue; // Remove NP1~NP9, only choose Total error.

    //if (m_isZmumu && m_doSys && sysName != "" &&  sysName != "MUON_EFF_SYS__1down" && sysName != "MUON_EFF_SYS__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1down") continue;
    //if (m_isZee   && m_doSys && sysName != "" &&  sysName != "EL_EFF_ID_TotalCorrUncertainty__1down" && sysName != "EL_EFF_ID_TotalCorrUncertainty__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1down") continue
    //if (m_isZnunu && m_doSys && sysName != "" &&  sysName != "JET_EtaIntercalibration_Modelling__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1down") continue;



    // Number of Interactions
    if (sysName == ""){
      addHist(hMap1D, "avg_interaction"+sysName, 40, 0., 40.);
    }


    if (m_isZnunu) {
      h_channel = "h_znunu_";

      // For Ratio plot
      addHist(hMap1D, "Znunu_MET_mono"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "Znunu_MET_search"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "Znunu_Mjj_search"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, "Znunu_DeltaPhiAll"+sysName, nbinDPhi, binsDPhi);

      // For publication
      addHist(hMap1D, h_channel+"monojet_met"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_met"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_mjj"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_dPhijj"+sysName, nbinDPhi, binsDPhi);

      // Number of Interactions
      addHist(hMap1D, h_channel+"monojet_avg_interaction"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"vbf_avg_interaction"+sysName, 40, 0., 40.);

      if (sysName == ""){
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
        // Jets
        addHist(hMap1D, h_channel+"monojet_njet"+sysName, 40, 0., 40.);
        addHist(hMap1D, h_channel+"monojet_jet_pt"+sysName, 60, 0., 3000.);
        addHist(hMap1D, h_channel+"monojet_jet_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"monojet_jet_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"monojet_jet_rap"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"monojet_dPhimetjet"+sysName, 16, 0., 3.2);
        addHist(hMap1D, h_channel+"monojet_dPhiMinmetjet"+sysName, 16, 0., 3.2);
        ////////////////////
        // VBF phasespace //
        ////////////////////
        // Jets
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
        ////////////////////
        // SM1 phasespace //
        ////////////////////
        // For publication
        addHist(hMap1D, h_channel+"sm1_met"+sysName, nbinMET, binsMET);
        // Jets
        addHist(hMap1D, h_channel+"sm1_njet"+sysName, 40, 0., 40.);
        addHist(hMap1D, h_channel+"sm1_jet_pt"+sysName, 60, 0., 3000.);
        addHist(hMap1D, h_channel+"sm1_jet_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_jet_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_jet_rap"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_dPhimetjet"+sysName, 16, 0., 3.2);
        addHist(hMap1D, h_channel+"sm1_dPhiMinmetjet"+sysName, 16, 0., 3.2);

      }

    }



    if (m_isZmumu) {
      h_channel = "h_zmumu_";

      // For Ratio plot
      addHist(hMap1D, "Zmumu_MET_mono"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "Zmumu_MET_search"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "Zmumu_Mjj_search"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, "Zmumu_DeltaPhiAll"+sysName, nbinDPhi, binsDPhi);

      // For publication
      addHist(hMap1D, h_channel+"monojet_met_emulmet"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_met_emulmet"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_mjj"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_dPhijj"+sysName, nbinDPhi, binsDPhi);

      // Number of Interactions
      addHist(hMap1D, h_channel+"monojet_avg_interaction"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"vbf_avg_interaction"+sysName, 40, 0., 40.);

      if (sysName == ""){
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
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
        // Jets
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
        ////////////////////
        // SM1 phasespace //
        ////////////////////
        // For publication
        addHist(hMap1D, h_channel+"sm1_met_emulmet"+sysName, nbinMET, binsMET);
        // Jets
        addHist(hMap1D, h_channel+"sm1_njet"+sysName, 40, 0., 40.);
        addHist(hMap1D, h_channel+"sm1_jet_pt"+sysName, 60, 0., 3000.);
        addHist(hMap1D, h_channel+"sm1_jet_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_jet_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_jet_rap"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_dPhimetjet"+sysName, 16, 0., 3.2);
        addHist(hMap1D, h_channel+"sm1_dPhiMinmetjet"+sysName, 16, 0., 3.2);
        // Leptons
        addHist(hMap1D, h_channel+"sm1_lepton1_pt"+sysName, 30, 0., 1500.);
        addHist(hMap1D, h_channel+"sm1_lepton2_pt"+sysName, 30, 0., 1500.);
        addHist(hMap1D, h_channel+"sm1_lepton1_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_lepton2_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_lepton1_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_lepton2_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_mll"+sysName, 150, 0., 300.);

        //////////////////////////////////
        // MET Trigger Efficiency Study //
        //////////////////////////////////
        addHist(hMap1D, h_channel+"vbf_eff_study_met_emulmet"+sysName, 250, 0., 500.);
        addHist(hMap1D, h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70"+sysName, 250, 0., 500.);
        addHist(hMap1D, h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70_tclcw"+sysName, 250, 0., 500.);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met130"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met130"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70_tclcw"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70_tclcw"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met150"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met150"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70_tclcw"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70_tclcw"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met200"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met200"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70"+sysName, nbinDPhi, binsDPhi);
        addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70_tclcw"+sysName, nbinMjj, binsMjj);
        addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70_tclcw"+sysName, nbinDPhi, binsDPhi);
      }


      // Multijet Background study (Method 1)
      //
      ////////////////////////
      // Monojet phasespace //
      ////////////////////////
      addHist(hMap1D, h_channel+"monojet_multijet_study_mll_all_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"monojet_multijet_study_met_emulmet_all_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"monojet_multijet_study_mll_os_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"monojet_multijet_study_met_emulmet_os_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"monojet_multijet_study_mll_ss_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"monojet_multijet_study_met_emulmet_ss_lep"+sysName, nbinMET, binsMET);
      ////////////////////
      // VBF phasespace //
      ////////////////////
      addHist(hMap1D, h_channel+"vbf_multijet_study_mll_all_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"vbf_multijet_study_met_emulmet_all_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mjj_all_lep"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_multijet_study_dPhijj_all_lep"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mll_os_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"vbf_multijet_study_met_emulmet_os_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mjj_os_lep"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_multijet_study_dPhijj_os_lep"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mll_ss_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"vbf_multijet_study_met_emulmet_ss_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mjj_ss_lep"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_multijet_study_dPhijj_ss_lep"+sysName, nbinDPhi, binsDPhi);


      // Multijet Background study (Method 2)
      //
      // Nominal cut
      if (sysName == "") { // for Data and MC without systematics
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
        addHist(hMap1D, h_channel+"monojet_qcd_method2_nominal_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_nominal_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_nominal_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        ////////////////////
        // VBF phasespace //
        ////////////////////
        addHist(hMap1D, h_channel+"vbf_qcd_method2_nominal_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_nominal_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_nominal_cut_met500_inf_mll"+sysName, 150, 0., 300.);
      }
      //
      // Reverse cut
      if (m_isData) { // only for Data
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
        // Case 1
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case1_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case1_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case1_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 2
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case2_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case2_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case2_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 3
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case3_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case3_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case3_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        ////////////////////
        // VBF phasespace //
        ////////////////////
        // Case 1
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case1_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case1_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case1_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 2
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case2_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case2_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case2_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 3
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case3_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case3_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case3_cut_met500_inf_mll"+sysName, 150, 0., 300.);
      }

    } // m_isZmumu




    if (m_isZee) {
      h_channel = "h_zee_";

      // For Ratio plot
      addHist(hMap1D, "Zee_MET_mono"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "Zee_MET_search"+sysName, nbinMET, binsMET);
      addHist(hMap1D, "Zee_Mjj_search"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, "Zee_DeltaPhiAll"+sysName, nbinDPhi, binsDPhi);

      // For publication
      addHist(hMap1D, h_channel+"monojet_met_emulmet"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_met_emulmet"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_mjj"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_dPhijj"+sysName, nbinDPhi, binsDPhi);

      // Number of Interactions
      addHist(hMap1D, h_channel+"monojet_avg_interaction"+sysName, 40, 0., 40.);
      addHist(hMap1D, h_channel+"vbf_avg_interaction"+sysName, 40, 0., 40.);

      if (sysName == ""){
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
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
        // Jets
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
        ////////////////////
        // SM1 phasespace //
        ////////////////////
        // For publication
        addHist(hMap1D, h_channel+"sm1_met_emulmet"+sysName, nbinMET, binsMET);
        // Jets
        addHist(hMap1D, h_channel+"sm1_njet"+sysName, 40, 0., 40.);
        addHist(hMap1D, h_channel+"sm1_jet_pt"+sysName, 60, 0., 3000.);
        addHist(hMap1D, h_channel+"sm1_jet_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_jet_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_jet_rap"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_dPhimetjet"+sysName, 16, 0., 3.2);
        addHist(hMap1D, h_channel+"sm1_dPhiMinmetjet"+sysName, 16, 0., 3.2);
        // Leptons
        addHist(hMap1D, h_channel+"sm1_lepton1_pt"+sysName, 30, 0., 1500.);
        addHist(hMap1D, h_channel+"sm1_lepton2_pt"+sysName, 30, 0., 1500.);
        addHist(hMap1D, h_channel+"sm1_lepton1_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_lepton2_phi"+sysName, 32, -3.2, 3.2);
        addHist(hMap1D, h_channel+"sm1_lepton1_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_lepton2_eta"+sysName, 25, -5., 5.);
        addHist(hMap1D, h_channel+"sm1_mll"+sysName, 150, 0., 300.);
      }


      // Multijet Background study (Method 1)
      //
      ////////////////////////
      // Monojet phasespace //
      ////////////////////////
      addHist(hMap1D, h_channel+"monojet_multijet_study_mll_all_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"monojet_multijet_study_met_emulmet_all_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"monojet_multijet_study_mll_os_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"monojet_multijet_study_met_emulmet_os_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"monojet_multijet_study_mll_ss_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"monojet_multijet_study_met_emulmet_ss_lep"+sysName, nbinMET, binsMET);
      ////////////////////
      // VBF phasespace //
      ////////////////////
      addHist(hMap1D, h_channel+"vbf_multijet_study_mll_all_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"vbf_multijet_study_met_emulmet_all_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mjj_all_lep"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_multijet_study_dPhijj_all_lep"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mll_os_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"vbf_multijet_study_met_emulmet_os_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mjj_os_lep"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_multijet_study_dPhijj_os_lep"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mll_ss_lep"+sysName, 150, 0., 300.);
      addHist(hMap1D, h_channel+"vbf_multijet_study_met_emulmet_ss_lep"+sysName, nbinMET, binsMET);
      addHist(hMap1D, h_channel+"vbf_multijet_study_mjj_ss_lep"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_multijet_study_dPhijj_ss_lep"+sysName, nbinDPhi, binsDPhi);


      // Multijet Background study (Method 2)
      //
      // Nominal cut
      if (sysName == "") { // for Data and MC without systematics
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
        addHist(hMap1D, h_channel+"monojet_qcd_method2_nominal_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_nominal_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_nominal_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        ////////////////////
        // VBF phasespace //
        ////////////////////
        addHist(hMap1D, h_channel+"vbf_qcd_method2_nominal_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_nominal_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_nominal_cut_met500_inf_mll"+sysName, 150, 0., 300.);
      }
      //
      // Reverse cut
      if (m_isData) { // only for Data
        ////////////////////////
        // Monojet phasespace //
        ////////////////////////
        // Case 1
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case1_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case1_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case1_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 2
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case2_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case2_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case2_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 3
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case3_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case3_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"monojet_qcd_method2_case3_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        ////////////////////
        // VBF phasespace //
        ////////////////////
        // Case 1
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case1_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case1_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case1_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 2
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case2_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case2_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case2_cut_met500_inf_mll"+sysName, 150, 0., 300.);
        // Case 3
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case3_cut_met200_300_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case3_cut_met300_500_mll"+sysName, 150, 0., 300.);
        addHist(hMap1D, h_channel+"vbf_qcd_method2_case3_cut_met500_inf_mll"+sysName, 150, 0., 300.);
      }

    } // m_isZee




    //////////////////////////////////
    // MET Trigger Efficiency Study //
    //////////////////////////////////
    if (m_isWmunu && sysName == "") {
      h_channel = "h_wmunu_";
      // MET Trigger Efficiency Study
      addHist(hMap1D, h_channel+"vbf_eff_study_met_emulmet"+sysName, 250, 0., 500.);
      addHist(hMap1D, h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70"+sysName, 250, 0., 500.);
      addHist(hMap1D, h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70_tclcw"+sysName, 250, 0., 500.);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met130"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met130"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70_tclcw"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70_tclcw"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met150"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met150"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70_tclcw"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70_tclcw"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met200"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met200"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70"+sysName, nbinDPhi, binsDPhi);
      addHist(hMap1D, h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70_tclcw"+sysName, nbinMjj, binsMjj);
      addHist(hMap1D, h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70_tclcw"+sysName, nbinDPhi, binsDPhi);
    }





    // For cutflow comparison with Emily
    if (m_isZee || m_isZmumu) {
      if (m_isEmilyCutflow && sysName == ""){
        ////////////////////////////
        // Before Overlap Removal //
        ////////////////////////////
        addHist(hMap1D, "NTauBefore"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NEleBefore"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NMuBefore"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NMuZBefore"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NJetBefore"+sysName, 20, 0., 20.);
        ////////////////////////////
        // After Overlap Removal //
        ////////////////////////////
        addHist(hMap1D, "NTauAfter"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NEleAfter"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NMuAfter"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NMuZAfter"+sysName, 20, 0., 20.);
        addHist(hMap1D, "NJetAfter"+sysName, 20, 0., 20.);

        ///////////////////////////////////
        // Basic distribution comparison //
        ///////////////////////////////////
        addHist(hMap1D, "Emily_Zee_MET_mono"+sysName, 150,   0.,  1500.);
        addHist(hMap1D, "Emily_Zee_MET_search"+sysName, 30,    0.0,  1500.);
        addHist(hMap1D, "Emily_Zee_Mjj_search"+sysName, 80,    0.0,  4000.);
        addHist(hMap1D, "Emily_Zee_DeltaPhiAll"+sysName, 100, 0, TMath::Pi());
        addHist(hMap1D, "Emily_Zmumu_MET_mono"+sysName, 150,   0.,  1500.);
        addHist(hMap1D, "Emily_Zmumu_MET_search"+sysName, 30,    0.0,  1500.);
        addHist(hMap1D, "Emily_Zmumu_Mjj_search"+sysName, 80,    0.0,  4000.);
        addHist(hMap1D, "Emily_Zmumu_DeltaPhiAll"+sysName, 100, 0, TMath::Pi());
      }
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


  float print_mcWeight_origin = 1.;

  float mcWeight = 1.;
  if (!m_isData) {
    mcWeight = eventInfo->mcEventWeight();
    print_mcWeight_origin = mcWeight;

  }
  if (((mcChannelNumber > 363100 && mcChannelNumber < 363500) || (mcChannelNumber > 304010 && mcChannelNumber < 304025)) && std::abs(mcWeight) > 100.) mcWeight = 1.;

  float print_mcWeight_nohigh = mcWeight;
  double print_sherpaReweightValue = 1.;

  // Sherpa v2.2 V+jets NJet reweight 
  if (!m_isData) {
    // See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/CentralMC15ProductionList#NEW_Sherpa_v2_2_V_jets_NJet_rewe
    if (mcChannelNumber > 363100 &&  mcChannelNumber < 363500) { // only reweight W and Z samples
      double sherpaReweightValue = m_PMGSherpa22VJetsWeightTool->getWeight();
      print_sherpaReweightValue = sherpaReweightValue;
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
  const xAOD::Vertex* primVertex = 0;
  for (const auto &vtx : *vertices) {
    if (vtx->vertexType() == xAOD::VxType::PriVtx) {
      primVertex = vtx;
    }
  }

  //--------------
  // Preseletion
  //--------------
  if (vertices->size() < 1 || !primVertex) {
    Info("execute()", "WARNING: no primary vertex found! Skipping event.");
    return EL::StatusCode::SUCCESS;
  }
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




  ////////////////////////////
  // Create copy containers //
  ////////////////////////////

  xAOD::JetContainer* m_goodJet = new xAOD::JetContainer(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Jet>

  xAOD::MuonContainer* m_goodMuon = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);
  xAOD::MuonContainer* m_goodMuonForZ = new xAOD::MuonContainer(SG::VIEW_ELEMENTS); // only For Z->mumu selections (goodMuonForZ are the non-isolated muons)
  xAOD::MuonContainer* m_baselineMuon = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);

  xAOD::ElectronContainer* m_goodElectron = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
  xAOD::ElectronContainer* m_baselineElectron = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);

  xAOD::TauJetContainer* m_goodTau = new xAOD::TauJetContainer(SG::VIEW_ELEMENTS);

  xAOD::PhotonContainer* m_goodPhoton = new xAOD::PhotonContainer(SG::VIEW_ELEMENTS);

  // Create a MissingETContainer with its aux store for each systematic
  xAOD::MissingETContainer* m_met = new xAOD::MissingETContainer();
  xAOD::MissingETAuxContainer* m_metAux = new xAOD::MissingETAuxContainer();
  m_met->setStore(m_metAux);






  //-----------------------
  // Systematics Start
  //-----------------------
  // loop over recommended systematics
  for (const auto &sysList : m_sysList){
    std::string sysName = (sysList).name();
    if ((!m_doSys || m_isData) && sysName != "") continue;

    if (m_doSys && (sysName.find("TAUS_")!=std::string::npos || sysName.find("PH_")!=std::string::npos )) continue;
    if (m_isZmumu && !m_isZee && !m_isZnunu && m_doSys && ((sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos))) continue;
    if (m_isZee && !m_isZmumu && !m_isZnunu && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos))) continue;
    if (m_isZnunu && !m_isZmumu && !m_isZee && m_doSys && ((sysName.find("MUON_")!=std::string::npos || sysName.find("MUONS_")!=std::string::npos || sysName.find("EL_")!=std::string::npos || sysName.find("EG_")!=std::string::npos)) ) continue;

    //if (m_isZee && m_doSys && sysName.find("CorrUncertaintyNP")!=std::string::npos) continue; // Remove NP1~NP9, only choose Total error.

    //if (m_isZmumu && m_doSys && sysName != "" &&  sysName != "MUON_EFF_SYS__1down" && sysName != "MUON_EFF_SYS__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1down") continue;
    //if (m_isZee   && m_doSys && sysName != "" &&  sysName != "EL_EFF_ID_TotalCorrUncertainty__1down" && sysName != "EL_EFF_ID_TotalCorrUncertainty__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1down") continue
    //if (m_isZnunu && m_doSys && sysName != "" &&  sysName != "JET_EtaIntercalibration_Modelling__1up" &&  sysName != "JET_EtaIntercalibration_Modelling__1down") continue;

    // Print the list of systematics
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




    float print_puweight = 1.;

    //---------------------
    // Pile-up reweighting
    //---------------------
    if (!m_isData) {
      if ( mcChannelNumber == 363121 || mcChannelNumber == 363351 || // One of the Ztautau or Wtaunu (from Valentinos)
          (mcChannelNumber >= 361063 && mcChannelNumber <= 361068) || mcChannelNumber == 361088 || mcChannelNumber == 361089 ) // Diboson samples
      { // These samples have missing mu values and the pileup reweighting tool doesn't like that and crashes.
        mcEventWeight = mcWeight;
      }
      else {
        float pu_weight = m_prwTool->getCombinedWeight(*eventInfo); // Get Pile-up weight
        print_puweight = pu_weight;
        mcEventWeight = mcWeight * pu_weight;
      }
    }



    ///////////////////////////
    // Clear copy containers //
    ///////////////////////////

    m_goodJet->clear();
    m_goodMuon->clear();
    m_goodMuonForZ->clear();
    m_baselineMuon->clear();
    m_goodElectron->clear();
    m_baselineElectron->clear();
    m_goodTau->clear();
    m_goodPhoton->clear();
    m_met->clear();




    // Average Interaction
    float m_AverageInteractionsPerCrossing = 0.;
    if (!m_isData) {
      m_AverageInteractionsPerCrossing = eventInfo->averageInteractionsPerCrossing();
    }
    else {
      m_AverageInteractionsPerCrossing = m_prwTool->getCorrectedMu(*eventInfo, true);
    }

    if (sysName == ""){
      hMap1D["avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight);
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
      // TauOverlappingElectronLLHDecorator
      m_tauOverlappingElectronLLHDecorator->decorate(*taujet);
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

    // iterate over our shallow copy
    for (const auto& jets : *jetSC) { // C++11 shortcut
      //Info("execute()", "  original jet pt = %.2f GeV", jets->pt() * 0.001);

      // According to https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15

      // JES calibration
      if ( !m_jetCalibration->applyCalibration(*jets).isSuccess() ){
        Error("execute()", "Failed to apply calibration to Jet objects. Exiting." );
        return EL::StatusCode::FAILURE;
      }

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
      float newjvt = m_jvtag->updateJvt(*jets);
      acc_jvt(*jets) = newjvt;

      //Info("execute()", "  corrected jet pt = %.2f GeV", jets->pt() * 0.001);
      //Info("execute()", "  updated jet jvt = %.2f ", newjvt);

      dec_signal(*jets) = false;
      selectDec(*jets) = false; // To select objects for Overlap removal

    } // end for loop over shallow copied jets

    // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
    // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
    // The method is defined in the header file xAODBase/IParticleHelpers.h
    bool setLinksJet = xAOD::setOriginalObjectLink(*m_jets,*jetSC);
    if(!setLinksJet) {
      Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
      return StatusCode::FAILURE;
    }



    // -----------------
    // Select Good Jet
    // -----------------
    /// Creating New Hard Object Containers
    // [For jet identification] filter the Jet container m_jets, placing selected jets into m_goodJet

    bool isBadJet = false;

    // iterate over our shallow copy
    for (const auto& jets : *jetSC) { // C++11 shortcut

      // Jet Signal Selection
      if (IsSignalJet(*jets)) {

        m_goodJet->push_back( jets );
      }
    } // end for loop over shallow copied jets



    //-----------------------------------------------
    // Define Good Leptons and Calculate Scale Factor
    //-----------------------------------------------

    ///////////////
    // Good Muon //
    ///////////////
    // iterate over our shallow copy
    for (const auto& muon : *muonSC) { // C++11 shortcut
      // Muon Selection for VBF study
      if (dec_signal(*muon)) {
        m_goodMuon->push_back( muon );
        //Info("execute()", "  Good muon pt = %.2f GeV", (muon->pt() * 0.001));
      }
      if (dec_signal_forZ(*muon)) {
        muon->auxdata<bool>("brem") = false; // For overlap removal with electron
        m_goodMuonForZ->push_back( muon );
      }
      if (dec_baseline(*muon)) {
        m_baselineMuon->push_back( muon ); // For QCD multijet study
      }
    } // end for loop over shallow copied muons

    ///////////////////
    // Good Electron //
    ///////////////////
    // iterate over our shallow copy
    for (const auto& electron : *elecSC) { // C++11 shortcut
      // Electron Selection for VBF study
      if (dec_signal(*electron)) {
          m_goodElectron->push_back( electron );
      }
      if (dec_baseline(*electron)) {
        m_baselineElectron->push_back( electron ); // For QCD multijet study
      }
    } // end for loop over shallow copied electrons

    //////////////
    // Good Tau //
    //////////////
    // iterate over our shallow copy
    for (const auto& taujet : *tauSC) { // C++11 shortcut
      // Tau Selection for VBF study
      if (dec_signal(*taujet)) {
        m_goodTau->push_back( taujet );
        //Info("execute()", "  Good tau pt = %.2f GeV", (taujet->pt() * 0.001));
      }
    } // end for loop over shallow copied taus

    /////////////////
    // Good Photon //
    /////////////////
    // iterate over our shallow copy
    for (const auto& photon : *photSC) { // C++11 shortcut
      // Photon Selection for VBF study
      if (dec_signal(*photon)) {
        m_goodPhoton->push_back( photon );
      }
    } // end for loop over shallow copied photons





    /////////////////////////////////
    // Sort Good Muon and Electron //
    /////////////////////////////////
    /*
    // Muon
    if (m_goodMuon->size() > 1) std::sort(m_goodMuon->begin(), m_goodMuon->end(), DescendingPt());
    if (m_goodMuonForZ->size() > 1) std::sort(m_goodMuonForZ->begin(), m_goodMuonForZ->end(), DescendingPt());
    if (m_baselineMuon->size() > 1) std::sort(m_baselineMuon->begin(), m_baselineMuon->end(), DescendingPt());
    // Electron
    if (m_goodElectron->size() > 1) std::sort(m_goodElectron->begin(), m_goodElectron->end(), DescendingPt());
    if (m_baselineElectron->size() > 1) std::sort(m_baselineElectron->begin(), m_baselineElectron->end(), DescendingPt());
    */

    if (m_goodMuonForZ->size() > 1) std::partial_sort(m_goodMuonForZ->begin(), m_goodMuonForZ->begin()+2, m_goodMuonForZ->end(), DescendingPt());
    if (m_goodMuon->size() > 1) std::partial_sort(m_goodMuon->begin(), m_goodMuon->begin()+2, m_goodMuon->end(), DescendingPt());
    if (m_baselineMuon->size() > 1) std::partial_sort(m_baselineMuon->begin(), m_baselineMuon->begin()+2, m_baselineMuon->end(), DescendingPt());
    if (m_goodElectron->size() > 1) std::partial_sort(m_goodElectron->begin(), m_goodElectron->begin()+2, m_goodElectron->end(), DescendingPt());
    if (m_baselineElectron->size() > 1) std::partial_sort(m_baselineElectron->begin(), m_baselineElectron->begin()+2, m_baselineElectron->end(), DescendingPt());




    //----------------------------------------------------
    // Decorate overlapped objects using official OR Tool
    //----------------------------------------------------

    if ( !m_orTool->removeOverlaps(m_goodElectron, m_goodMuon, m_goodJet, m_goodTau, m_goodPhoton).isSuccess() ){
      Error("execute()", "Failed to apply the overlap removal to all objects. Exiting." );
      return EL::StatusCode::FAILURE;
    }
/*
    // Now, dump all of the results
    if (sysName == "") {

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
*/



    // Cutflow comparison with Emily
    ////////////////////////////
    // Before Overlap Removal //
    ////////////////////////////
    if (m_isEmilyCutflow && sysName == "") {
      if ( (m_isZee || m_isZmumu) ){
        hMap1D["NTauBefore"+sysName]->Fill(m_goodTau->size(),1.0);
        hMap1D["NEleBefore"+sysName]->Fill(m_goodElectron->size(),1.0);
        hMap1D["NMuBefore"+sysName]->Fill(m_goodMuon->size(),1.0);
        hMap1D["NMuZBefore"+sysName]->Fill(m_goodMuonForZ->size(),1.0);
        hMap1D["NJetBefore"+sysName]->Fill(m_goodJet->size(),1.0);
        /*
        if (m_goodMuonForZ->size() != m_goodMuon->size()){
          Info("execute()", "============================");
          Info("execute()", " Event # = %llu", eventInfo->eventNumber());
          Info("execute()", "======= goodMuonForZ =======");
          //int muonZCount = 0;
          for (const auto& muon : *m_goodMuonForZ) {
            //muonZCount++;
            //Info("execute()", " muonZ # : %i", muonCount);
            Info("execute()", " muonZ pt = %.4f GeV", muon->pt() * 0.001);
            //Info("execute()", " muonZ eta = %.3f", muon->eta());
            //Info("execute()", " muonZi phi = %.3f", muon->phi());
          }
          Info("execute()", "========= goodMuon =========");
          //int muonCount = 0;
          for (const auto& muon : *m_goodMuon) {
            //muonCount++;
            //Info("execute()", " muon # : %i", muonCount);
            Info("execute()", " muon pt = %.4f GeV", muon->pt() * 0.001);
            //Info("execute()", " muon eta = %.3f", muon->eta());
            //Info("execute()", " muon phi = %.3f", muon->phi());
          }
        }
        if (eventInfo->eventNumber() == 52 || eventInfo->eventNumber() == 798 || eventInfo->eventNumber() == 1300){
          Info("execute()", "============================");
          Info("execute()", " Event # = %llu", eventInfo->eventNumber());
          Info("execute()", "======= goodMuonForZ =======");
          //int muonZCount = 0;
          for (const auto& muon : *m_goodMuonForZ) {
            Info("execute()", " muonZ pt = %f GeV", muon->pt() * 0.001);
          }
          Info("execute()", "========= goodMuon =========");
          for (const auto& muon : *m_goodMuon) {
            Info("execute()", " muon pt = %f GeV", muon->pt() * 0.001);
          }
        }

        / Tau Study
        if (m_goodTau->size() > 0){
          Info("execute()", "=====================================");
          Info("execute()", " Event # = %llu", eventInfo->eventNumber());
          int tauCount = 0;
          for (const auto& tau : *m_goodTau) {
            tauCount++;
            Info("execute()", " tau # : %i", tauCount);
            Info("execute()", " tau pt = %.3f GeV", tau->pt() * 0.001);
            //Info("execute()", " tau charge = %.1f", tau->charge());
            Info("execute()", " tau eta = %.3f", tau->eta());
            Info("execute()", " tau phi = %.3f", tau->phi());
          }
        }
        */
      }
    }





    ////////////////////////////////////////////
    // Overlap removal manually for VBF study //
    ////////////////////////////////////////////
    // For Z->mumu good muons
    // This is done before the erasing of the overlapping objects is done
    int i=0;
    for (const auto& muon : *m_goodMuonForZ) {
      for (const auto& jet : *m_goodJet) { // C++11 shortcut
        if (deltaR(jet->eta(), muon->eta(), jet->phi(), muon->phi()) < m_ORJETdeltaR) {

          int ntrks = 0;
          float sumpt = 0.;
          std::vector<int> ntrks_vec = jet->auxdata<std::vector<int> >("NumTrkPt500");
          std::vector<float> sumpt_vec = jet->auxdata<std::vector<float> >("SumPtTrkPt500");
          if (ntrks_vec.size() > 0) {
            ntrks = ntrks_vec[primVertex->index()];
            sumpt = sumpt_vec[primVertex->index()];
          }
          if (ntrks < 5 || (muon->pt()/jet->pt() > 0.5 && muon->pt()/sumpt > 0.7)) {
            muon->auxdata<bool>("brem") = true;
          }   else if (muon->pt() < 20000.) {
            m_goodMuonForZ->erase(m_goodMuonForZ->begin()+i);
            i--;
            break;
          } //else  muon->auxdata<bool>("overlap") = true;
        }
      } // break to here
      i++;
    }







    //////////////////////////////////////////////
    // Overlap removal officially for VBF study //
    //////////////////////////////////////////////
    m_goodMuon->erase(std::remove_if(std::begin(*m_goodMuon), std::end(*m_goodMuon), [](xAOD::Muon* mu) {return overlapAcc(*mu);}), std::end(*m_goodMuon));
    m_goodElectron->erase(std::remove_if(std::begin(*m_goodElectron), std::end(*m_goodElectron), [](xAOD::Electron* elec) {return overlapAcc(*elec);}), std::end(*m_goodElectron));
    m_goodJet->erase(std::remove_if(std::begin(*m_goodJet), std::end(*m_goodJet), [](xAOD::Jet* jet) {return overlapAcc(*jet);}), std::end(*m_goodJet));
    m_goodTau->erase(std::remove_if(std::begin(*m_goodTau), std::end(*m_goodTau), [](xAOD::TauJet* tau) {return overlapAcc(*tau);}), std::end(*m_goodTau));
    m_goodPhoton->erase(std::remove_if(std::begin(*m_goodPhoton), std::end(*m_goodPhoton), [](xAOD::Photon* phot) {return overlapAcc(*phot);}), std::end(*m_goodPhoton));




    // loop round electrons and remove if it is close to a muon that has bremmed and likely faked an electron
    int j=0;
    for (const auto &electron : *m_goodElectron) {
      for (const auto &muon : *m_goodMuonForZ) {
        if (deltaR(electron->caloCluster()->etaBE(2),muon->eta(),electron->phi(),muon->phi()) < 0.3 && muon->auxdata<bool>("brem") == true) {
          m_goodElectron->erase(m_goodElectron->begin()+j);
          j--;
          break;
        }
      } // break to here
      j++;
    }


    // loop round taus and remove if it is close to a muon that has bremmed and likely faked an electron
    int jj=0;
    for (const auto &tau : *m_goodTau) {
      for (const auto &muon : *m_goodMuonForZ) {
        if (deltaR(tau->eta(),muon->eta(),tau->phi(),muon->phi()) < 0.3) {
          m_goodTau->erase(m_goodTau->begin()+jj);
          jj--;
          break;
        }
      } // break to here
      jj++;
    }

    // Cutflow comparison with Emily
    ////////////////////////////
    // After Overlap Removal //
    ////////////////////////////
    if (m_isEmilyCutflow && sysName == "") {
      if ( (m_isZee || m_isZmumu) ){
        hMap1D["NTauAfter"+sysName]->Fill(m_goodTau->size(),1.0);
        hMap1D["NEleAfter"+sysName]->Fill(m_goodElectron->size(),1.0);
        hMap1D["NMuAfter"+sysName]->Fill(m_goodMuon->size(),1.0);
        hMap1D["NMuZAfter"+sysName]->Fill(m_goodMuonForZ->size(),1.0);
        hMap1D["NJetAfter"+sysName]->Fill(m_goodJet->size(),1.0);
        /*
        if (m_goodTau->size() > 0){
          Info("execute()", "=====================================");
          Info("execute()", " Event # = %llu", eventInfo->eventNumber());
          int tauCount = 0;
          for (const auto& tau : *m_goodTau) {
            tauCount++;
            Info("execute()", " tau # : %i", tauCount);
            Info("execute()", " tau pt = %.3f GeV", tau->pt() * 0.001);
            //Info("execute()", " tau charge = %.1f", tau->charge());
            Info("execute()", " tau eta = %.3f", tau->eta());
            Info("execute()", " tau phi = %.3f", tau->phi());
          }
        }
        */
      }
    }






    //------------------
    // Bad Jet Decision 
    //------------------
    // iterate over our shallow copy
    for (const auto& jets : *m_goodJet) { // C++11 shortcut

      // Veto Jet (cleaning Jet)
      if (IsBadJet(*jets)) isBadJet = true;

    } // end for loop over shallow copied jets


    //------------------------------------
    // Event Cleaning (Jet cleaning tool)
    //------------------------------------
    //if (isBadJet) return EL::StatusCode::SUCCESS;
    if (isBadJet){

      //////////////////////////////////
      // Delete shallow copy containers
      //////////////////////////////////

      // The containers created by the shallow copy are owned by you. Remember to delete them
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


      continue; // escape from the systematic loop
    }
    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Jet Cleaning");
    if (m_useArrayCutflow) m_eventCutflow[4]+=1;





    //////////////////////
    // Sort Good Jets //
    //////////////////////
    m_goodJet->erase(std::remove_if(std::begin(*m_goodJet), std::end(*m_goodJet), [](xAOD::Jet* jet) {return (jet->rapidity() > 4.4);}), std::end(*m_goodJet));
    std::sort(m_goodJet->begin(), m_goodJet->end(), DescendingPt());





    //==============//
    // MET building //
    //==============//


    // do CST or TST
    //std::string softTerm = "SoftClus";
    std::string softTerm = "PVSoftTrk";


    // Real MET
    float MET = -9e9;
    float MET_phi = -9e9;
    // Wenu MET
    float emulMET_Wenu = -9e9;
    float emulMET_Wenu_phi = -9e9;
    // Zee MET
    float emulMET_Zee = -9e9;
    float emulMET_Zee_phi = -9e9;
    // Wmunu MET
    float emulMET_Wmunu = -9e9;
    float emulMET_Wmunu_phi = -9e9;
    // Zmumu MET
    float emulMET_Zmumu = -9e9;
    float emulMET_Zmumu_phi = -9e9;


    //=============================
    // Create MissingETContainers 
    //=============================


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



    // It is necessary to reset the selected objects before every MET calculation
    m_met->clear();
    m_metMap->resetObjSelectionFlags();



    //===========================
    // For rebuild the real MET
    //===========================

    // Electron
    //-----------------
    /// Creat New Hard Object Containers
    // [For MET building] filter the Electron container m_electrons, placing selected electrons into m_MetElectrons
    ConstDataVector<xAOD::ElectronContainer> m_MetElectrons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Electron>

    // iterate over our shallow copy
    for (const auto& electron : *m_goodElectron) { // C++11 shortcut
      // For MET rebuilding
      m_MetElectrons.push_back( electron );
    } // end for loop over shallow copied electrons
    //const xAOD::ElectronContainer* p_MetElectrons = m_MetElectrons.asDataVector();

    // For real MET
    m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
        xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
        m_met,                                      //filling this met container
        m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
        m_metMap);                                  //and this association map


    /*
    // Photon
    //-----------------
    /// Creat New Hard Object Containers
    // [For MET building] filter the Photon container m_photons, placing selected photons into m_MetPhotons
    ConstDataVector<xAOD::PhotonContainer> m_MetPhotons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Photon>

    // iterate over our shallow copy
    for (const auto& photon : *m_goodPhoton) { // C++11 shortcut
      // For MET rebuilding
      m_MetPhotons.push_back( photon );
    } // end for loop over shallow copied photons

    // For real MET
    m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
        xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
        m_met,                                    //filling this met container
        m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
        m_metMap);                                //and this association map

    */

    // TAUS
    //-----------------
    /// Creat New Hard Object Containers
    // [For MET building] filter the TauJet container m_taus, placing selected taus into m_MetTaus
    ConstDataVector<xAOD::TauJetContainer> m_MetTaus(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::TauJet>

    // iterate over our shallow copy
    for (const auto& taujet : *m_goodTau) { // C++11 shortcut
      // For MET rebuilding
      m_MetTaus.push_back( taujet );
    } // end for loop over shallow copied taus

    // For real MET
    m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
        xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
        m_met,                                 //filling this met container
        m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
        m_metMap);                             //and this association map


    // Muon
    //-----------------
    /// Creat New Hard Object Containers
    // [For MET building] filter the Muon container m_muons, placing selected muons into m_MetMuons
    ConstDataVector<xAOD::MuonContainer> m_MetMuons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Muon>

    // iterate over our shallow copy
    for (const auto& muon : *m_goodMuon) { // C++11 shortcut
      // For MET rebuilding
      m_MetMuons.push_back( muon );
    } // end for loop over shallow copied muons
    // For real MET
    m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
        xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
        m_met,                                  //filling this met container
        m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
        m_metMap);                              //and this association map


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




    /////////////////////////////
    // Soft term uncertainties //
    /////////////////////////////
    if (!m_isData) {
      // Get the track soft term (For real MET)
      xAOD::MissingET* softTrkmet = (*m_met)[softTerm];
      if (m_metSystTool->applyCorrection(*softTrkmet) != CP::CorrectionCode::Ok) {
        Error("execute()", "METSystematicsTool returns Error CorrectionCode");
      }
    }


    ///////////////
    // MET Build //
    ///////////////
    //m_metMaker->rebuildTrackMET("RefJetTrk", softTerm, m_met, jetSC, m_metCore, m_metMap, true);

    //this builds the final track or cluster met sums, using systematic varied container
    //In the future, you will be able to run both of these on the same container to easily output CST and TST

    // For real MET
    m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());




    ///////////////////
    // Fill real MET //
    ///////////////////

    //MET_ex = ((*m_met)["Final"]->mpx());
    //MET_ey = ((*m_met)["Final"]->mpy());
    MET = ((*m_met)["Final"]->met());
    //SumET = ((*m_met)["Final"]->sumet());
    MET_phi = ((*m_met)["Final"]->phi());








    //======================================================================
    // For rebuild the emulated MET for Wenu (by marking Electron invisible)
    //======================================================================

    if (m_isWenu) {

      // It is necessary to reset the selected objects before every MET calculation
      m_met->clear();
      m_metMap->resetObjSelectionFlags();


      // Electron
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Electron container m_electrons, placing selected electrons into m_MetElectrons
      //
      // For emulated MET (No electrons)
      // Make a empty container for invisible electrons
      ConstDataVector<xAOD::ElectronContainer> m_EmptyElectrons(SG::VIEW_ELEMENTS);
      m_EmptyElectrons.clear();
      m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
          xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
          m_met,                           //filling this met container
          m_EmptyElectrons.asDataVector(),            //using these metElectrons that accepted our cuts
          m_metMap);                       //and this association map


      /*
      // Photon
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Photon container m_photons, placing selected photons into m_MetPhotons
      // For emulated MET marking electrons invisible
      m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
          xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
          m_met,                         //filling this met container
          m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
          m_metMap);                     //and this association map
      */


      // TAUS
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the TauJet container m_taus, placing selected taus into m_MetTaus
      // For emulated MET marking electrons invisible
      m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
          xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
          m_met,                      //filling this met container
          m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
          m_metMap);                  //and this association map


      // Muon
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Muon container m_muons, placing selected muons into m_MetMuons
      //
      // For emulated MET (No electrons)
      m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
          xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
          m_met,                       //filling this met container
          m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
          m_metMap);                   //and this association map



      // JET
      //-----------------
      //Now time to rebuild jetMet and get the soft term
      //This adds the necessary soft term for both CST and TST
      //these functions create an xAODMissingET object with the given names inside the container

      // For emulated MET marking electrons invisible
      m_metMaker->rebuildJetMET("RefJet",          //name of jet met
          "SoftClus",           //name of soft cluster term met
          "PVSoftTrk",          //name of soft track term met
          m_met,     //adding to this new met container
          jetSC,                //using this jet collection to calculate jet met
          m_metCore, //core met container
          m_metMap,  //with this association map
          true);                //apply jet jvt cut




      /////////////////////////////
      // Soft term uncertainties //
      /////////////////////////////
      if (!m_isData) {
        // Get the track soft term for Wenu (For emulated MET marking electrons invisible)
        xAOD::MissingET* softTrkmet = (*m_met)[softTerm];
        if (m_metSystTool->applyCorrection(*softTrkmet) != CP::CorrectionCode::Ok) {
          Error("execute()", "METSystematicsTool returns Error CorrectionCode");
        }
      }



      ///////////////
      // MET Build //
      ///////////////
      // For emulated MET for Wenu marking electrons invisible
      m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());



      /////////////////////////////////////////////////////////////////
      // Fill emulated MET for Wenu (by marking electrons invisible) //
      /////////////////////////////////////////////////////////////////
      //emulMET_Wenu_ex = ((*m_met)["Final"]->mpx());
      //emulMET_Wenu_ey = ((*m_met)["Final"]->mpy());
      emulMET_Wenu = ((*m_met)["Final"]->met());
      //emulSumET_Wenu = ((*m_met)["Final"]->sumet());
      emulMET_Wenu_phi = ((*m_met)["Final"]->phi());



    } // m_isWenu








    //=====================================================================
    // For rebuild the emulated MET for Zee (by marking Electron invisible)
    //=====================================================================

    if (m_isZee) {

      // It is necessary to reset the selected objects before every MET calculation
      m_met->clear();
      m_metMap->resetObjSelectionFlags();


      // Electron
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Electron container m_electrons, placing selected electrons into m_MetElectrons
      //
      // For emulated MET (No electrons)
      // Make a empty container for invisible electrons
      /*
         ConstDataVector<xAOD::ElectronContainer> m_EmptyElectrons(SG::VIEW_ELEMENTS);
         m_EmptyElectrons.clear();
         m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
         xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
         m_met,                           //filling this met container
         m_EmptyElectrons.asDataVector(),            //using these metElectrons that accepted our cuts
         m_metMap);                       //and this association map
         */
      // Make a container for invisible electrons
      ConstDataVector<xAOD::ElectronContainer> m_invisibleElectrons(SG::VIEW_ELEMENTS);
      for (const auto& electron : *m_goodElectron) { // C++11 shortcut
        m_invisibleElectrons.push_back( electron );
      }
      // Mark electrons invisible (No electrons)
      m_metMaker->markInvisible(m_invisibleElectrons.asDataVector(), m_metMap);

      // Not adding Photon, Tau, Muon objects as we veto on additional leptons and photons might be an issue for muon FSR


      // JET
      //-----------------
      //Now time to rebuild jetMet and get the soft term
      //This adds the necessary soft term for both CST and TST
      //these functions create an xAODMissingET object with the given names inside the container

      // For emulated MET marking electrons invisible
      m_metMaker->rebuildJetMET("RefJet",          //name of jet met
          "SoftClus",           //name of soft cluster term met
          "PVSoftTrk",          //name of soft track term met
          m_met,     //adding to this new met container
          jetSC,                //using this jet collection to calculate jet met
          m_metCore, //core met container
          m_metMap,  //with this association map
          true);                //apply jet jvt cut



      /////////////////////////////
      // Soft term uncertainties //
      /////////////////////////////
      if (!m_isData) {
        // Get the track soft term for Zee (For emulated MET marking electrons invisible)
        xAOD::MissingET* softTrkmet = (*m_met)[softTerm];
        if (m_metSystTool->applyCorrection(*softTrkmet) != CP::CorrectionCode::Ok) {
          Error("execute()", "METSystematicsTool returns Error CorrectionCode");
        }
      }


      ///////////////
      // MET Build //
      ///////////////
      // For emulated MET for Zee marking electrons invisible
      m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());


      ////////////////////////////////////////////////////////////////
      // Fill emulated MET for Zee (by marking electrons invisible) //
      ////////////////////////////////////////////////////////////////
      //emulMET_Zee_ex = ((*m_met)["Final"]->mpx());
      //emulMET_Zee_ey = ((*m_met)["Final"]->mpy());
      emulMET_Zee = ((*m_met)["Final"]->met());
      //emulSumET_Zee = ((*m_met)["Final"]->sumet());
      emulMET_Zee_phi = ((*m_met)["Final"]->phi());


    } // m_isZee





    //===================================================================
    // For rebuild the emulated MET for Wmunu (by marking Muon invisible)
    //===================================================================

    if (m_isWmunu) {

      // It is necessary to reset the selected objects before every MET calculation
      m_met->clear();
      m_metMap->resetObjSelectionFlags();

      // Electron
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Electron container m_electrons, placing selected electrons into m_MetElectrons
      // For emulated MET (No muons)
      m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
          xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
          m_met,                             //filling this met container
          m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
          m_metMap);                         //and this association map


      /*
      // Photon
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Photon container m_photons, placing selected photons into m_MetPhotons
      // For emulated MET marking muons invisible
      m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
          xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
          m_met,                           //filling this met container
          m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
          m_metMap);                       //and this association map
      */


      // TAUS
      //-----------------
      //
      /// Creat New Hard Object Containers
      // [For MET building] filter the TauJet container m_taus, placing selected taus into m_MetTaus
      // For emulated MET marking muons invisible
      m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
          xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
          m_met,                        //filling this met container
          m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
          m_metMap);                    //and this association map



      /*
      // Muon
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Muon container m_muons, placing selected muons into m_MetMuons
      //
      // For emulated MET (No muons)
      // Make a empty container for invisible electrons
      ConstDataVector<xAOD::MuonContainer> m_EmptyMuons(SG::VIEW_ELEMENTS);
      m_EmptyMuons.clear();
      m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
          xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
          m_met,                         //filling this met container
          m_EmptyMuons.asDataVector(),            //using these metMuons that accepted our cuts
          m_metMap);                     //and this association map
      */
      // Make a container for invisible muons
      ConstDataVector<xAOD::MuonContainer> m_invisibleMuons(SG::VIEW_ELEMENTS);
      for (const auto& muon : *m_goodMuon) { // C++11 shortcut
        m_invisibleMuons.push_back( muon );
      }
      // Mark muons invisible
      m_metMaker->markInvisible(m_invisibleMuons.asDataVector(), m_metMap);


      // JET
      //-----------------
      //Now time to rebuild jetMet and get the soft term
      //This adds the necessary soft term for both CST and TST
      //these functions create an xAODMissingET object with the given names inside the container

      // For emulated MET marking muons invisible
      m_metMaker->rebuildJetMET("RefJet",          //name of jet met
          "SoftClus",           //name of soft cluster term met
          "PVSoftTrk",          //name of soft track term met
          m_met,       //adding to this new met container
          jetSC,                //using this jet collection to calculate jet met
          m_metCore,   //core met container
          m_metMap,    //with this association map
          true);                //apply jet jvt cut




      /////////////////////////////
      // Soft term uncertainties //
      /////////////////////////////
      if (!m_isData) {
        // Get the track soft term for Wmunu (For emulated MET marking muons invisible)
        xAOD::MissingET* softTrkmet = (*m_met)[softTerm];
        if (m_metSystTool->applyCorrection(*softTrkmet) != CP::CorrectionCode::Ok) {
          Error("execute()", "METSystematicsTool returns Error CorrectionCode");
        }
      }



      ///////////////
      // MET Build //
      ///////////////
      // For emulated MET for Wmunu marking muons invisible
      m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());


      //////////////////////////////////////////////////////////////
      // Fill emulated MET for Wmunu (by marking muons invisible) //
      //////////////////////////////////////////////////////////////
      //emulMET_Wmunu_ex = ((*m_met)["Final"]->mpx());
      //emulMET_Wmunu_ey = ((*m_met)["Final"]->mpy());
      emulMET_Wmunu = ((*m_met)["Final"]->met());
      //emulSumET_Wmunu = ((*m_met)["Final"]->sumet());
      emulMET_Wmunu_phi = ((*m_met)["Final"]->phi());


    } // m_isWmunu








    //===================================================================
    // For rebuild the emulated MET for Zmumu (by marking Muon invisible)
    //===================================================================

    if (m_isZmumu) {

      // It is necessary to reset the selected objects before every MET calculation
      m_met->clear();
      m_metMap->resetObjSelectionFlags();


      // Not adding Electron, Photon, Tau objects as we veto on additional leptons and photons might be an issue for muon FSR

      // Muon
      //-----------------
      /// Creat New Hard Object Containers
      // [For MET building] filter the Muon container m_muons, placing selected muons into m_MetMuons
      //
      // For emulated MET (No muons)
      // Make a empty container for invisible electrons
      ConstDataVector<xAOD::MuonContainer> m_EmptyMuons(SG::VIEW_ELEMENTS);
      m_EmptyMuons.clear();
      m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
          xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
          m_met,                         //filling this met container
          m_EmptyMuons.asDataVector(),            //using these metMuons that accepted our cuts
          m_metMap);                     //and this association map
      // Make a container for invisible muons
      ConstDataVector<xAOD::MuonContainer> m_invisibleMuonsForZ(SG::VIEW_ELEMENTS);
      for (const auto& muon : *m_goodMuonForZ) { // C++11 shortcut
        m_invisibleMuonsForZ.push_back( muon );
      }
      // Mark muons invisible
      m_metMaker->markInvisible(m_invisibleMuonsForZ.asDataVector(), m_metMap);

      met::addGhostMuonsToJets(*m_muons, *jetSC);


      // JET
      //-----------------
      //Now time to rebuild jetMet and get the soft term
      //This adds the necessary soft term for both CST and TST
      //these functions create an xAODMissingET object with the given names inside the container

      // For emulated MET marking muons invisible
      m_metMaker->rebuildJetMET("RefJet",          //name of jet met
          "SoftClus",           //name of soft cluster term met
          "PVSoftTrk",          //name of soft track term met
          m_met,       //adding to this new met container
          jetSC,                //using this jet collection to calculate jet met
          m_metCore,   //core met container
          m_metMap,    //with this association map
          true);                //apply jet jvt cut




      /////////////////////////////
      // Soft term uncertainties //
      /////////////////////////////
      if (!m_isData) {
        // Get the track soft term for Zmumu (For emulated MET marking muons invisible)
        xAOD::MissingET* softTrkmet = (*m_met)[softTerm];
        if (m_metSystTool->applyCorrection(*softTrkmet) != CP::CorrectionCode::Ok) {
          Error("execute()", "METSystematicsTool returns Error CorrectionCode");
        }
      }



      ///////////////
      // MET Build //
      ///////////////
      // For emulated MET for Zmumu marking muons invisible
      m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());



      //////////////////////////////////////////////////////////////
      // Fill emulated MET for Zmumu (by marking muons invisible) //
      //////////////////////////////////////////////////////////////
      //emulMET_Zmumu_ex = ((*m_met)["Final"]->mpx());
      //emulMET_Zmumu_ey = ((*m_met)["Final"]->mpy());
      emulMET_Zmumu = ((*m_met)["Final"]->met());
      //emulSumET_Zmumu = ((*m_met)["Final"]->sumet());
      emulMET_Zmumu_phi = ((*m_met)["Final"]->phi());


    } // m_isZmumu








    //-------------------------------------
    // Define Monojet and DiJet Properties
    //-------------------------------------

    // Monojet
    float monojet_pt = 0;
    float monojet_phi = 0;
    float monojet_eta = 0;
    float monojet_rapidity = 0;
    float dPhiMonojetMet = 0;
    float dPhiMonojetMet_Zmumu = 0;
    float dPhiMonojetMet_Wmunu = 0;
    float dPhiMonojetMet_Zee = 0;
    float dPhiMonojetMet_Wenu = 0;
    // Dijet
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
    // SM1
    float sm1jet_pt = 0;
    float sm1jet_phi = 0;
    float sm1jet_eta = 0;
    float sm1jet_rapidity = 0;
    float dPhiSM1jetMet = 0;
    float dPhiSM1jetMet_Zmumu = 0;
    float dPhiSM1jetMet_Zee = 0;


    float goodJet_ht = 0;
    float dPhiJet1Met = 0;
    float dPhiJet2Met = 0;
    float dPhiJet3Met = 0;
    float dPhiJet1Met_Zmumu = 0;
    float dPhiJet2Met_Zmumu = 0;
    float dPhiJet3Met_Zmumu = 0;
    float dPhiJet1Met_Wmunu = 0;
    float dPhiJet2Met_Wmunu = 0;
    float dPhiJet3Met_Wmunu = 0;
    float dPhiJet1Met_Zee = 0;
    float dPhiJet2Met_Zee = 0;
    float dPhiJet3Met_Zee = 0;
    float dPhiJet1Met_Wenu = 0;
    float dPhiJet2Met_Wenu = 0;
    float dPhiJet3Met_Wenu = 0;

    float mjj = 0;
    bool pass_monoJet = false; // Select monoJet
    bool pass_diJet = false; // Select DiJet
    bool pass_CJV = true; // Central Jet Veto (CJV)
    bool pass_sm1Jet = false; // Select monoJet
    bool pass_dPhijetmet = true; // deltaPhi(Jet_i,MET)
    bool pass_dPhijetmet_Zmumu = true; // deltaPhi(Jet_i,MET_Zmumu)
    bool pass_dPhijetmet_Wmunu = true; // deltaPhi(Jet_i,MET_Wmunu)
    bool pass_dPhijetmet_Zee = true; // deltaPhi(Jet_i,MET_Zee)
    bool pass_dPhijetmet_Wenu = true; // deltaPhi(Jet_i,MET_Wenu)
    float dPhiMinjetmet = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)
    float dPhiMinjetmet_Zmumu = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)
    float dPhiMinjetmet_Wmunu = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)
    float dPhiMinjetmet_Zee = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)
    float dPhiMinjetmet_Wenu = 10.; // initialize with 10. to obtain minimum value of deltaPhi(Jet_i,MET)


    ///////////////////////
    // Monojet Selection //
    ///////////////////////
    if (m_goodJet->size() > 0) {

      monojet_pt = m_goodJet->at(0)->pt();
      monojet_phi = m_goodJet->at(0)->phi();
      monojet_eta = m_goodJet->at(0)->eta();
      monojet_rapidity = m_goodJet->at(0)->rapidity();


      // Define Monojet
      if ( monojet_pt > m_monoJetPtCut ){
        if ( fabs(monojet_eta) < m_monoJetEtaCut){
          if ( m_jetCleaningTight->accept( *m_goodJet->at(0) ) ){ //Tight Leading Jet 
            pass_monoJet = true;
            //Info("execute()", "  Leading jet pt = %.2f GeV", monojet_pt * 0.001);
          }
        }
      }

      // deltaPhi(monojet,MET) decision
      // For Znunu
      if (m_isZnunu){
        dPhiMonojetMet = deltaPhi(monojet_phi, MET_phi);
      }
      // For Zmumu
      if (m_isZmumu){
        dPhiMonojetMet_Zmumu = deltaPhi(monojet_phi, emulMET_Zmumu_phi);
      }
      // For Wmunu
      if (m_isWmunu){
        dPhiMonojetMet_Wmunu = deltaPhi(monojet_phi, emulMET_Wmunu_phi);
      }
      // For Zee
      if (m_isZee){
        dPhiMonojetMet_Zee = deltaPhi(monojet_phi, emulMET_Zee_phi);
      }
      // For Wenu
      if (m_isWenu){
        dPhiMonojetMet_Wenu = deltaPhi(monojet_phi, emulMET_Wenu_phi);
      }

    } // MonoJet selection 



    /////////////////////
    // DiJet Selection //
    /////////////////////
    if (m_goodJet->size() > 1) {

      jet1 = m_goodJet->at(0)->p4();
      jet2 = m_goodJet->at(1)->p4();
      jet1_pt = m_goodJet->at(0)->pt();
      jet2_pt = m_goodJet->at(1)->pt();
      jet1_phi = m_goodJet->at(0)->phi();
      jet2_phi = m_goodJet->at(1)->phi();
      jet1_eta = m_goodJet->at(0)->eta();
      jet2_eta = m_goodJet->at(1)->eta();
      jet1_rapidity = m_goodJet->at(0)->rapidity();
      jet2_rapidity = m_goodJet->at(1)->rapidity();
      auto dijet = jet1 + jet2;
      mjj = dijet.M();

      //Info("execute()", "  jet1 = %.2f GeV, jet2 = %.2f GeV", jet1_pt * 0.001, jet2_pt * 0.001);
      //Info("execute()", "  mjj = %.2f GeV", mjj * 0.001);

      // Define Dijet
      if ( jet1_pt > m_diJet1PtCut && jet2_pt > m_diJet2PtCut ){
        if ( fabs(jet1_rapidity) < m_diJetRapCut && fabs(jet2_rapidity) < m_diJetRapCut ){
          if ( m_jetCleaningTight->accept( *m_goodJet->at(0) ) ){ //Tight Leading Jet 
            pass_diJet = true;
          }
        }
      }

      // deltaPhi(Jet1,MET) or deltaPhi(Jet2,MET) decision
      // For Znunu
      if (m_isZnunu){
        dPhiJet1Met = deltaPhi(jet1_phi, MET_phi);
        dPhiJet2Met = deltaPhi(jet2_phi, MET_phi);
      }
      // For Zmumu
      if (m_isZmumu){
        dPhiJet1Met_Zmumu = deltaPhi(jet1_phi, emulMET_Zmumu_phi);
        dPhiJet2Met_Zmumu = deltaPhi(jet2_phi, emulMET_Zmumu_phi);
      }
      // For Wmunu
      if (m_isWmunu){
        dPhiJet1Met_Wmunu = deltaPhi(jet1_phi, emulMET_Wmunu_phi);
        dPhiJet2Met_Wmunu = deltaPhi(jet2_phi, emulMET_Wmunu_phi);
      }
      // For Zee
      if (m_isZee){
        dPhiJet1Met_Zee = deltaPhi(jet1_phi, emulMET_Zee_phi);
        dPhiJet2Met_Zee = deltaPhi(jet2_phi, emulMET_Zee_phi);
      }
      // For Wenu
      if (m_isWenu){
        dPhiJet1Met_Wenu = deltaPhi(jet1_phi, emulMET_Wenu_phi);
        dPhiJet2Met_Wenu = deltaPhi(jet2_phi, emulMET_Wenu_phi);
      }

    } // DiJet selection 



    ///////////////////
    // SM1 Selection //
    ///////////////////
    if (m_goodJet->size() > 0) {

      sm1jet_pt = m_goodJet->at(0)->pt();
      sm1jet_phi = m_goodJet->at(0)->phi();
      sm1jet_eta = m_goodJet->at(0)->eta();
      sm1jet_rapidity = m_goodJet->at(0)->rapidity();

      // Define SM1jet
      if ( sm1jet_pt > m_sm1JetPtCut ){
        if ( fabs(sm1jet_eta) < m_sm1JetEtaCut){ // Threshold of leading jet in SM1
          pass_sm1Jet = true;
          //Info("execute()", "  Leading jet pt = %.2f GeV", sm1jet_pt * 0.001);
        }
      }

      // loop over the jets in the Good Jets Container
      if (pass_sm1Jet) {
        for (const auto& jet : *m_goodJet) {
          uint sm1jet_index = m_goodJet->at(0)->index(); // index for leading jet in SM1
          if (jet->index() != sm1jet_index) { // For subleading jets in SM1
            if (jet->pt() > 30000.) pass_sm1Jet = false; // Subleading jets in SM1 should not be greater than 30GeV
          }
        }
      }

      /* // Test SM1 jets
      if (pass_sm1Jet) {
        Info("execute()", " [SM1] Event # = %llu", eventInfo->eventNumber());
        for (const auto& jet : *m_goodJet) {
          Info("execute()", " [SM1] Jet_pt = %.2f GeV ", jet->pt()*0.001);
        }
      }
      */


      // deltaPhi(sm1jet,MET) decision

      // For Znunu
      if (m_isZnunu){
        dPhiSM1jetMet = deltaPhi(sm1jet_phi, MET_phi);
      }
      // For Zmumu
      if (m_isZmumu){
        dPhiSM1jetMet_Zmumu = deltaPhi(sm1jet_phi, emulMET_Zmumu_phi);
      }
      // For Zee
      if (m_isZee){
        dPhiSM1jetMet_Zee = deltaPhi(sm1jet_phi, emulMET_Zee_phi);
      }

    } // SM1 selection 





    // For jet3
    if (m_goodJet->size() > 2) {
      jet3_pt = m_goodJet->at(2)->pt();
      jet3_phi = m_goodJet->at(2)->phi();
      jet3_eta = m_goodJet->at(2)->eta();
      jet3_rapidity = m_goodJet->at(2)->rapidity();
      // deltaPhi(Jet3,MET)
      dPhiJet3Met = deltaPhi(jet3_phi, MET_phi);
      dPhiJet3Met_Zmumu = deltaPhi(jet3_phi, emulMET_Zmumu_phi);
      dPhiJet3Met_Wmunu = deltaPhi(jet3_phi, emulMET_Wmunu_phi);
      dPhiJet3Met_Zee = deltaPhi(jet3_phi, emulMET_Zee_phi);
      dPhiJet3Met_Wenu = deltaPhi(jet3_phi, emulMET_Wenu_phi);
    }


    // Define deltaPhi(Jet_i,MET) cut and Central Jet Veto (CJV)
    if (m_goodJet->size() > 0) {

      // loop over the jets in the Good Jets Container
      for (const auto& jet : *m_goodJet) {
        float good_jet_pt = jet->pt();
        float good_jet_rapidity = jet->rapidity();
        float good_jet_phi = jet->phi();

        // Calculate dPhi(Jet_i,MET) and dPhi_min(Jet_i,MET)
        if (m_goodJet->at(0) == jet || m_goodJet->at(1) == jet || m_goodJet->at(2) == jet || m_goodJet->at(3) == jet){ // apply cut only to leading jet1, jet2, jet3 and jet4
          // For Znunu
          if (m_isZnunu){
            float dPhijetmet = deltaPhi(good_jet_phi,MET_phi);
            //Info("execute()", " [Znunu] Event # = %llu", eventInfo->eventNumber());
            //Info("execute()", " [Znunu] dPhi = %.2f", dPhijetmet);
            if ( good_jet_pt > 30000. && fabs(good_jet_rapidity) < 4.4 && dPhijetmet < 0.4 ) pass_dPhijetmet = false;
            dPhiMinjetmet = std::min(dPhiMinjetmet, dPhijetmet);
            //Info("execute()", " [Znunu] dPhi_min = %.2f", dPhiMinjetmet);
          }
          // For Zmumu
          if (m_isZmumu){
            float dPhijetmet_Zmumu = deltaPhi(good_jet_phi,emulMET_Zmumu_phi);
            //Info("execute()", " [Zmumu] Event # = %llu", eventInfo->eventNumber());
            //Info("execute()", " [Zmumu] dPhi = %.2f", dPhijetmet_Zmumu);
            if ( good_jet_pt > 30000. && fabs(good_jet_rapidity) < 4.4 && dPhijetmet_Zmumu < 0.4 ) pass_dPhijetmet_Zmumu = false;
            dPhiMinjetmet_Zmumu = std::min(dPhiMinjetmet_Zmumu, dPhijetmet_Zmumu);
            //Info("execute()", " [Zmumu] dPhi_min = %.2f", dPhiMinjetmet_Zmumu);
          }
          // For Wmunu
          if (m_isWmunu){
            float dPhijetmet_Wmunu = deltaPhi(good_jet_phi,emulMET_Wmunu_phi);
            //Info("execute()", " [Wmunu] Event # = %llu", eventInfo->eventNumber());
            //Info("execute()", " [Wmunu] dPhi = %.2f", dPhijetmet_Wmunu);
            if ( good_jet_pt > 30000. && fabs(good_jet_rapidity) < 4.4 && dPhijetmet_Wmunu < 0.4 ) pass_dPhijetmet_Wmunu = false;
            dPhiMinjetmet_Wmunu = std::min(dPhiMinjetmet_Wmunu, dPhijetmet_Wmunu);
            //Info("execute()", " [Wmunu] dPhi_min = %.2f", dPhiMinjetmet_Wmunu);
          }
          // For Zee
          if (m_isZee){
            float dPhijetmet_Zee = deltaPhi(good_jet_phi,emulMET_Zee_phi);
            //Info("execute()", " [Zee] Event # = %llu", eventInfo->eventNumber());
            //Info("execute()", " [Zee] dPhi = %.2f", dPhijetmet_Zee);
            if ( good_jet_pt > 30000. && fabs(good_jet_rapidity) < 4.4 && dPhijetmet_Zee < 0.4 ) pass_dPhijetmet_Zee = false;
            dPhiMinjetmet_Zee = std::min(dPhiMinjetmet_Zee, dPhijetmet_Zee);
            //Info("execute()", " [Zee] dPhi_min = %.2f", dPhiMinjetmet_Zee);
          }
          // For Wenu
          if (m_isWenu){
            float dPhijetmet_Wenu = deltaPhi(good_jet_phi,emulMET_Wenu_phi);
            //Info("execute()", " [Wenu] Event # = %llu", eventInfo->eventNumber());
            //Info("execute()", " [Wenu] dPhi = %.2f", dPhijetmet_Wenu);
            if ( good_jet_pt > 30000. && fabs(good_jet_rapidity) < 4.4 && dPhijetmet_Wenu < 0.4 ) pass_dPhijetmet_Wenu = false;
            dPhiMinjetmet_Wenu = std::min(dPhiMinjetmet_Wenu, dPhijetmet_Wenu);
            //Info("execute()", " [Wenu] dPhi_min = %.2f", dPhiMinjetmet_Wenu);
          }
        }

        // Central Jet Veto (CJV)
        if ( m_goodJet->size() > 2 && pass_diJet ){
          if (m_goodJet->at(0) != jet && m_goodJet->at(1) != jet){
            //cout << "m_goodJet->at(0) = " << m_goodJet->at(0) << " jet = " << jet << endl;
            if (good_jet_pt > m_CJVptCut && fabs(good_jet_rapidity) < m_diJetRapCut) {
              if ( (jet1_rapidity > jet2_rapidity) && (good_jet_rapidity < jet1_rapidity && good_jet_rapidity > jet2_rapidity)){
                pass_CJV = false;
              }
              if ( (jet1_rapidity < jet2_rapidity) && (good_jet_rapidity > jet1_rapidity && good_jet_rapidity < jet2_rapidity)){
                pass_CJV = false;
              }
              /* //Valentinos' way (same result as mine)
              float rapLow  = std::min(jet1_rapidity, jet2_rapidity);
              float rapHigh = std::max(jet1_rapidity, jet2_rapidity);
              if (good_jet_rapidity > rapLow && good_jet_rapidity < rapHigh) pass_CJV = false;
              */
            }
          }
        }

        //Info("execute()", "  Znunu Signal Jet pt = %.2f GeV, eta = %.2f", good_pt_jet * 0.001, good_eta_jet);
        goodJet_ht += good_jet_pt;
      } // Jet loop

    } // End deltaPhi(Jet_i,MET) cut and Central Jet Veto (CJV)




    //----------------------------------
    // Define Zmumu and Wmunu Selection
    //----------------------------------

    // Zmumu
    float mll_muon = 0.;
    float muon1_pt = 0.;
    float muon2_pt = 0.;
    float muon1_charge = 0.;
    float muon2_charge = 0.;
    bool pass_dimuonPtCut = false; // dimuon pT cut
    bool pass_OSmuon = false; // Opposite sign change muon
    bool pass_SSmuon = false; // Same sign change muon
    int numExtra = 0;

    if (m_isZmumu) {

      // For Zmumu Selection
      if (m_goodMuonForZ->size() > 1) {

        TLorentzVector muon1 = m_goodMuonForZ->at(0)->p4();
        TLorentzVector muon2 = m_goodMuonForZ->at(1)->p4();
        muon1_pt = m_goodMuonForZ->at(0)->pt();
        muon2_pt = m_goodMuonForZ->at(1)->pt();
        muon1_charge = m_goodMuonForZ->at(0)->charge();
        muon2_charge = m_goodMuonForZ->at(1)->charge();
        auto Zmass_muon = muon1 + muon2;
        mll_muon = Zmass_muon.M();

        //Info("execute()", "  muon1 = %.2f GeV, muon2 = %.2f GeV", muon1_pt * 0.001, muon2_pt * 0.001);
        //Info("execute()", "  mll (Zmumu) = %.2f GeV", mll_muon * 0.001);

        if ( muon1_pt >  m_LeadLepPtCut && muon2_pt > m_SubLeadLepPtCut ) pass_dimuonPtCut = true;
        if ( muon1_charge * muon2_charge < 0 ) pass_OSmuon = true;
        if ( muon1_charge * muon2_charge > 0 ) pass_SSmuon = true;
        //Info("execute()", "  muon1 charge = %f, muon2 charge = %f, pass_OSmuon = %d, pass_SSmuon = %d", muon1_charge, muon2_charge, pass_OSmuon, pass_SSmuon);

        uint index1 = m_goodMuonForZ->at(0)->index();
        uint index2 = m_goodMuonForZ->at(1)->index();
        for (const auto &muon : *m_goodMuon) {
          if (muon->index() != index1 && muon->index() != index2) numExtra++;
          //std::cout << muon->index() << " " << index1 << " " << index2 << std::endl;
        }

      } // Zmumu selection loop

    }


    // Wmunu
    bool pass_Wmunu = false;
    float mT_muon = 0.;

    if (m_isWmunu) {
      // Wmunu Selection
      if (m_goodMuon->size() == 1) {
        float muon_pt = m_goodMuon->at(0)->pt();
        float muon_phi = m_goodMuon->at(0)->phi();
        mT_muon = TMath::Sqrt( 2. * muon_pt * MET * ( 1. - TMath::Cos(muon_phi - MET_phi) ) );

        if ( muon_pt > 25000. ){
          pass_Wmunu = true;
        }

      } //Wmunu Selection loop

    }


    //-------------------------------
    // Define Zee and Wenu Selection
    // ------------------------------

    // Zee
    float mll_electron = 0.;
    float electron1_pt = 0.;
    float electron2_pt = 0.;
    float electron1_charge = 0.;
    float electron2_charge = 0.;
    bool pass_dielectronPtCut = false; // dielectron pT cut
    bool pass_OSelectron = false; // Opposite sign change electron
    bool pass_SSelectron = false; // Same sign change electron

    if (m_isZee){

      // Zee Selection
      if (m_goodElectron->size() > 1) {

        TLorentzVector electron1 = m_goodElectron->at(0)->p4();
        TLorentzVector electron2 = m_goodElectron->at(1)->p4();
        electron1_pt = m_goodElectron->at(0)->pt();
        electron2_pt = m_goodElectron->at(1)->pt();
        electron1_charge = m_goodElectron->at(0)->charge();
        electron2_charge = m_goodElectron->at(1)->charge();
        auto Zmass_electron = electron1 + electron2;
        mll_electron = Zmass_electron.M();

        //Info("execute()", "  electron1 = %.2f GeV, electron2 = %.2f GeV", electron1_pt * 0.001, electron2_pt * 0.001);
        //Info("execute()", "  mll (Zee) = %.2f GeV", mll_electron * 0.001);


        if ( electron1_pt >  m_LeadLepPtCut && electron2_pt > m_SubLeadLepPtCut ) pass_dielectronPtCut = true;
        if ( electron1_charge * electron2_charge < 0 ) pass_OSelectron = true;
        if ( electron1_charge * electron2_charge > 0 ) pass_SSelectron = true;

      } // Zee selection loop

    }


    // Wenu
    bool pass_Wenu = false;
    float mT_electron = 0.;

    if (m_isWenu){

      // Wenu Selection
      if (m_goodElectron->size() == 1) {
        float electron_pt = m_goodElectron->at(0)->pt();
        float electron_phi = m_goodElectron->at(0)->phi();
        mT_electron = TMath::Sqrt( 2. * electron_pt * MET * ( 1. - TMath::Cos(electron_phi - MET_phi) ) );

        if ( electron_pt > 25000. ){
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
    //Info("execute()", "  The number of Isolated track counted = %i (N_SignalMuon = %lu, N_SignalElec = %lu)", Nisotrk, m_goodMuon->size(), m_goodElectron->size() );
    }
    */



    //-----------
    // VBF study 
    //-----------

    //-------------------------------
    // Z -> nunu + JET EVENT SELECTION
    //-------------------------------

    if (m_isZnunu){
      h_channel = "h_znunu_";
      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]MET Trigger");
        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[5]+=1;
        if ( MET > m_metCut ) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]MET cut");
          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[6]+=1;
          if (m_goodElectron->size() == 0) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]Electron Veto");
            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[7]+=1;
            if ( m_goodMuon->size() == 0) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]Muon Veto");
              if (sysName == "" && m_useArrayCutflow) m_eventCutflow[8]+=1;
              if (m_goodTau->size() == 0) {
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]Tau Veto");
                if (sysName == "" && m_useArrayCutflow) m_eventCutflow[9]+=1;
                if ( m_goodJet->size() > 0 ) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu]At least One Jets");
                  if (sysName == "" && m_useArrayCutflow) m_eventCutflow[10]+=1;

                  ////////////////////////
                  // MonoJet phasespace //
                  ////////////////////////
                  if ( pass_monoJet ) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu, monojet]MonoJet");
                    if ( pass_dPhijetmet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu, monojet]dPhi(jet_i,MET) cut");

                      // Fill histogram
                      // For Ratio plot (Blind MET and Mjj for Ratio)
                      if (MET < m_METblindcut) {
                        hMap1D["Znunu_MET_mono"+sysName]->Fill(MET * 0.001, mcEventWeight);
                      }
                      // For publication
                      hMap1D[h_channel+"monojet_met"+sysName]->Fill(MET * 0.001, mcEventWeight);
                      // Average Interaction
                      hMap1D[h_channel+"monojet_avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight);
                      if (sysName == ""){
                        // Jets
                        hMap1D[h_channel+"monojet_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight);
                        hMap1D[h_channel+"monojet_jet_pt"+sysName]->Fill(monojet_pt * 0.001, mcEventWeight);
                        hMap1D[h_channel+"monojet_jet_phi"+sysName]->Fill(monojet_phi, mcEventWeight);
                        hMap1D[h_channel+"monojet_jet_eta"+sysName]->Fill(monojet_eta, mcEventWeight);
                        hMap1D[h_channel+"monojet_jet_rap"+sysName]->Fill(monojet_rapidity, mcEventWeight);
                        hMap1D[h_channel+"monojet_dPhimetjet"+sysName]->Fill(dPhiMonojetMet, mcEventWeight);
                        hMap1D[h_channel+"monojet_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet, mcEventWeight);
                      }

                    } // pass dPhijetmet
                  } // pass monojet

                  ////////////////////
                  // VBF phasespace //
                  ////////////////////
                  if ( pass_diJet ) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu, VBF]DiJet");
                    if (sysName == "" && m_useArrayCutflow) m_eventCutflow[11]+=1;
                    if ( mjj > m_mjjCut ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu, VBF]mjj cut");
                      if (sysName == "" && m_useArrayCutflow) m_eventCutflow[12]+=1;
                      if ( pass_CJV ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu, VBF]CJV cut");
                        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[13]+=1;
                        if ( pass_dPhijetmet ) {
                          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Znunu, VBF]dPhi(jet_i,MET) cut");
                          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[14]+=1;
                          // Fill histogram
                          // For Ratio plot (Blind MET and Mjj for Ratio)
                          if (MET < m_METblindcut && mjj < m_Mjjblindcut) {
                            hMap1D["Znunu_MET_search"+sysName]->Fill(MET * 0.001, mcEventWeight);
                            hMap1D["Znunu_Mjj_search"+sysName]->Fill(mjj * 0.001, mcEventWeight);
                            hMap1D["Znunu_DeltaPhiAll"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight);
                          }
                          // For publication
                          hMap1D[h_channel+"vbf_met"+sysName]->Fill(MET * 0.001, mcEventWeight);
                          hMap1D[h_channel+"vbf_mjj"+sysName]->Fill(mjj * 0.001, mcEventWeight);
                          hMap1D[h_channel+"vbf_dPhijj"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight);
                          // Average Interaction
                          hMap1D[h_channel+"vbf_avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight);
                          if (sysName == ""){
                            // Jets
                            hMap1D[h_channel+"vbf_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight);
                            hMap1D[h_channel+"vbf_jet1_pt"+sysName]->Fill(jet1_pt * 0.001, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet2_pt"+sysName]->Fill(jet2_pt * 0.001, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet1_phi"+sysName]->Fill(jet1_phi, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet2_phi"+sysName]->Fill(jet2_phi, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet1_eta"+sysName]->Fill(jet1_eta, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet2_eta"+sysName]->Fill(jet2_eta, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet1_rap"+sysName]->Fill(jet1_rapidity, mcEventWeight);
                            hMap1D[h_channel+"vbf_jet2_rap"+sysName]->Fill(jet2_rapidity, mcEventWeight);
                            hMap1D[h_channel+"vbf_dRjj"+sysName]->Fill(deltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight);
                            hMap1D[h_channel+"vbf_dPhimetj1"+sysName]->Fill(dPhiJet1Met, mcEventWeight);
                            hMap1D[h_channel+"vbf_dPhimetj2"+sysName]->Fill(dPhiJet2Met, mcEventWeight);
                            hMap1D[h_channel+"vbf_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet, mcEventWeight);
                            // For jet3
                            if (m_goodJet->size() > 2){
                              hMap1D[h_channel+"vbf_jet3_pt"+sysName]->Fill(jet3_pt * 0.001, mcEventWeight);
                              hMap1D[h_channel+"vbf_jet3_phi"+sysName]->Fill(jet3_phi, mcEventWeight);
                              hMap1D[h_channel+"vbf_jet3_eta"+sysName]->Fill(jet3_eta, mcEventWeight);
                              hMap1D[h_channel+"vbf_jet3_rap"+sysName]->Fill(jet3_rapidity, mcEventWeight);
                              hMap1D[h_channel+"vbf_dPhimetj3"+sysName]->Fill(dPhiJet3Met, mcEventWeight);
                            }
                          }

                        } // pass diPhijetMET
                      } // pass CJV
                    } // mjj Cut
                  } // pass diJet

                } // at least 1 jet
              } // Tau veto
            } // Muon veto
          } // Electron veto
        } // MET cut
      } // HLT_xe70
    } // m_isZnunu


          /*
             Info("execute()", "=====================================");
             Info("execute()", " Event # = %llu", eventInfo->eventNumber());
             Info("execute()", " Good Event number = %i", m_eventCutflow[6]);
             Info("execute()", " MET = %.3f GeV", MET * 0.001);
             Info("execute()", " RefElectron = %.3f GeV", ((*m_met)["RefElectron"]->met()) * 0.001);
             Info("execute()", " RefPhoton = %.3f GeV", ((*m_met)["RefPhoton"]->met()) * 0.001);
             Info("execute()", " RefTau = %.3f GeV", ((*m_met)["RefTau"]->met()) * 0.001);
             Info("execute()", " RefMuon = %.3f GeV", ((*m_met)["RefMuon"]->met()) * 0.001);
             Info("execute()", " RefJet = %.3f GeV", ((*m_met)["RefJet"]->met()) * 0.001);
             Info("execute()", " SoftClus = %.3f GeV", ((*m_met)["SoftClus"]->met()) * 0.001);
             Info("execute()", " PVSoftTrk = %.3f GeV", ((*m_met)["PVSoftTrk"]->met()) * 0.001);
             Info("execute()", " # of good jets = %lu", m_goodJet->size());
             if (m_goodJet->size() > 0){
             int jetCount = 0;
             for (const auto& jet : *m_goodJet) {
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
             Info("execute()", " electron eta = %.3f", electron->eta());
             Info("execute()", " electron phi = %.3f", electron->phi());
             }
             }
             if (m_goodMuon->size() > 0){
             int muCount = 0;
             for (const auto& muon : *m_goodMuon) {
             muCount++;
             Info("execute()", " muon # : %i", muCount);
             Info("execute()", " muon pt = %.3f GeV", muon->pt() * 0.001);
             Info("execute()", " muon eta = %.3f", muon->eta());
             Info("execute()", " muon phi = %.3f", muon->phi());
             }
             }
             if (m_goodTau->size() > 0){
             int tauCount = 0;
             for (const auto& tau : *m_goodTau) {
             tauCount++;
             Info("execute()", " tau # : %i", tauCount);
             Info("execute()", " tau pt = %.3f GeV", tau->pt() * 0.001);
             Info("execute()", " tau eta = %.3f", tau->eta());
             Info("execute()", " tau phi = %.3f", tau->phi());
             }
             }
             */



    //---------------------------------
    // Z -> mumu + JET EVENT SELECTION
    //---------------------------------

    if (m_isZmumu){
      h_channel = "h_zmumu_";
      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]MET Trigger");
        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[16]+=1;
        if ( emulMET_Zmumu > m_metCut ) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]MET cut");
          if (sysName == "" && m_useArrayCutflow) m_eventCutflow[17]+=1;
          if (m_goodElectron->size() == 0) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]Electron Veto");
            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[18]+=1;
            if ( m_goodMuonForZ->size() > 1) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]At least Two Muons");
              if (sysName == "" && m_useArrayCutflow) m_eventCutflow[19]+=1;
              if (m_goodTau->size() == 0) {
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]Tau Veto");
                if (sysName == "" && m_useArrayCutflow) m_eventCutflow[20]+=1;
                if ( pass_dimuonPtCut && pass_OSmuon && numExtra == 0 && mll_muon > m_mllMin && mll_muon < m_mllMax ){
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]mll cut");
                  if (sysName == "" && m_useArrayCutflow) m_eventCutflow[21]+=1;
                  if ( m_goodJet->size() > 0 ) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]At least One Jets");
                    if (sysName == "" && m_useArrayCutflow) m_eventCutflow[22]+=1;

                    ////////////////////////
                    // MonoJet phasespace //
                    ////////////////////////
                    if ( pass_monoJet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, monojet]MonoJet");
                      if ( pass_dPhijetmet_Zmumu ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, monojet]dPhi(jet_i,MET) cut");

                        // Calculate muon SF for Zmumu
                        float mcEventWeight_Zmumu = 1.;
                        if (!m_isData) {
                          double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                          //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                          mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                        }

                        // Fill histogram
                        // For Ratio plot (Blind MET and Mjj for Ratio)
                        if (emulMET_Zmumu < m_METblindcut) {
                          hMap1D["Zmumu_MET_mono"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                        }
                        // For publication
                        hMap1D[h_channel+"monojet_met_emulmet"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                        // Average Interaction
                        hMap1D[h_channel+"monojet_avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight_Zmumu);
                        if (sysName == ""){
                          // Jets
                          hMap1D[h_channel+"monojet_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_jet_pt"+sysName]->Fill(monojet_pt * 0.001, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_jet_phi"+sysName]->Fill(monojet_phi, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_jet_eta"+sysName]->Fill(monojet_eta, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_jet_rap"+sysName]->Fill(monojet_rapidity, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_dPhimetjet"+sysName]->Fill(dPhiMonojetMet_Zmumu, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_Zmumu, mcEventWeight_Zmumu);
                          // Leptons
                          hMap1D[h_channel+"monojet_lepton1_pt"+sysName]->Fill(m_goodMuonForZ->at(0)->pt() * 0.001, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_lepton2_pt"+sysName]->Fill(m_goodMuonForZ->at(1)->pt() * 0.001, mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_lepton1_phi"+sysName]->Fill(m_goodMuonForZ->at(0)->phi(), mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_lepton2_phi"+sysName]->Fill(m_goodMuonForZ->at(1)->phi(), mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_lepton1_eta"+sysName]->Fill(m_goodMuonForZ->at(0)->eta(), mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_lepton2_eta"+sysName]->Fill(m_goodMuonForZ->at(1)->eta(), mcEventWeight_Zmumu);
                          hMap1D[h_channel+"monojet_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                        }

                      } // pass dPhijetmet_Zmumu
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
                          if ( pass_dPhijetmet_Zmumu ) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu, VBF]dPhi(jet_i,MET) cut");
                            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[26]+=1;

                            // Calculate muon SF for Zmumu
                            float mcEventWeight_Zmumu = 1.;
                            if (!m_isData) {
                              double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                              //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                              mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                            }

                            // Fill histogram
                            // For Ratio plot (Blind MET and Mjj for Ratio)
                            if (emulMET_Zmumu < m_METblindcut && mjj < m_Mjjblindcut) {
                              hMap1D["Zmumu_MET_search"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                              hMap1D["Zmumu_Mjj_search"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zmumu);
                              hMap1D["Zmumu_DeltaPhiAll"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                            }
                            // For publication
                            hMap1D[h_channel+"vbf_met_emulmet"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_mjj"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zmumu);
                            hMap1D[h_channel+"vbf_dPhijj"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                            // Average Interaction
                            hMap1D[h_channel+"vbf_avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight_Zmumu);
                            if (sysName == ""){
                              // Jets
                              hMap1D[h_channel+"vbf_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet1_pt"+sysName]->Fill(jet1_pt * 0.001, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet2_pt"+sysName]->Fill(jet2_pt * 0.001, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet1_phi"+sysName]->Fill(jet1_phi, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet2_phi"+sysName]->Fill(jet2_phi, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet1_eta"+sysName]->Fill(jet1_eta, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet2_eta"+sysName]->Fill(jet2_eta, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet1_rap"+sysName]->Fill(jet1_rapidity, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_jet2_rap"+sysName]->Fill(jet2_rapidity, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_dRjj"+sysName]->Fill(deltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_dPhimetj1"+sysName]->Fill(dPhiJet1Met_Zmumu, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_dPhimetj2"+sysName]->Fill(dPhiJet2Met_Zmumu, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_Zmumu, mcEventWeight_Zmumu);
                              // For jet3
                              if (m_goodJet->size() > 2){
                                hMap1D[h_channel+"vbf_jet3_pt"+sysName]->Fill(jet3_pt * 0.001, mcEventWeight_Zmumu);
                                hMap1D[h_channel+"vbf_jet3_phi"+sysName]->Fill(jet3_phi, mcEventWeight_Zmumu);
                                hMap1D[h_channel+"vbf_jet3_eta"+sysName]->Fill(jet3_eta, mcEventWeight_Zmumu);
                                hMap1D[h_channel+"vbf_jet3_rap"+sysName]->Fill(jet3_rapidity, mcEventWeight_Zmumu);
                                hMap1D[h_channel+"vbf_dPhimetj3"+sysName]->Fill(dPhiJet3Met_Zmumu, mcEventWeight_Zmumu);
                              }
                              // Leptons
                              hMap1D[h_channel+"vbf_lepton1_pt"+sysName]->Fill(m_goodMuonForZ->at(0)->pt() * 0.001, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_lepton2_pt"+sysName]->Fill(m_goodMuonForZ->at(1)->pt() * 0.001, mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_lepton1_phi"+sysName]->Fill(m_goodMuonForZ->at(0)->phi(), mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_lepton2_phi"+sysName]->Fill(m_goodMuonForZ->at(1)->phi(), mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_lepton1_eta"+sysName]->Fill(m_goodMuonForZ->at(0)->eta(), mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_lepton2_eta"+sysName]->Fill(m_goodMuonForZ->at(1)->eta(), mcEventWeight_Zmumu);
                              hMap1D[h_channel+"vbf_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                            }


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
        if ( emulMET_Wmunu > m_metCut ) {
          if (m_goodElectron->size() == 0) {
            if ( m_goodMuon->size() > 0 ) {
              if (m_goodTau->size() == 0) {
                if ( pass_Wmunu && m_goodMuon->size() == 1 && mT_muon > 30000. && mT_muon < 100000. ){
                  if ( m_goodJet->size() > 1 ) {
                    if ( pass_diJet ) {
                      if ( mjj > m_mjjCut ) {
                        if ( pass_dPhijetmet_Wmunu ) {
                          if ( pass_CJV ) {
                            // Calculate muon SF for Wmunu
                            float mcEventWeight_Wmunu = 1.;
                            if (!m_isData) {
                              double totalMuonSF_Wmunu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSF, m_ttvaSF);
                              //Info("execute()", " Wmunu Total Muon SF = %.3f ", totalMuonSF_Wmunu);
                              mcEventWeight_Wmunu = mcEventWeight * totalMuonSF_Wmunu;
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
      }
    }




    //-------------------------------
    // Z -> ee + JET EVENT SELECTION
    //-------------------------------

    if (m_isZee){
      h_channel = "h_zee_";
      if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]Electron Trigger");
        if (sysName == "" && m_useArrayCutflow) m_eventCutflow[28]+=1;
        if ( emulMET_Zee > m_metCut ) {
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
                if ( pass_dielectronPtCut && pass_OSelectron && m_goodElectron->size() == 2 && mll_electron > m_mllMin && mll_electron < m_mllMax ) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]mll cut");
                  if (sysName == "" && m_useArrayCutflow) m_eventCutflow[33]+=1;
                  if ( m_goodJet->size() > 0 ) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee]At least One Jets");
                    if (sysName == "" && m_useArrayCutflow) m_eventCutflow[34]+=1;

                    ////////////////////////
                    // MonoJet phasespace //
                    ////////////////////////
                    if ( pass_monoJet ) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, monojet]MonoJet");
                      if ( pass_dPhijetmet_Zee ) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, monojet]dPhi(jet_i,MET) cut");

                        // Calculate electron SF
                        float mcEventWeight_Zee = 1.;
                        if (!m_isData) {
                          float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                          //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                          mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                        }

                        // Fill histogram
                        // For Ratio plot (Blind MET and Mjj)
                        if (emulMET_Zee < m_METblindcut) {
                          hMap1D["Zee_MET_mono"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                        }
                        // For publication
                        hMap1D[h_channel+"monojet_met_emulmet"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                        // Average Interaction
                        hMap1D[h_channel+"monojet_avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight_Zee);
                        if (sysName == ""){
                          // Jets
                          hMap1D[h_channel+"monojet_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_jet_pt"+sysName]->Fill(monojet_pt * 0.001, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_jet_phi"+sysName]->Fill(monojet_phi, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_jet_eta"+sysName]->Fill(monojet_eta, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_jet_rap"+sysName]->Fill(monojet_rapidity, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_dPhimetjet"+sysName]->Fill(dPhiMonojetMet_Zee, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_Zee, mcEventWeight_Zee);
                          // Leptons
                          hMap1D[h_channel+"monojet_lepton1_pt"+sysName]->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_lepton2_pt"+sysName]->Fill(m_goodElectron->at(1)->pt() * 0.001, mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_lepton1_phi"+sysName]->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_lepton2_phi"+sysName]->Fill(m_goodElectron->at(1)->phi(), mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_lepton1_eta"+sysName]->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_lepton2_eta"+sysName]->Fill(m_goodElectron->at(1)->eta(), mcEventWeight_Zee);
                          hMap1D[h_channel+"monojet_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                        }

                      } // pass dPhijetmet_Zee
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
                          if ( pass_dPhijetmet_Zee ) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zee, VBF]dPhi(jet_i,MET) cut");
                            if (sysName == "" && m_useArrayCutflow) m_eventCutflow[38]+=1;

                            // Calculate electron SF
                            float mcEventWeight_Zee = 1.;
                            if (!m_isData) {
                              float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                              //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                              mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                            }

                            // Fill histogram
                            // For Ratio plot (Blind MET and Mjj)
                            if (emulMET_Zee < m_METblindcut && mjj < m_Mjjblindcut) {
                              hMap1D["Zee_MET_search"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                              hMap1D["Zee_Mjj_search"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zee);
                              hMap1D["Zee_DeltaPhiAll"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                            }
                            // For publication
                            hMap1D[h_channel+"vbf_met_emulmet"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_mjj"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zee);
                            hMap1D[h_channel+"vbf_dPhijj"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                            // Average Interaction
                            hMap1D[h_channel+"vbf_avg_interaction"+sysName]->Fill(m_AverageInteractionsPerCrossing, mcEventWeight_Zee);
                            if (sysName == ""){
                              // Jets
                              hMap1D[h_channel+"vbf_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet1_pt"+sysName]->Fill(jet1_pt * 0.001, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet2_pt"+sysName]->Fill(jet2_pt * 0.001, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet1_phi"+sysName]->Fill(jet1_phi, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet2_phi"+sysName]->Fill(jet2_phi, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet1_eta"+sysName]->Fill(jet1_eta, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet2_eta"+sysName]->Fill(jet2_eta, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet1_rap"+sysName]->Fill(jet1_rapidity, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_jet2_rap"+sysName]->Fill(jet2_rapidity, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_dRjj"+sysName]->Fill(deltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_dPhimetj1"+sysName]->Fill(dPhiJet1Met_Zee, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_dPhimetj2"+sysName]->Fill(dPhiJet2Met_Zee, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_Zee, mcEventWeight_Zee);
                              // For jet3
                              if (m_goodJet->size() > 2){
                                hMap1D[h_channel+"vbf_jet3_pt"+sysName]->Fill(jet3_pt * 0.001, mcEventWeight_Zee);
                                hMap1D[h_channel+"vbf_jet3_phi"+sysName]->Fill(jet3_phi, mcEventWeight_Zee);
                                hMap1D[h_channel+"vbf_jet3_eta"+sysName]->Fill(jet3_eta, mcEventWeight_Zee);
                                hMap1D[h_channel+"vbf_jet3_rap"+sysName]->Fill(jet3_rapidity, mcEventWeight_Zee);
                                hMap1D[h_channel+"vbf_dPhimetj3"+sysName]->Fill(dPhiJet3Met_Zee, mcEventWeight_Zee);
                              }
                              // Leptons
                              hMap1D[h_channel+"vbf_lepton1_pt"+sysName]->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_lepton2_pt"+sysName]->Fill(m_goodElectron->at(1)->pt() * 0.001, mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_lepton1_phi"+sysName]->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_lepton2_phi"+sysName]->Fill(m_goodElectron->at(1)->phi(), mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_lepton1_eta"+sysName]->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_lepton2_eta"+sysName]->Fill(m_goodElectron->at(1)->eta(), mcEventWeight_Zee);
                              hMap1D[h_channel+"vbf_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                            }

                          } // pass dPhijetmet_Zee
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
        if ( emulMET_Wenu > m_metCut ) {
          if (m_goodElectron->size() > 0) {
            if ( m_goodMuon->size() == 0 ) {
              if (m_goodTau->size() == 0) {
                if ( pass_Wenu && m_goodElectron->size() == 1 && mT_electron > 30000. && mT_electron < 100000. ){
                  if ( m_goodJet->size() > 1 ) {
                    if ( pass_diJet ) {
                      if ( pass_dPhijetmet_Wenu){
                        if ( mjj > m_mjjCut ) {
                          if ( pass_CJV ) {

                            // Calculate electron SF
                            float mcEventWeight_Wenu = 1.;
                            if (!m_isData) {
                              double totalElectronSF_Wenu = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                              //Info("execute()", " Wenu Total Electron SF = %.3f ", totalElectronSF_Wenu);
                              mcEventWeight_Wenu = mcEventWeight * totalElectronSF_Wenu;
                            }
                            /*
                            // Fill histogram
                            // MET
                            h_wenu_met->Fill(MET * 0.001, mcEventWeight_Wenu);
                            h_wenu_emulmet_Wenu->Fill(emulMET_Wenu * 0.001, mcEventWeight_Wenu);
                            // Jets
                            h_wenu_njet->Fill(m_goodJet->size(), mcEventWeight_Wenu);
                            h_wenu_jet1_pt->Fill(jet1_pt * 0.001, mcEventWeight_Wenu);
                            h_wenu_jet2_pt->Fill(jet2_pt * 0.001, mcEventWeight_Wenu);
                            h_wenu_jet1_phi->Fill(jet1_phi, mcEventWeight_Wenu);
                            h_wenu_jet2_phi->Fill(jet2_phi, mcEventWeight_Wenu);
                            h_wenu_jet1_eta->Fill(jet1_eta, mcEventWeight_Wenu);
                            h_wenu_jet2_eta->Fill(jet2_eta, mcEventWeight_Wenu);
                            h_wenu_jet1_rap->Fill(jet1_rapidity, mcEventWeight_Wenu);
                            h_wenu_jet2_rap->Fill(jet2_rapidity, mcEventWeight_Wenu);
                            h_wenu_mjj->Fill(mjj * 0.001, mcEventWeight_Wenu);
                            h_wenu_dPhijj->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Wenu);
                            h_wenu_dRjj->Fill(deltaR(jet1_eta, jet2_eta, jet1_phi, jet2_phi), mcEventWeight_Wenu);
                            h_wenu_dPhimetj1->Fill(dPhiJet1Met_Wenu, mcEventWeight_Wenu);
                            h_wenu_dPhimetj2->Fill(dPhiJet2Met_Wenu, mcEventWeight_Wenu);
                            h_wenu_dPhiMinmetjet->Fill(dPhiMinjetmet_Wenu, mcEventWeight_Wenu);
                            // For jet3
                            if (m_goodJet->size() > 2){
                            h_wenu_jet3_pt->Fill(jet3_pt * 0.001, mcEventWeight_Wenu);
                            h_wenu_jet3_phi->Fill(jet3_phi, mcEventWeight_Wenu);
                            h_wenu_jet3_eta->Fill(jet3_eta, mcEventWeight_Wenu);
                            h_wenu_jet3_rap->Fill(jet3_rapidity, mcEventWeight_Wenu);
                            h_wenu_dPhimetj3->Fill(dPhiJet3Met_Wenu, mcEventWeight_Wenu);
                            }
                            // Leptons
                            h_wenu_electron_pt->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Wenu);
                            h_wenu_electron_phi->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Wenu);
                            h_wenu_electron_eta->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Wenu);
                            h_wenu_mT->Fill(mT_electron * 0.001, mcEventWeight_Wenu);
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
    }




    //----------------------------------------
    // Z -> mumu + JET MET Trigger Efficiency
    //----------------------------------------

    ////////////////////
    // VBF phasespace //
    ////////////////////

    if (m_isZmumu && sysName == ""){
      h_channel = "h_zmumu_";

      if ( m_trigDecisionTool->isPassed("HLT_mu20_iloose_L1MU15") || m_trigDecisionTool->isPassed("HLT_mu50") ) { // pass muon trigger to avoid bias
        if ( m_goodMuonForZ->size() > 1 && m_goodElectron->size() == 0 && m_goodTau->size() == 0 ) { // Letopn veto
          if (numExtra == 0 && pass_dimuonPtCut && pass_OSmuon && mll_muon > m_mllMin && mll_muon < m_mllMax) {
            if ( m_goodJet->size() > 0 ) {
              if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zmumu) {
                // Fill histogram
                // MET Trigger efficiency (for turn-on curve)
                hMap1D[h_channel+"vbf_eff_study_met_emulmet"+sysName]->Fill(emulMET_Zmumu * 0.001, 1.);
                if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                  hMap1D[h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70"+sysName]->Fill(emulMET_Zmumu * 0.001, 1.);
                }
                if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                  hMap1D[h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70_tclcw"+sysName]->Fill(emulMET_Zmumu * 0.001, 1.);
                }
                // MET Trigger efficiency for mjj and dPhi(j1,j2)
                // MET > 130 GeV
                if ( emulMET_Zmumu > 130000. ) {
                  hMap1D[h_channel+"vbf_eff_study_mjj_met130"+sysName]->Fill(mjj * 0.001, 1.);
                  hMap1D[h_channel+"vbf_eff_study_dPhijj_met130"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                    hMap1D[h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70"+sysName]->Fill(mjj * 0.001, 1.);
                    hMap1D[h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  }
                  if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                    hMap1D[h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70_tclcw"+sysName]->Fill(mjj * 0.001, 1.);
                    hMap1D[h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70_tclcw"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  }
                }
                // MET > 150 GeV
                if ( emulMET_Zmumu > 150000. ) {
                  hMap1D[h_channel+"vbf_eff_study_mjj_met150"+sysName]->Fill(mjj * 0.001, 1.);
                  hMap1D[h_channel+"vbf_eff_study_dPhijj_met150"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                    hMap1D[h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70"+sysName]->Fill(mjj * 0.001, 1.);
                    hMap1D[h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  }
                  if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                    hMap1D[h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70_tclcw"+sysName]->Fill(mjj * 0.001, 1.);
                    hMap1D[h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70_tclcw"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  }
                }
                // MET > 200 GeV
                if ( emulMET_Zmumu > 200000. ) {
                  hMap1D[h_channel+"vbf_eff_study_mjj_met200"+sysName]->Fill(mjj * 0.001, 1.);
                  hMap1D[h_channel+"vbf_eff_study_dPhijj_met200"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                    hMap1D[h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70"+sysName]->Fill(mjj * 0.001, 1.);
                    hMap1D[h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  }
                  if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                    hMap1D[h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70_tclcw"+sysName]->Fill(mjj * 0.001, 1.);
                    hMap1D[h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70_tclcw"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                  }
                }

              } // VBF cut
            } // At least 1 jet
          } // Dimuon cut
        } // Lepton veto
      } // Single muon trigger

    } // end cutflow



    //-----------------------------------------
    // W -> munu + JET MET Trigger Efficiency
    //-----------------------------------------

    if (m_isWmunu && sysName == ""){
      h_channel = "h_wmunu_";

      if ( m_trigDecisionTool->isPassed("HLT_mu20_iloose_L1MU15") || m_trigDecisionTool->isPassed("HLT_mu50") ) { // pass muon trigger to avoid bias
        if (m_goodElectron->size() == 0) {
          if ( m_goodMuon->size() > 0 ) {
            if (m_goodTau->size() == 0) {
              if ( pass_Wmunu && m_goodMuon->size() == 1 && mT_muon > 30000. && mT_muon < 100000. ){
                if ( m_goodJet->size() > 1 ) {
                  if ( pass_diJet ) {
                    if ( mjj > m_mjjCut ) {
                      if ( pass_dPhijetmet_Wmunu ) {
                        if ( pass_CJV ) {
                          // Fill histogram
                          // MET Trigger efficiency (for turn-on curve)
                          hMap1D[h_channel+"vbf_eff_study_met_emulmet"+sysName]->Fill(emulMET_Zmumu * 0.001, 1.);
                          if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                            hMap1D[h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70"+sysName]->Fill(emulMET_Zmumu * 0.001, 1.);
                          }
                          if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                            hMap1D[h_channel+"vbf_eff_study_met_emulmet_pass_HLT_xe70_tclcw"+sysName]->Fill(emulMET_Zmumu * 0.001, 1.);
                          }
                          // MET Trigger efficiency for mjj and dPhi(j1,j2)
                          // MET > 130 GeV
                          if ( emulMET_Wmunu > 130000. ) {
                            hMap1D[h_channel+"vbf_eff_study_mjj_met130"+sysName]->Fill(mjj * 0.001, 1.);
                            hMap1D[h_channel+"vbf_eff_study_dPhijj_met130"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                              hMap1D[h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70"+sysName]->Fill(mjj * 0.001, 1.);
                              hMap1D[h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            }
                            if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                              hMap1D[h_channel+"vbf_eff_study_mjj_met130_pass_HLT_xe70_tclcw"+sysName]->Fill(mjj * 0.001, 1.);
                              hMap1D[h_channel+"vbf_eff_study_dPhijj_met130_pass_HLT_xe70_tclcw"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            }
                          }
                          // MET > 150 GeV
                          if ( emulMET_Wmunu > 150000. ) {
                            hMap1D[h_channel+"vbf_eff_study_mjj_met150"+sysName]->Fill(mjj * 0.001, 1.);
                            hMap1D[h_channel+"vbf_eff_study_dPhijj_met150"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                              hMap1D[h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70"+sysName]->Fill(mjj * 0.001, 1.);
                              hMap1D[h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            }
                            if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                              hMap1D[h_channel+"vbf_eff_study_mjj_met150_pass_HLT_xe70_tclcw"+sysName]->Fill(mjj * 0.001, 1.);
                              hMap1D[h_channel+"vbf_eff_study_dPhijj_met150_pass_HLT_xe70_tclcw"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            }
                          }
                          // MET > 200 GeV
                          if ( emulMET_Wmunu > 200000. ) {
                            hMap1D[h_channel+"vbf_eff_study_mjj_met200"+sysName]->Fill(mjj * 0.001, 1.);
                            hMap1D[h_channel+"vbf_eff_study_dPhijj_met200"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                              hMap1D[h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70"+sysName]->Fill(mjj * 0.001, 1.);
                              hMap1D[h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
                            }
                            if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") ) {
                              hMap1D[h_channel+"vbf_eff_study_mjj_met200_pass_HLT_xe70_tclcw"+sysName]->Fill(mjj * 0.001, 1.);
                              hMap1D[h_channel+"vbf_eff_study_dPhijj_met200_pass_HLT_xe70_tclcw"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), 1.);
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
      }
    }





    //-----------------------------------------------------
    // Z -> mumu + JET Multijet Background study (Method 1)
    //-----------------------------------------------------

    if (m_isZmumu) {
      h_channel = "h_zmumu_";

      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
        if ( m_goodMuonForZ->size() > 1 && m_goodElectron->size() == 0 && m_goodTau->size() == 0 ) { // Letopn veto
          if (numExtra == 0 && pass_dimuonPtCut) {
            if ( emulMET_Zmumu > m_metCut ) {
              if ( m_goodJet->size() > 0 ) {

                ////////////////////////
                // MonoJet phasespace //
                ////////////////////////
                if ( pass_monoJet && pass_dPhijetmet_Zmumu ) {

                  // Calculate muon SF for Zmumu
                  float mcEventWeight_Zmumu = 1.;
                  if (!m_isData) {
                    double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                    //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                    mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                  }

                  // Fill histogram
                  // All charge muon
                  hMap1D[h_channel+"monojet_multijet_study_mll_all_lep"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                  if (mll_muon > m_mllMin && mll_muon < m_mllMax){
                    hMap1D[h_channel+"monojet_multijet_study_met_emulmet_all_lep"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                  }
                  // Opposite sign charge muon
                  if ( pass_OSmuon ) {
                    hMap1D[h_channel+"monojet_multijet_study_mll_os_lep"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                    if (mll_muon > m_mllMin && mll_muon < m_mllMax){
                      hMap1D[h_channel+"monojet_multijet_study_met_emulmet_os_lep"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                    }
                  }
                  // Same sign charge muon
                  if ( pass_SSmuon ) {
                    hMap1D[h_channel+"monojet_multijet_study_mll_ss_lep"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                    if (mll_muon > m_mllMin && mll_muon < m_mllMax){
                      hMap1D[h_channel+"monojet_multijet_study_met_emulmet_ss_lep"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                    }
                  }
                } // monojet cut

                ////////////////////
                // VBF phasespace //
                ////////////////////
                if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zmumu) {

                  // Calculate muon SF for Zmumu
                  float mcEventWeight_Zmumu = 1.;
                  if (!m_isData) {
                    double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                    //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                    mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                  }

                  // Fill histogram
                  // All charge muon
                  hMap1D[h_channel+"vbf_multijet_study_mll_all_lep"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                  if (mll_muon > m_mllMin && mll_muon < m_mllMax){
                    hMap1D[h_channel+"vbf_multijet_study_met_emulmet_all_lep"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                    hMap1D[h_channel+"vbf_multijet_study_mjj_all_lep"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zmumu);
                    hMap1D[h_channel+"vbf_multijet_study_dPhijj_all_lep"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                  }
                  // Opposite sign charge muon
                  if ( pass_OSmuon ) {
                    hMap1D[h_channel+"vbf_multijet_study_mll_os_lep"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                    if (mll_muon > m_mllMin && mll_muon < m_mllMax){
                      hMap1D[h_channel+"vbf_multijet_study_met_emulmet_os_lep"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                      hMap1D[h_channel+"vbf_multijet_study_mjj_os_lep"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zmumu);
                      hMap1D[h_channel+"vbf_multijet_study_dPhijj_os_lep"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                    }
                  }
                  // Same sign charge muon
                  if ( pass_SSmuon ) {
                    hMap1D[h_channel+"vbf_multijet_study_mll_ss_lep"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
                    if (mll_muon > m_mllMin && mll_muon < m_mllMax){
                      hMap1D[h_channel+"vbf_multijet_study_met_emulmet_ss_lep"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                      hMap1D[h_channel+"vbf_multijet_study_mjj_ss_lep"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zmumu);
                      hMap1D[h_channel+"vbf_multijet_study_dPhijj_ss_lep"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                    }
                  }
                } // VBF cut


              } // At least 1 jet
            } // MET cut
          } // Dimuon cut
        } // Lepton veto
      } // MET Trigger

    } // end multijet study



    //---------------------------------------------------
    // Z -> ee + JET Multijet Background study (Method 1)
    //---------------------------------------------------

    if (m_isZee) {
      h_channel = "h_zee_";

      if ( m_goodJet->size() > 0 ) {
        if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){
          if ( m_goodMuon->size() == 0 && m_goodTau->size() == 0 ) { // Letopn veto

            //-------------------------------------------------
            // OS and SS charge for multijet background study
            //-------------------------------------------------

            if (m_goodElectron->size() == 2 && pass_dielectronPtCut) {
              if ( emulMET_Zee > m_metCut ) {

                ////////////////////////
                // MonoJet phasespace //
                ////////////////////////
                if ( pass_monoJet && pass_dPhijetmet_Zee ) {

                  // Calculate electron SF
                  float mcEventWeight_Zee = 1.;
                  if (!m_isData) {
                    float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                    //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                    mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                  }

                  // Fill histogram
                  // All charge electron
                  hMap1D[h_channel+"monojet_multijet_study_mll_all_lep"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                  if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                    hMap1D[h_channel+"monojet_multijet_study_met_emulmet_all_lep"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                  }
                  // Opposite sign charge electron
                  if ( pass_OSelectron ) {
                    hMap1D[h_channel+"monojet_multijet_study_mll_os_lep"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                    if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                      hMap1D[h_channel+"monojet_multijet_study_met_emulmet_os_lep"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                    }
                  }
                  // Same sign charge electron
                  if ( pass_SSelectron ) {
                    hMap1D[h_channel+"monojet_multijet_study_mll_ss_lep"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                    if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                      hMap1D[h_channel+"monojet_multijet_study_met_emulmet_ss_lep"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                    }
                  }
                } // monojet cut

                ////////////////////
                // VBF phasespace //
                ////////////////////
                if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zee) {

                  // Calculate electron SF
                  float mcEventWeight_Zee = 1.;
                  if (!m_isData) {
                    float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                    //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                    mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                  }

                  // Fill histogram
                  // All charge electron
                  hMap1D[h_channel+"vbf_multijet_study_mll_all_lep"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                  if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                    hMap1D[h_channel+"vbf_multijet_study_met_emulmet_all_lep"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                    hMap1D[h_channel+"vbf_multijet_study_mjj_all_lep"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zee);
                    hMap1D[h_channel+"vbf_multijet_study_dPhijj_all_lep"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                  }
                  // Opposite sign charge electron
                  if ( pass_OSelectron ) {
                    hMap1D[h_channel+"vbf_multijet_study_mll_os_lep"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                    if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                      hMap1D[h_channel+"vbf_multijet_study_met_emulmet_os_lep"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                      hMap1D[h_channel+"vbf_multijet_study_mjj_os_lep"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zee);
                      hMap1D[h_channel+"vbf_multijet_study_dPhijj_os_lep"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                    }
                  }
                  // Same sign charge electron
                  if ( pass_SSelectron ) {
                    hMap1D[h_channel+"vbf_multijet_study_mll_ss_lep"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
                    if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                      hMap1D[h_channel+"vbf_multijet_study_met_emulmet_ss_lep"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                      hMap1D[h_channel+"vbf_multijet_study_mjj_ss_lep"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zee);
                      hMap1D[h_channel+"vbf_multijet_study_dPhijj_ss_lep"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                    }
                  }


                } // VBF cut
              } // MET cut
            } // Dielectron cut
          } // Lepton veto
        } // Electron Trigger
      } // At least 1 jet

    } // end multijet study



    //-----------------------------------------------------
    // Z -> mumu + JET Multijet Background study (Method 2)
    //-----------------------------------------------------

    ///////////////////////////////////////
    // Nominal Muon cut with MC and Data //
    ///////////////////////////////////////

    if (m_isZmumu && sysName == "") {
      h_channel = "h_zmumu_";

      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {

        if ( m_goodJet->size() > 0 && m_goodElectron->size() == 0 && m_goodTau->size() == 0 && m_goodMuonForZ->size() > 1 ) {

          if ( pass_dimuonPtCut && pass_OSmuon && numExtra == 0 ) {

            ////////////////////////
            // MonoJet phasespace //
            ////////////////////////
            if ( pass_monoJet && pass_dPhijetmet_Zmumu ) {

              // Calculate muon SF for Zmumu
              float mcEventWeight_Zmumu = 1.;
              if (!m_isData) {
                double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
              }

              // Fill histogram
              // 200 < MET < 300 GeV
              if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                hMap1D[h_channel+"monojet_qcd_method2_nominal_cut_met200_300_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
              }
              // 300 < MET < 500 GeV
              if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                hMap1D[h_channel+"monojet_qcd_method2_nominal_cut_met300_500_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
              }
              // MET > 500 GeV
              if ( emulMET_Zmumu > 500000.  ) {
                hMap1D[h_channel+"monojet_qcd_method2_nominal_cut_met500_inf_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
              }

            } // Monojet


            ////////////////////
            // VBF phasespace //
            ////////////////////
            if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zmumu) {

              // Calculate muon SF for Zmumu
              float mcEventWeight_Zmumu = 1.;
              if (!m_isData) {
                double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
              }

              // Fill histogram
              // 200 < MET < 300 GeV
              if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                hMap1D[h_channel+"vbf_qcd_method2_nominal_cut_met200_300_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
              }
              // 300 < MET < 500 GeV
              if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                hMap1D[h_channel+"vbf_qcd_method2_nominal_cut_met300_500_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
              }
              // MET > 500 GeV
              if ( emulMET_Zmumu > 500000.  ) {
                hMap1D[h_channel+"vbf_qcd_method2_nominal_cut_met500_inf_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);
              }


            } // VBF

          } // Norminal cut

        } // At least 1 jet && electron and tau veto
      } // MET Trigger
    } // m_isZmumu with MC and Data


    ////////////////////////////////
    // Reverse Muon cut with Data //
    ////////////////////////////////

    if (m_isZmumu && m_isData) {
      h_channel = "h_zmumu_";

      if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {

        if ( m_goodJet->size() > 0 && m_goodElectron->size() == 0 && m_goodTau->size() == 0 ) {

          if ( pass_dimuonPtCut && m_baselineMuon->size() == 2 ) {

            // S-S or O-S charge decision
            float muon_OS = false;
            if ( m_baselineMuon->at(0)->charge() * m_baselineMuon->at(1)->charge() < 0 ) muon_OS = true;

            // Mll calculation
            TLorentzVector baselineMuon1 = m_baselineMuon->at(0)->p4();
            TLorentzVector baselineMuon2 = m_baselineMuon->at(1)->p4();
            auto Zmass_baselineMuon = baselineMuon1 + baselineMuon2;
            float mll_baselineMuon = Zmass_baselineMuon.M();

            // Reverse cut decisions
            bool m_case1_cut = false;
            bool m_case2_cut = false;
            bool m_case3_cut = false;
            // Case 1: Fail ID, Fail Iso, Pass OS
            if ( !m_loosemuonSelection->accept(*m_baselineMuon->at(0)) || !m_loosemuonSelection->accept(*m_baselineMuon->at(1)) ) { // Fail ID
              if ( !m_IsoToolVBF->accept(*m_baselineMuon->at(0)) || !m_IsoToolVBF->accept(*m_baselineMuon->at(1)) ) { // Fail Iso
                if ( muon_OS ) { // Pass OS
                  m_case1_cut = true;
                }
              }
            }
            // Case 2: Fail ID, Pass Iso, Fail OS
            if ( !m_loosemuonSelection->accept(*m_baselineMuon->at(0)) || !m_loosemuonSelection->accept(*m_baselineMuon->at(1)) ) { // Fail ID
              if ( m_IsoToolVBF->accept(*m_baselineMuon->at(0)) && m_IsoToolVBF->accept(*m_baselineMuon->at(1)) ) { // Pass Iso
                if ( !muon_OS ) { // Fail OS
                  m_case2_cut = true;
                }
              }
            }
            // Case 3: Pass ID, Fail Iso, Fail OS
            if ( m_loosemuonSelection->accept(*m_baselineMuon->at(0)) && m_loosemuonSelection->accept(*m_baselineMuon->at(1)) ) { // Pass ID
              if ( !m_IsoToolVBF->accept(*m_baselineMuon->at(0)) || !m_IsoToolVBF->accept(*m_baselineMuon->at(1)) ) { // Fail Iso
                if ( !muon_OS ) { // Fail OS
                  m_case3_cut = true;
                }
              }
            }


            ////////////////////////
            // MonoJet phasespace //
            ////////////////////////
            if ( pass_monoJet && pass_dPhijetmet_Zmumu ) {

              // Calculate muon SF for Zmumu
              float mcEventWeight_Zmumu = 1.;
              if (!m_isData) {
                double totalMuonSF_Zmumu = GetTotalMuonSF(*m_baselineMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
              }

              // Fill histogram
              // Case 1
              if ( m_case1_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case1_cut_met200_300_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case1_cut_met300_500_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // MET > 500 GeV
                if ( emulMET_Zmumu > 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case1_cut_met500_inf_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
              }
              // Case 2
              if ( m_case2_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case2_cut_met200_300_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case2_cut_met300_500_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // MET > 500 GeV
                if ( emulMET_Zmumu > 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case2_cut_met500_inf_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
              }
              // Case 3
              if ( m_case3_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case3_cut_met200_300_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case3_cut_met300_500_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // MET > 500 GeV
                if ( emulMET_Zmumu > 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case3_cut_met500_inf_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
              }

            } // monojet


            ////////////////////
            // VBF phasespace //
            ////////////////////
            if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zmumu) {

              // Calculate muon SF for Zmumu
              float mcEventWeight_Zmumu = 1.;
              if (!m_isData) {
                double totalMuonSF_Zmumu = GetTotalMuonSF(*m_baselineMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
              }

              // Fill histogram
              // Case 1
              if ( m_case1_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case1_cut_met200_300_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case1_cut_met300_500_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // MET > 500 GeV
                if ( emulMET_Zmumu > 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case1_cut_met500_inf_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
              }
              // Case 2
              if ( m_case2_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case2_cut_met200_300_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case2_cut_met300_500_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // MET > 500 GeV
                if ( emulMET_Zmumu > 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case2_cut_met500_inf_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
              }
              // Case 3
              if ( m_case3_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zmumu > 200000. && emulMET_Zmumu < 300000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case3_cut_met200_300_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zmumu > 300000. && emulMET_Zmumu < 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case3_cut_met300_500_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
                // MET > 500 GeV
                if ( emulMET_Zmumu > 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case3_cut_met500_inf_mll"+sysName]->Fill(mll_baselineMuon * 0.001, mcEventWeight_Zmumu);
                }
              }

            } // VBF

          } // Reverse cut
        } // At least 1 jet && electron and tau veto
      } // MET trigger
    } // m_isZmumu with Data






    //---------------------------------------------------
    // Z -> ee + JET Multijet Background study (Method 2)
    //---------------------------------------------------

    ///////////////////////////////////////////
    // Nominal Electron cut with MC and Data //
    ///////////////////////////////////////////

    if (m_isZee && sysName == "") {
      h_channel = "h_zee_";

      if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){

        if ( m_goodJet->size() > 0 && m_goodMuon->size() == 0 && m_goodTau->size() == 0 ) {

          if ( pass_dielectronPtCut && pass_OSelectron && m_goodElectron->size() == 2 ) {

            ////////////////////////
            // MonoJet phasespace //
            ////////////////////////
            if ( pass_monoJet && pass_dPhijetmet_Zee ) {

              // Calculate electron SF
              float mcEventWeight_Zee = 1.;
              if (!m_isData) {
                float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
              }

              // Fill histogram
              // 200 < MET < 300 GeV
              if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                hMap1D[h_channel+"monojet_qcd_method2_nominal_cut_met200_300_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
              }
              // 300 < MET < 500 GeV
              if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                hMap1D[h_channel+"monojet_qcd_method2_nominal_cut_met300_500_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
              }
              // MET > 500 GeV
              if ( emulMET_Zee > 500000.  ) {
                hMap1D[h_channel+"monojet_qcd_method2_nominal_cut_met500_inf_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
              }

            } // Monojet


            ////////////////////
            // VBF phasespace //
            ////////////////////
            if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zee) {

              // Calculate electron SF
              float mcEventWeight_Zee = 1.;
              if (!m_isData) {
                float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
              }

              // Fill histogram
              // 200 < MET < 300 GeV
              if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                hMap1D[h_channel+"vbf_qcd_method2_nominal_cut_met200_300_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
              }
              // 300 < MET < 500 GeV
              if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                hMap1D[h_channel+"vbf_qcd_method2_nominal_cut_met300_500_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
              }
              // MET > 500 GeV
              if ( emulMET_Zee > 500000.  ) {
                hMap1D[h_channel+"vbf_qcd_method2_nominal_cut_met500_inf_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);
              }


            } // VBF

          } // Norminal cut

        } // At least 1 jet && muon and tau veto
      } // Electron trigger
    } // m_isZee with MC and Data


    ////////////////////////////////////
    // Reverse Electron cut with Data //
    ////////////////////////////////////

    if (m_isZee && m_isData) {
      h_channel = "h_zee_";

      if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){

        if ( m_goodJet->size() > 0 && m_goodMuon->size() == 0 && m_goodTau->size() == 0 ) {

          if ( pass_dielectronPtCut && m_baselineElectron->size() == 2 ) {

            // S-S or O-S charge decision
            float elec_OS = false;
            if ( m_baselineElectron->at(0)->charge() * m_baselineElectron->at(1)->charge() < 0 ) elec_OS = true;

            // Mll calculation
            TLorentzVector baselineElectron1 = m_baselineElectron->at(0)->p4();
            TLorentzVector baselineElectron2 = m_baselineElectron->at(1)->p4();
            auto Zmass_baselineElectron = baselineElectron1 + baselineElectron2;
            float mll_baselineElectron = Zmass_baselineElectron.M();

            // Reverse cut decisions
            bool m_case1_cut = false;
            bool m_case2_cut = false;
            bool m_case3_cut = false;
            // Case 1: Fail ID, Fail Iso, Pass OS
            if ( !m_LHToolLoose2015->accept(*m_baselineElectron->at(0)) || !m_LHToolLoose2015->accept(*m_baselineElectron->at(1)) ) { // Fail ID
              if ( !m_IsoToolVBF->accept(*m_baselineElectron->at(0)) || !m_IsoToolVBF->accept(*m_baselineElectron->at(1)) ) { // Fail Iso
                if ( elec_OS ) { // Pass OS
                  m_case1_cut = true;
                }
              }
            }
            // Case 2: Fail ID, Pass Iso, Fail OS
            if ( !m_LHToolLoose2015->accept(*m_baselineElectron->at(0)) || !m_LHToolLoose2015->accept(*m_baselineElectron->at(1)) ) { // Fail ID
              if ( m_IsoToolVBF->accept(*m_baselineElectron->at(0)) && m_IsoToolVBF->accept(*m_baselineElectron->at(1)) ) { // Pass Iso
                if ( !elec_OS ) { // Fail OS
                  m_case2_cut = true;
                }
              }
            }
            // Case 3: Pass ID, Fail Iso, Fail OS
            if ( m_LHToolLoose2015->accept(*m_baselineElectron->at(0)) && m_LHToolLoose2015->accept(*m_baselineElectron->at(1)) ) { // Pass ID
              if ( !m_IsoToolVBF->accept(*m_baselineElectron->at(0)) || !m_IsoToolVBF->accept(*m_baselineElectron->at(1)) ) { // Fail Iso
                if ( !elec_OS ) { // Fail OS
                  m_case3_cut = true;
                }
              }
            }


            ////////////////////////
            // MonoJet phasespace //
            ////////////////////////
            if ( pass_monoJet && pass_dPhijetmet_Zee ) {

              // Calculate electron SF
              float mcEventWeight_Zee = 1.;
              if (!m_isData) {
                float totalElectronSF_Zee = GetTotalElectronSF(*m_baselineElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
              }

              // Fill histogram
              // Case 1
              if ( m_case1_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case1_cut_met200_300_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case1_cut_met300_500_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // MET > 500 GeV
                if ( emulMET_Zee > 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case1_cut_met500_inf_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
              }
              // Case 2
              if ( m_case2_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case2_cut_met200_300_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case2_cut_met300_500_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // MET > 500 GeV
                if ( emulMET_Zee > 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case2_cut_met500_inf_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
              }
              // Case 3
              if ( m_case3_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case3_cut_met200_300_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case3_cut_met300_500_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // MET > 500 GeV
                if ( emulMET_Zee > 500000.  ) {
                  hMap1D[h_channel+"monojet_qcd_method2_case3_cut_met500_inf_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
              }

            } // monojet


            ////////////////////
            // VBF phasespace //
            ////////////////////
            if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zee) {

              // Calculate electron SF
              float mcEventWeight_Zee = 1.;
              if (!m_isData) {
                float totalElectronSF_Zee = GetTotalElectronSF(*m_baselineElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
              }

              // Fill histogram
              // Case 1
              if ( m_case1_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case1_cut_met200_300_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case1_cut_met300_500_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // MET > 500 GeV
                if ( emulMET_Zee > 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case1_cut_met500_inf_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
              }
              // Case 2
              if ( m_case2_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case2_cut_met200_300_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case2_cut_met300_500_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // MET > 500 GeV
                if ( emulMET_Zee > 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case2_cut_met500_inf_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
              }
              // Case 3
              if ( m_case3_cut ) { 
                // 200 < MET < 300 GeV
                if ( emulMET_Zee > 200000. && emulMET_Zee < 300000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case3_cut_met200_300_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // 300 < MET < 500 GeV
                if ( emulMET_Zee > 300000. && emulMET_Zee < 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case3_cut_met300_500_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
                // MET > 500 GeV
                if ( emulMET_Zee > 500000.  ) {
                  hMap1D[h_channel+"vbf_qcd_method2_case3_cut_met500_inf_mll"+sysName]->Fill(mll_baselineElectron * 0.001, mcEventWeight_Zee);
                }
              }

            } // VBF

          } // Reverse cut

        } // At least 1 jet && muon and tau veto
      } // Electron trigger
    } // m_isZee with Data



    //------------------------
    // Z -> nunu + JET in SM1
    //------------------------


    if (m_isZnunu && sysName == "") {
      h_channel = "h_znunu_";

      if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") && MET > 150000. ) {

        if ( m_goodJet->size() > 0 && m_goodElectron->size() == 0 && m_goodMuon->size() == 0 && m_goodTau->size() == 0 ) {

          if ( pass_sm1Jet && pass_dPhijetmet ) {

            // Fill histogram
            // For publication
            hMap1D[h_channel+"sm1_met"+sysName]->Fill(MET * 0.001, mcEventWeight);
            // Jets
            hMap1D[h_channel+"sm1_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight);
            hMap1D[h_channel+"sm1_jet_pt"+sysName]->Fill(sm1jet_pt * 0.001, mcEventWeight);
            hMap1D[h_channel+"sm1_jet_phi"+sysName]->Fill(sm1jet_phi, mcEventWeight);
            hMap1D[h_channel+"sm1_jet_eta"+sysName]->Fill(sm1jet_eta, mcEventWeight);
            hMap1D[h_channel+"sm1_jet_rap"+sysName]->Fill(sm1jet_rapidity, mcEventWeight);
            hMap1D[h_channel+"sm1_dPhimetjet"+sysName]->Fill(dPhiSM1jetMet, mcEventWeight);
            hMap1D[h_channel+"sm1_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet, mcEventWeight);

          } // sm1jet
        } // Veto
      } // MET trigger and MET cut
    } // m_isZnunu



    //------------------------
    // Z -> mumu + JET in SM1
    //------------------------


    if (m_isZmumu && sysName == "") {
      h_channel = "h_zmumu_";

      if ( m_trigDecisionTool->isPassed("HLT_xe70_tc_lcw") && emulMET_Zmumu > 150000. ) {

        if ( m_goodJet->size() > 0 && m_goodElectron->size() == 0 && m_goodTau->size() == 0 && m_goodMuonForZ->size() > 1 ) {

          if ( pass_dimuonPtCut && pass_OSmuon && numExtra == 0 && mll_muon > m_mllMin && mll_muon < m_mllMax ) {

            if ( pass_sm1Jet && pass_dPhijetmet_Zmumu ) {

              // Calculate muon SF for Zmumu
              float mcEventWeight_Zmumu = 1.;
              if (!m_isData) {
                double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
              }

              // Fill histogram
              // For publication
              hMap1D[h_channel+"sm1_met_emulmet"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
              // Jets
              hMap1D[h_channel+"sm1_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_jet_pt"+sysName]->Fill(sm1jet_pt * 0.001, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_jet_phi"+sysName]->Fill(sm1jet_phi, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_jet_eta"+sysName]->Fill(sm1jet_eta, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_jet_rap"+sysName]->Fill(sm1jet_rapidity, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_dPhimetjet"+sysName]->Fill(dPhiSM1jetMet_Zmumu, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_Zmumu, mcEventWeight_Zmumu);
              // Leptons
              hMap1D[h_channel+"sm1_lepton1_pt"+sysName]->Fill(m_goodMuonForZ->at(0)->pt() * 0.001, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_lepton2_pt"+sysName]->Fill(m_goodMuonForZ->at(1)->pt() * 0.001, mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_lepton1_phi"+sysName]->Fill(m_goodMuonForZ->at(0)->phi(), mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_lepton2_phi"+sysName]->Fill(m_goodMuonForZ->at(1)->phi(), mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_lepton1_eta"+sysName]->Fill(m_goodMuonForZ->at(0)->eta(), mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_lepton2_eta"+sysName]->Fill(m_goodMuonForZ->at(1)->eta(), mcEventWeight_Zmumu);
              hMap1D[h_channel+"sm1_mll"+sysName]->Fill(mll_muon * 0.001, mcEventWeight_Zmumu);

            } // sm1jet
          } // dimuon
        } // Veto
      } // MET trigger and MET cut
    } // m_isZmumu




    //----------------------
    // Z -> ee + JET in SM1
    //----------------------

    if (m_isZee && sysName == "") {
      h_channel = "h_zee_";

      if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")) {

        if ( emulMET_Zee > 150000. ) {

          if ( m_goodJet->size() > 0 && m_goodMuon->size() == 0 && m_goodTau->size() == 0 ) {

            if ( pass_dielectronPtCut && pass_OSelectron && m_goodElectron->size() == 2 && mll_electron > m_mllMin && mll_electron < m_mllMax ) {

              if ( pass_sm1Jet && pass_dPhijetmet_Zee ) {

                // Calculate electron SF
                float mcEventWeight_Zee = 1.;
                if (!m_isData) {
                  float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                  //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                  mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                }

                // Fill histogram
                // For publication
                hMap1D[h_channel+"sm1_met_emulmet"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                // Jets
                hMap1D[h_channel+"sm1_njet"+sysName]->Fill(m_goodJet->size(), mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_jet_pt"+sysName]->Fill(sm1jet_pt * 0.001, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_jet_phi"+sysName]->Fill(sm1jet_phi, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_jet_eta"+sysName]->Fill(sm1jet_eta, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_jet_rap"+sysName]->Fill(sm1jet_rapidity, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_dPhimetjet"+sysName]->Fill(dPhiSM1jetMet_Zee, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_dPhiMinmetjet"+sysName]->Fill(dPhiMinjetmet_Zee, mcEventWeight_Zee);
                // Leptons
                hMap1D[h_channel+"sm1_lepton1_pt"+sysName]->Fill(m_goodElectron->at(0)->pt() * 0.001, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_lepton2_pt"+sysName]->Fill(m_goodElectron->at(1)->pt() * 0.001, mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_lepton1_phi"+sysName]->Fill(m_goodElectron->at(0)->phi(), mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_lepton2_phi"+sysName]->Fill(m_goodElectron->at(1)->phi(), mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_lepton1_eta"+sysName]->Fill(m_goodElectron->at(0)->eta(), mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_lepton2_eta"+sysName]->Fill(m_goodElectron->at(1)->eta(), mcEventWeight_Zee);
                hMap1D[h_channel+"sm1_mll"+sysName]->Fill(mll_electron * 0.001, mcEventWeight_Zee);

              } // sm1jet
            } // dilepton
          } // Veto
        } // MET cut
      } // Electron trigger
    } // m_isZee








    //-------------------------------
    // Z -> mumu + JET EVENT Cutflow
    //-------------------------------

    if (m_isZmumu && m_isEmilyCutflow && sysName == ""){

      if ( (m_goodJet->size() > 0 && monojet_pt > 100000.) || (m_goodJet->size() > 1 && jet1_pt > 55000. && jet2_pt > 45000.) ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Skim cuts");
        if (m_goodMuonForZ->size() > 1) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]At least Two Muon");
          if (pass_OSmuon) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Opposite sign charge");
            if (pass_dimuonPtCut) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Dimuon pT cut");
              //if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]MET Trigger");
                if (mll_muon > m_mllMin && mll_muon < m_mllMax) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Zmass window");
                  /*
                  // MET test
                  if (emulMET_Zmumu < m_metCut) {
                    Info("execute()", "  emulMET_Zmumu = %.2f GeV", emulMET_Zmumu * 0.001);
                    int jetCount = 0;
                    for (const auto& jet : *m_goodJet) {
                      jetCount++;
                      Info("execute()", " jet # : %i", jetCount);
                      Info("execute()", " jet pt = %.3f GeV", jet->pt() * 0.001);
                      Info("execute()", " jet eta = %.3f GeV", jet->eta());
                      Info("execute()", " jet phi = %.3f GeV", jet->phi());
                    }
                  }
                  */
                  if (emulMET_Zmumu > m_metCut) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]MET cut");
                    /*
                       if (m_goodTau->size() > 0){
                       Info("execute()", "=====================================");
                       Info("execute()", " Event # = %llu", eventInfo->eventNumber());
                       int tauCount = 0;
                       for (const auto& tau : *m_goodTau) {
                       tauCount++;
                       Info("execute()", " tau # : %i", tauCount);
                       Info("execute()", " tau pt = %.3f GeV", tau->pt() * 0.001);
                    //Info("execute()", " tau charge = %.1f", tau->charge());
                    //Info("execute()", " tau eta = %.3f", tau->eta());
                    //Info("execute()", " tau phi = %.3f", tau->phi());
                    }
                    }
                    */
                    if (numExtra == 0) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Exact two muon");
                      if (m_goodElectron->size() == 0) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Electron veto");
                        if (m_goodTau->size() == 0) {
                          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Tau veto");
                          ////////////////////////
                          // MonoJet phasespace //
                          ////////////////////////
                          if (pass_monoJet && pass_dPhijetmet_Zmumu) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]Monojet cut");
/*
                            Info("execute()", "================================================");
                            Info("execute()", " Event # = %llu in monojet phasespace", eventInfo->eventNumber());
                            Info("execute()", "------------------------------------------------");
                            Info("execute()", " Original mcWeight = %f", print_mcWeight_origin);
                            Info("execute()", " No high weight  mcWeight = %f", print_mcWeight_nohigh);
                            Info("execute()", " PU Weight = %f", print_puweight);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight) = %f", print_mcWeight_nohigh*print_puweight);
                            Info("execute()", " sherpaReweight = %lf", print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight * sherpa_weight) = %f", print_mcWeight_nohigh*print_puweight*print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight = %f", mcEventWeight);
*/
                            // Calculate muon SF for Zmumu
                            float mcEventWeight_Zmumu = 1.;
                            if (!m_isData) {
                              double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                              //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                              mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                            }

                            hMap1D["Emily_Zmumu_MET_mono"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                          }
                          ////////////////////
                          // VBF phasespace //
                          ////////////////////
                          if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zmumu) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zmumu]VBF cut");
/*
                            Info("execute()", "================================================");
                            Info("execute()", " Event # = %llu in VBF phasespace", eventInfo->eventNumber());
                            Info("execute()", "------------------------------------------------");
                            Info("execute()", " Original mcWeight = %f", print_mcWeight_origin);
                            Info("execute()", " No high weight  mcWeight = %f", print_mcWeight_nohigh);
                            Info("execute()", " PU Weight = %f", print_puweight);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight) = %f", print_mcWeight_nohigh*print_puweight);
                            Info("execute()", " sherpaReweight = %lf", print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight * sherpa_weight) = %f", print_mcWeight_nohigh*print_puweight*print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight = %f", mcEventWeight);
*/
                            // Calculate muon SF for Zmumu
                            float mcEventWeight_Zmumu = 1.;
                            if (!m_isData) {
                              double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoMuonSFforZ, m_ttvaSF);
                              //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                              mcEventWeight_Zmumu = mcEventWeight * totalMuonSF_Zmumu;
                            }

                            hMap1D["Emily_Zmumu_MET_search"+sysName]->Fill(emulMET_Zmumu * 0.001, mcEventWeight_Zmumu);
                            hMap1D["Emily_Zmumu_Mjj_search"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zmumu);
                            hMap1D["Emily_Zmumu_DeltaPhiAll"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zmumu);
                          }
                        } // Tau veto
                      } // Electron veto
                    } // Exact two muon
                  } // MET cut
                } // Zmass window
              //} // MET trigger
            } // Dimuon pT cut
          } // Opposite sign change cut
        } // At least 2 muons
      } // Skim cuts

    } // end cutflow





    //-------------------------------
    // Z -> ee + JET EVENT Cutflow
    //-------------------------------

    if (m_isZee && m_isEmilyCutflow && sysName == ""){

      if ( (m_goodJet->size() > 0 && monojet_pt > 100000.) || (m_goodJet->size() > 1 && jet1_pt > 55000. && jet2_pt > 45000.) ) {
        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Skim cuts");
        if (m_goodElectron->size() > 1) {
          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]At least Two Electron");
             /*
             Info("execute()", "=====================================");
             Info("execute()", " Event # = %llu", eventInfo->eventNumber());
             int eleCount = 0;
             for (const auto& electron : *m_goodElectron) {
             eleCount++;
             Info("execute()", " electron # : %i", eleCount);
             Info("execute()", " electron pt = %.3f GeV", electron->pt() * 0.001);
             Info("execute()", " electron eta = %.3f", electron->eta());
             Info("execute()", " electron phi = %.3f", electron->phi());
             }
             */
          if (pass_OSelectron) {
            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Opposite sign charge");
            if (pass_dielectronPtCut) {
              if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Dielectron pT cut");
              if ((!m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM18VH")) || (m_isData && m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH")) || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose")){
                if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Electron Trigger");
                if (mll_electron > m_mllMin && mll_electron < m_mllMax) {
                  if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Zmass window");
                  //Info("execute()", "  # Electron = %llu, # Muon = %llu, # Tau = %llu", m_goodElectron->size(), m_goodMuon->size(), m_goodTau->size());
                  if (emulMET_Zee > m_metCut) {
                    if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]MET cut");
                    if (m_goodMuon->size() == 0) {
                      if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Muon veto");
                      if (m_goodElectron->size() == 2) {
                        if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Exact two electrons");
                        if (m_goodTau->size() == 0) {
                          if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Tau veto");
                          ////////////////////////
                          // MonoJet phasespace //
                          ////////////////////////
                          if (pass_monoJet && pass_dPhijetmet_Zee) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]Monojet cut");
/*
                            Info("execute()", "================================================");
                            Info("execute()", " Event # = %llu in monojet phasespace", eventInfo->eventNumber());
                            Info("execute()", "------------------------------------------------");
                            Info("execute()", " Original mcWeight = %f", print_mcWeight_origin);
                            Info("execute()", " No high weight  mcWeight = %f", print_mcWeight_nohigh);
                            Info("execute()", " PU Weight = %f", print_puweight);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight) = %f", print_mcWeight_nohigh*print_puweight);
                            Info("execute()", " sherpaReweight = %lf", print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight * sherpa_weight) = %f", print_mcWeight_nohigh*print_puweight*print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight = %f", mcEventWeight);
*/
                            // Calculate electron SF
                            float mcEventWeight_Zee = 1.;
                            if (!m_isData) {
                              float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                              //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                              mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                            }

                            hMap1D["Emily_Zee_MET_mono"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                          }
                          ////////////////////
                          // VBF phasespace //
                          ////////////////////
                          if (pass_diJet && mjj > m_mjjCut && pass_CJV && pass_dPhijetmet_Zee) {
                            if (sysName == "" && m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Emily, Zee]VBF cut");
/*
                            Info("execute()", "================================================");
                            Info("execute()", " Event # = %llu in VBF phasespace", eventInfo->eventNumber());
                            Info("execute()", "------------------------------------------------");
                            Info("execute()", " Original mcWeight = %f", print_mcWeight_origin);
                            Info("execute()", " No high weight  mcWeight = %f", print_mcWeight_nohigh);
                            Info("execute()", " PU Weight = %f", print_puweight);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight) = %f", print_mcWeight_nohigh*print_puweight);
                            Info("execute()", " sherpaReweight = %lf", print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight (mcWeight * pu_weight * sherpa_weight) = %f", print_mcWeight_nohigh*print_puweight*print_sherpaReweightValue);
                            Info("execute()", " mcEventWeight = %f", mcEventWeight);
*/
                            // Calculate electron SF
                            float mcEventWeight_Zee = 1.;
                            if (!m_isData) {
                              float totalElectronSF_Zee = GetTotalElectronSF(*m_goodElectron, m_recoSF, m_idSF, m_isoElectronSF, m_trigSF);
                              //Info("execute()", " Zee Total Electron SF = %.3f ", totalElectronSF_Zee);
                              mcEventWeight_Zee = mcEventWeight * totalElectronSF_Zee;
                            }

                            hMap1D["Emily_Zee_MET_search"+sysName]->Fill(emulMET_Zee * 0.001, mcEventWeight_Zee);
                            hMap1D["Emily_Zee_Mjj_search"+sysName]->Fill(mjj * 0.001, mcEventWeight_Zee);
                            hMap1D["Emily_Zee_DeltaPhiAll"+sysName]->Fill(deltaPhi(jet1_phi, jet2_phi), mcEventWeight_Zee);
                          }
                        } // Tau veto
                      } // Exact two electrons
                    } // Muon veto
                  } // MET cut
                } // Zmass window
              } // Electron trigger
            } // Dielectron pT cut
          } // Opposite sign change cut
        } // At least 2 muons
      } // Skim cuts

    } // end Cutflow



    //////////////////////////////////
    // Delete shallow copy containers
    //////////////////////////////////

    // The containers created by the shallow copy are owned by you. Remember to delete them

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



  //////////////////////////
  // Delete copy containers
  //////////////////////////

  // Deep copies. Clearing containers deletes contents including AuxStore.
  delete m_goodJet;

  // VBF study
  delete m_goodMuonForZ;
  delete m_goodMuon;
  delete m_baselineMuon;
  delete m_goodElectron;
  delete m_baselineElectron;
  delete m_goodTau;
  delete m_goodPhoton;

  delete m_met;
  delete m_metAux;


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


/*
    // print out the number of Overlap removal
    Info("finalize()", "======================================================");
    Info("finalize()", "Number overlap electrons:    %i / %i", nOverlapElectrons, nInputElectrons);
    Info("finalize()", "Number overlap muons:    %i / %i", nOverlapMuons, nInputMuons);
    Info("finalize()", "Number overlap jets:    %i / %i", nOverlapJets, nInputJets);
    Info("finalize()", "Number overlap taus:    %i / %i", nOverlapTaus, nInputTaus);
    Info("finalize()", "Number overlap photons:    %i / %i", nOverlapPhotons, nInputPhotons);
    Info("finalize()", "======================================================");
*/
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
      const xAOD::EventInfo* eventInfo, const xAOD::Vertex* primVertex){

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
    if (std::abs(muEta) >= 2.47) return EL::StatusCode::SUCCESS;

    //  if (mu.muonType()=xAOD::Muon_v1::Combined) return EL::StatusCode::SUCCESS;
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return EL::StatusCode::SUCCESS;


    // Baseline Muon
    dec_baseline(mu) = true;
    selectDec(mu) = true; // To select objects for Overlap removal


    // Muon pt cut
    if (muPt <= 25. ) return EL::StatusCode::SUCCESS;
    // Muon eta cut
    if (std::abs(muEta) >= 2.4) return EL::StatusCode::SUCCESS;

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
      tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (std::abs(d0sig) > 3.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (std::abs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(mu)) return EL::StatusCode::SUCCESS;

    // Signal Muon
    dec_signal(mu) = true;

    return EL::StatusCode::SUCCESS;

  }





  EL::StatusCode ZinvxAODAnalysis :: passMuonVBF(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo, const xAOD::Vertex* primVertex){

    dec_baseline(mu) = false;
    dec_signal(mu) = false; // For m_goodMuon container
    dec_signal_forZ(mu) = false; // For m_goodMuonForZ container where muons are the non-isolated
    selectDec(mu) = false; // To select objects for Overlap removal

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSignal. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    //if ( mu.pt() < 4000. ) return EL::StatusCode::SUCCESS;

    // Combined (CB) or Segment-tagged (ST) muons (excluding Stand-alone (SA), Calorimeter-tagged (CaloTag) muons etc..)
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return EL::StatusCode::SUCCESS;

    // Muon Calibration
    if (!m_isData){
      if(m_muonCalibrationAndSmearingTool->applyCorrection(mu) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
        Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // Muon tranverse momentum
    if (mu.pt() < m_muonPtCut ) return EL::StatusCode::SUCCESS;

    /* // eta cut is included in MuonSelectionTool
    // Muon eta cut
    if (std::abs(mu.eta()) > m_muonEtaCut) return EL::StatusCode::SUCCESS;
    */

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    //if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
    //  tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    //else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (std::abs(d0sig) > 3.0) return EL::StatusCode::SUCCESS;
    // zo cut
    //float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    float z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (std::fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    dec_baseline(mu) = true;

    // MuonSelectionTool (Loose)
    if(!m_loosemuonSelection->accept(mu)) return EL::StatusCode::SUCCESS;

    dec_signal_forZ(mu) = true; // For m_goodMuonForZ container where muons are the non-isolated

    // Isolation requirement
    if (!m_IsoToolVBF->accept(mu)) return EL::StatusCode::SUCCESS;
    // Isolation for specific muon pT range
    //m_isoPtCut = true; if (mu.pt() > m_isoMuonPtMin && mu.pt() < m_isoMuonPtMax && !m_IsoToolVBF->accept(mu)) return EL::StatusCode::SUCCESS;


    dec_signal(mu) = true; // For m_goodMuon container
    selectDec(mu) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }




  EL::StatusCode ZinvxAODAnalysis :: passElectronSelection(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo, const xAOD::Vertex* primVertex){

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
    if ( std::abs(Eta) >= 2.47 || (std::abs(Eta) >= 1.37 && std::abs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

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
    if (std::abs(d0sig) > 5.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( elec.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (std::abs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(elec)) return EL::StatusCode::SUCCESS;

    // Signal Electron
    dec_signal(elec) = true;

    return EL::StatusCode::SUCCESS;

  }




  EL::StatusCode ZinvxAODAnalysis :: passElectronVBF(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo,
      const xAOD::Vertex* primVertex){

    dec_baseline(elec) = false;
    dec_signal(elec) = false; // For m_goodElectron container
    selectDec(elec) = false; // To select objects for Overlap removal


    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    //if ( elec.pt() < 4000. ) return EL::StatusCode::SUCCESS;
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
    if ( std::abs(Eta) > m_elecEtaCut || (std::abs(Eta) > 1.37 && std::abs(Eta) < 1.52)) return EL::StatusCode::SUCCESS;
    //if ( std::abs(Eta) > m_elecEtaCut ) return EL::StatusCode::SUCCESS;

    /// pT cut
    if (elec.pt() < m_elecPtCut ) return EL::StatusCode::SUCCESS; /// veto electron

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return EL::StatusCode::SUCCESS;

    // d0 / z0 cuts applied
    // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle *tp = elec.trackParticle() ; //your input track particle from the electron
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (std::abs(d0sig) > 5.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (std::fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    dec_baseline(elec) = true;

    // LH Electron identification
    //
    // LH Electron (Loose)
    if (!m_LHToolLoose2015->accept(elec)) return EL::StatusCode::SUCCESS;
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

    // Isolation Correction Tool
    if(m_isoCorrTool->applyCorrection(elec) == CP::CorrectionCode::Error) { 
      Error("execute()", "IsolationCorrectionTool returns Error CorrectionCode");
    }

    // Isolation requirement
    if (!m_IsoToolVBF->accept(elec)) return EL::StatusCode::SUCCESS;

    dec_signal(elec) = true; // For m_goodElectron container
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
    if ( std::abs(Eta) >= 2.37 || (std::abs(Eta) >= 1.37 && std::abs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

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
    dec_signal(phot) = false; // For m_goodPhoton container
    selectDec(phot) = false; // To select objects for Overlap removal


    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passPhotonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Photon author cuts
    //if ( !(phot.author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) )
    uint16_t author =  phot.author();
    if (!(author & xAOD::EgammaParameters::AuthorPhoton) && !(author & xAOD::EgammaParameters::AuthorAmbiguous))
      return EL::StatusCode::SUCCESS;

    // MC fudge tool
    if (!m_isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = phot.caloCluster()->eta();
    double Eta = phot.caloCluster()->etaBE(2);
    //if ( std::abs(Eta) >= m_photEtaCut || (std::abs(Eta) >= 1.37 && std::abs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;
    if ( std::abs(Eta) > m_photEtaCut ) return EL::StatusCode::SUCCESS;

    // pT cut
    if (phot.pt() < m_photPtCut) return EL::StatusCode::SUCCESS; /// veto photon
    //Info("execute()", "  Selected photon pt from new Photon Container = %.2f GeV", (phot.pt() * 0.001));

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return EL::StatusCode::SUCCESS;

    dec_baseline(phot) = true;

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;
    //if (!m_photonLooseIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsoToolVBF->accept(phot)) return EL::StatusCode::SUCCESS;


    dec_signal(phot) = true; // For m_goodPhoton container
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
    dec_signal(tau) = false; // For m_goodTau container
    selectDec(tau) = false; // To select objects for Overlap removal

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Tau Smearing (for MC)
    //if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !m_isData){ // it's MC!
    if( (bool) tau.auxdata<char>("IsTruthMatched") && !m_isData){ // it's MC!
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
    dec_signal(tau) = true; // For m_goodTau container
    selectDec(tau) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }




  bool ZinvxAODAnalysis :: IsBadJet(xAOD::Jet& jet) {

    //Info("execute()", "  corrected jet pt in IsBadJet function = %.2f GeV", jet.pt() );
    //Info("execute()", "  updated jet jvt in IsBadJet function = %.2f ", cacc_jvt(jet) );

    // Pile-up
    if ( cacc_jvt(jet) < 0.64 && std::abs(jet.eta()) < 2.4 && jet.pt() < 50000. ) return false;

    // pT cut
    if ( jet.pt() < m_jetPtCut || std::abs(jet.eta()) > m_jetEtaCut) return false;

    // Jet Cleaning Tool
    dec_bad(jet) = !m_jetCleaningLoose->accept( jet );

    return dec_bad(jet);

  }


  bool ZinvxAODAnalysis :: IsSignalJet(xAOD::Jet& jet) {

    // pT, eta cut
    if ( jet.pt() < m_jetPtCut || std::abs(jet.eta()) > m_jetEtaCut ) return false;

    //bool isgoodjet = !dec_bad(jet) && (cacc_jvt(jet) > 0.64 || std::abs(jet.eta()) > 2.4 || jet.pt() > 50000.);
    bool isgoodjet = cacc_jvt(jet) >= 0.64 || std::abs(jet.eta()) >= 2.4 || jet.pt() >= 50000.;

    dec_signal(jet) = isgoodjet; // For m_goodJet container
    selectDec(jet) = isgoodjet; // To select objects for Overlap removal

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

    //int muonCount = 0;

    for (const auto& muon : muons) {
      //muonCount++;
      //Info("execute()", "-------------------------------------------------");
      //Info("execute()", " muon # : %i, Good muon pt = %.2f GeV", muonCount, (muon->pt() * 0.001));
      sf *= GetGoodMuonSF(*muon, recoSF, isoSF, ttvaSF);
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

    //int elecCount = 0;

    for (const auto& electron : electrons) {
      //elecCount++;
      //Info("execute()", "-------------------------------------------------");
      //Info("execute()", " electron # : %i, Good electron pt = %.2f GeV", elecCount, (electron->pt() * 0.001));
      if (electrons.at(0) != electron) trig_SF = false; // Apply trigSF to only leading electron
      sf *= GetGoodElectronSF(*electron, recoSF, idSF, isoSF, trig_SF);
      //Info("execute()", "  Good electron pt = %.2f GeV, trig_SF = %d, Electron SF = %.2f", (electron->pt() * 0.001), trig_SF, sf);
    }

    //Info("execute()", "  GetTotalElectronSF: Total Electron SF = %.5f ", sf );
    return sf;

  }



  int ZinvxAODAnalysis :: NumIsoTracks(const xAOD::TrackParticleContainer* inTracks,
      const xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
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
      const xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
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
      const xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
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



  float ZinvxAODAnalysis :: deltaPhi(float phi1, float phi2) {

    float dPhi = std::fabs(phi1 - phi2);

    if(dPhi > TMath::Pi())
      dPhi = TMath::TwoPi() - dPhi;

    return dPhi;

  }



  float ZinvxAODAnalysis :: deltaR(float eta1, float eta2, float phi1, float phi2) {

    float dEta = eta1 - eta2;
    float dPhi = deltaPhi(phi1,phi2);

    return TMath::Sqrt(dEta*dEta + dPhi*dPhi);

  }

   

