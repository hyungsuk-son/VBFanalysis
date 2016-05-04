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
static const std::string bJetLabel = "isBJet";
//static SG::AuxElement::Accessor<char> overlapAcc("overlaps");
ort::inputAccessor_t selectAcc(inputLabel);
ort::inputDecorator_t selectDec(inputLabel);
ort::outputAccessor_t overlapAcc(outputLabel);
ort::inputDecorator_t bJetDec(bJetLabel);
ort::objLinkAccessor_t objLinkAcc("overlapObject");
//static const bool outputPassValue = false;
//static const std::string outputLabel = outputPassValue? "passOR" : "overlaps";
//static SG::AuxElement::Decorator<char> dec_overlap(outputLabel);

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

  m_useBitsetCutflow = true;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  h_jet_selection_pt = new TH1F("h_jet_selection_pt", "Jet Signal p_{T};p_{T} (GeV)", 250, 0, 500); // Jet pt [GeV]
  wk()->addOutput (h_jet_selection_pt);

  h_refFinal_ex = new TH1F("h_refFinal_ex", "RefFinal  Missing E_{x};E_{x} (GeV)", 150, -150,  150); // RefFinal MEx [GeV]
  h_refFinal_ey = new TH1F("h_refFinal_ey", "RefFinal Missing E_{y};E_{y} (GeV)", 150, -150,  150); // RefFinal MEy [GeV]
  h_refFinal_met = new TH1F("h_refFinal_met", "RefFinal |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // RefFinal MET [GeV]
  h_refFinal_sumet = new TH1F("h_refFinal_sumet", "RefFinal Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // RefFinal SumET [GeV]
  h_refFinal_phi = new TH1F("h_refFinal_phi", "RefFinal MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // RefFinal phi [GeV]
  wk()->addOutput (h_refFinal_ex);
  wk()->addOutput (h_refFinal_ey);
  wk()->addOutput (h_refFinal_met);
  wk()->addOutput (h_refFinal_sumet);
  wk()->addOutput (h_refFinal_phi);

  // Zinv study
  // Zvv
  h_zvv_offline_met = new TH1F("h_zvv_offline_met", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  wk()->addOutput (h_zvv_offline_met);




  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZinvxAODAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
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

  m_event = wk()->xaodEvent(); // you should have already added this as described before
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
  isData = true;
  // check if the event is MC
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isData = false; // can do something with this later
  }

  // count number of events
  m_eventCounter = 0;
  m_numCleanEvents = 0;

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
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "decorrelationModel", "FULL_v1" ));  // see below for options
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
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->setProperty( "MaxEta", 2.47 ));
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
  std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150712/";
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf"));
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
  //EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2015.conf"));
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2015_CutBL.conf"));
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
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("MuonWP","FixedCutLoose"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("ElectronWP","FixedCutLoose"));
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

  // Jet
  // JES Calibration (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JES_calibration_AN1)
  const std::string name = "ZinvxAODAnalysis"; //string describing the current thread, for logging
  TString jetAlgo = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo
  TString config = "JES_MC15Prerecommendation_April2015.config"; //Path to global config used to initialize the tool
  TString calibSeq = "JetArea_Residual_Origin_EtaJES_GSC"; //String describing the calibration sequence to apply
  if (isData) calibSeq += "_Insitu";
  //Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
  m_jetCalibration = new JetCalibrationTool(name, jetAlgo, config, calibSeq, isData);
  //Initialize the tool
  EL_RETURN_CHECK("initialize()",m_jetCalibration->initializeTool(name));

  // JES uncertainty (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JES_uncertainty)
  m_jetUncertaintiesTool = new JetUncertaintiesTool("JetUncertaintiesTool");
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("JetDefinition", "AntiKt4EMTopo"));
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("MCType", "MC15"));
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("ConfigFile", "JES_2015/Prerec/PrerecJES2015_AllNuisanceParameters_25ns.config"));
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
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("isMC", !isData));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("SystematicMode", "Simple")); //"Simple" provides one NP (smearing only in MC), "Full" provides 10NPs (smearing both on data and MC)
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->initialize());

  // Configure the JVT tool.
  m_jvtag = 0;
  m_jvtag = new JetVertexTaggerTool("jvtag");
  m_jvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
  EL_RETURN_CHECK("initialize()",m_jvtag->setProperty("JVTFileName","JetMomentTools/JVTlikelihood_20140805.root"));
  EL_RETURN_CHECK("initialize()",m_jvtag->initialize());

  // Initialize and configure the jet cleaning tool
  m_jetCleaning = new JetCleaningTool("JetCleaning");
  m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
  EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty( "CutLevel", "LooseBad")); // also "TightBad"
  EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty("DoUgly", false));
  EL_RETURN_CHECK("initialize()",m_jetCleaning->initialize());

  // Initialise MET tools
  m_metMaker = new met::METMaker("METMakerTool");
  m_metMaker->setProperty( "DoRemoveMuonJets", true);
  m_metMaker->setProperty( "DoSetMuonJetEMScale", true);
  //m_metMaker->setProperty( "DoMuonEloss", true);
  //m_metMaker->setProperty( "DoIsolMuonEloss", true);
  //m_metMaker->msg().setLevel( MSG::VERBOSE ); // or DEBUG or VERBOSE
  EL_RETURN_CHECK("initialize()",m_metMaker->initialize());

  // Initialize the harmonization reccommendation tools
  const bool doTaus = false, doPhotons = false;
  const bool boostedLeptons = false;
  //EL_RETURN_CHECK("initialize()",ORUtils::recommendedTools(m_toolBox, "OverlapRemovalTool", 
  //                                                        inputLabel, outputLabel, bJetLabel, 
  //                                                        boostedLeptons, outputPassValue, 
  //                                                        doTaus, doPhotons));
  EL_RETURN_CHECK("initialize()",ORUtils::harmonizedTools(m_toolBox, "OverlapRemovalTool", 
        inputLabel, outputLabel,
        outputPassValue, doTaus, doPhotons));
  // Set message level for all tools
  //m_toolBox->setMsgLevel(MSG::DEBUG);
  // Initialize all tools
  m_orTool = static_cast<ORUtils::OverlapRemovalTool*>(m_toolBox.getMasterTool());
  m_orTool->setName("ORTool");
  EL_RETURN_CHECK("initialize()",m_toolBox.initialize());


  // Initialize Cutflow
  if (m_useBitsetCutflow)
    m_BitsetCutflow = new BitsetCutflow(wk());


  // Initialize Cutflow count array
  for (int i=0; i<20; i++) {
    m_eventCutflow[i]=0;
  }  


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
  m_eventCutflow[0]+=1;

  //----------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection in execute. Exiting." );
    return EL::StatusCode::FAILURE;
  }



  // if data check if event passes GRL
  if(isData){ // it's data!
    if(!m_grl->passRunLB(*eventInfo)){
      return EL::StatusCode::SUCCESS; // go to next event
    }
  } // end if not MC
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("GRL");
  m_eventCutflow[1]+=1;


  //------------------------------------------------------------
  // Apply event cleaning to remove events due to 
  // problematic regions of the detector, and incomplete events.
  // Apply to data.
  //------------------------------------------------------------
  // reject event if:
  if(isData){
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
  m_eventCutflow[2]+=1;


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
  m_eventCutflow[3]+=1;



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

  /*
  ///////////////////
  // For MET study //
  ///////////////////
  // [For Muon identification] filter the Muon container m_muons, placing selected muons into m_signalMuon
  xAOD::MuonContainer* m_signalMuon = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);

  // iterate over our shallow copy
  for (const auto& muon : *muonSC) { // C++11 shortcut
    //Info("execute()", "  original muon pt = %.2f GeV", ((*muonSC_itr)->pt() * 0.001)); // just to print out something
    //if(((*muonSC_itr)->eta()) > 2.5)
    //Info("execute()", "  muon eta = %.2f ", ((*muonSC_itr)->eta())); // just to print out something
    //Info("execute()", "  corrected muon pt = %.2f GeV", ((*muonSC_itr)->pt() * 0.001));
    passMuonSelection(*muon, eventInfo);

    // Signal Muon Selection
    if (passMuonSignal(*muon, eventInfo, primVertex)) {
      m_signalMuon->push_back( muon );
      //Info("execute()", "  Signal muon pt = %.2f GeV", (muon->pt() * 0.001));
    }

  } // end for loop over shallow copied muons
  */

  ///////////////////
  // For VBF study //
  ///////////////////
  xAOD::MuonContainer* m_VBFmuon = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);

  // iterate over our shallow copy
  for (const auto& muon : *muonSC) { // C++11 shortcut

    // VBF Muon Selection
    if (passMuonVBF(*muon, eventInfo, primVertex)) {
      m_VBFmuon->push_back( muon );
      //Info("execute()", "  VBF muon pt = %.2f GeV", (muon->pt() * 0.001));
    }

  } // end for loop over shallow copied muons



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

  /*
  ///////////////////
  // For MET study //
  ///////////////////
  // [For Electron identification] filter the Electron container m_electrons, placing selected electrons into m_signalElectron
  xAOD::ElectronContainer* m_signalElectron = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);

  // iterate over our shallow copy
  for (const auto& electron : *elecSC) { // C++11 shortcut
    //Info("execute()", "  original electron pt = %.2f GeV", ((*elecSC_itr)->pt() * 0.001));
    //Info("execute()", "  corrected electron pt = %.2f GeV", ((*elecSC_itr)->pt() * 0.001));
    passElectronSelection(*electron, eventInfo);

    // Signal Electron Selection
    if (passElectronSignal(*electron, eventInfo, primVertex)) {
      m_signalElectron->push_back( electron );
      //Info("execute()", "  Signal electron pt = %.2f GeV", (electron->pt() * 0.001));
    }

  } // end for loop over shallow copied electrons
  */

  ///////////////////
  // For VBF study //
  ///////////////////
  xAOD::ElectronContainer* m_VBFelectron = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);

  // iterate over our shallow copy
  for (const auto& electron : *elecSC) { // C++11 shortcut

    // VBF Electron Selection
    if (passElectronVBF(*electron, eventInfo, primVertex)) {
      m_VBFelectron->push_back( electron );
      //Info("execute()", "  VBF electron pt = %.2f GeV", (electron->pt() * 0.001));
    }

  } // end for loop over shallow copied electrons



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

  /*
  ///////////////////
  // For MET study //
  ///////////////////
  // iterate over our shallow copy
  for (const auto& photon : *photSC) { // C++11 shortcut
    //Info("execute()", "  original photon pt = %.2f GeV", ((*photSC_itr)->pt() * 0.001));
    //Info("execute()", "  corrected photon pt = %.2f GeV", ((*photSC_itr)->pt() * 0.001));
    passPhotonSelection(*photon, eventInfo);
  } // end for loop over shallow copied photons
  */

  ///////////////////
  // For VBF study //
  ///////////////////
  xAOD::PhotonContainer* m_VBFphoton = new xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  // iterate over our shallow copy
  for (const auto& photon : *photSC) { // C++11 shortcut
    // VBF Tau Selection
    if (passPhotonVBF(*photon, eventInfo)) {
      m_VBFphoton->push_back( photon );
    }
  } // end for loop over shallow copied photons



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

  /*
  ///////////////////
  // For MET study //
  ///////////////////
  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    //Info("execute()", "  original tau pt = %.2f GeV", ((*tauSC_itr)->pt() * 0.001));
    //Info("execute()", "  corrected tau pt = %.2f GeV", ((*tauSC_itr)->pt() * 0.001));
    passTauSelection(*taujet, eventInfo);
  } // end for loop over shallow copied taus
  */

  ///////////////////
  // For VBF study //
  ///////////////////
  xAOD::TauJetContainer* m_VBFtau = new xAOD::TauJetContainer(SG::VIEW_ELEMENTS);

  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    // VBF Tau Selection
    if (passTauVBF(*taujet, eventInfo)) {
      m_VBFtau->push_back( taujet );
    }

  } // end for loop over shallow copied taus



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
    //Info("execute()", "  original jet pt = %.2f GeV", ((*jetSC_itr)->pt() * 0.001));

    // According to https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15

    // JES calibration
    if ( !m_jetCalibration->applyCalibration(*jets).isSuccess() ){
      Error("execute()", "Failed to apply calibration to Jet objects. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // JVT Tool
    float newjvt = m_jvtagup->updateJvt(*jets);
    acc_jvt(*jets) = newjvt;
    //Info("execute()", "  corrected jet pt = %.2f GeV", ((*jetSC_itr)->pt() * 0.001));

    double jetPt = (jets->pt()) * 0.001; /// GeV
    double jetEta = jets->eta();

    // JES correction
    if (!isData){
      if (jetPt > 15. ){
        if ( m_jetUncertaintiesTool->applyCorrection(*jets) != CP::CorrectionCode::Ok){ // apply correction and check return code
          Error("execute()", "Failed to apply JES correction to Jet objects. Exiting." );
          return EL::StatusCode::FAILURE;
        }
      }
    }

    // JER smearing
    if (!isData){
      if ( m_jerSmearingTool->applyCorrection(*jets) != CP::CorrectionCode::Ok){ // apply correction and check return code
        Error("execute()", "Failed to apply JER smearing. Exiting. ");
        return EL::StatusCode::FAILURE;
      }
    }

    dec_baseline(*jets) = false;
    selectDec(*jets) = false;

    // pT cut
    double jetPtCut = 20.0; /// GeV
    if (jetPt > jetPtCut) {
      dec_baseline(*jets) = true;
      selectDec(*jets) = true;
    }

  } // end for loop over shallow copied jets


  //----------------
  // Overlap Removal
  //----------------

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
      double jetPt = (jets->pt()) * 0.001; /// GeV
      h_jet_selection_pt->Fill( jetPt ); // GeV

      m_goodJet->push_back( jets );
    }
  } // end for loop over shallow copied jets



  //------------------------------------
  // Event Cleaning (Jet cleaning tool)
  //------------------------------------
  if (isBadJet) return EL::StatusCode::SUCCESS;
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Jet Cleaning");
  m_eventCutflow[4]+=1;



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
    //if (selectAcc(*electron)) {
      //if(!overlapAcc(*electron)){
      //double elecPt = (electron->pt()) * 0.001; /// GeV
      //Info("execute()", "  Selected MET electron pt from shallow copy = %.2f GeV", ((*elecSC_itr)->pt() * 0.001));
      //Info("execute()", "  Selected MET electron pt from new Electron Container = %.2f GeV", (electron->pt() * 0.001));  
      m_MetElectrons.push_back( electron );
      //}
    }
  } // end for loop over shallow copied electrons
  //const xAOD::ElectronContainer* p_MetElectrons = m_MetElectrons.asDataVector();
  m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
      xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
      m_met,                                      //filling this met container
      m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
      m_metMap);                                  //and this association map


  // Photon
  //-----------------

  /// Creat New Hard Object Containers
  // [For MET building] filter the Photon container m_photons, placing selected photons into m_MetPhotons
  ConstDataVector<xAOD::PhotonContainer> m_MetPhotons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Photon>

  // iterate over our shallow copy
  for (const auto& photon : *photSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*photon)) {
    //if (selectAcc(*photon)) {
      //if(!overlapAcc(*photon)){
      //double photPt = (photon->pt()) * 0.001; /// GeV
      //Info("execute()", "  Selected MET photon pt from shallow copy = %.2f GeV", ((*photSC_itr)->pt() * 0.001));
      //Info("execute()", "  Selected MET photon pt from new Photon Container = %.2f GeV", (phot->pt() * 0.001));  
      m_MetPhotons.push_back( photon );
      //}
    }
  } // end for loop over shallow copied photons
  m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
      xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
      m_met,                                    //filling this met container
      m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
      m_metMap);                                //and this association map


  // TAUS
  //-----------------

  /// Creat New Hard Object Containers
  // [For MET building] filter the TauJet container m_taus, placing selected taus into m_MetTaus
  ConstDataVector<xAOD::TauJetContainer> m_MetTaus(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::TauJet>

  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*taujet)) {
    //if (selectAcc(*taujet)) {
      //if(!overlapAcc(*taujet)){
      //double tauPt = (taujet->pt()) * 0.001; /// GeV
      //Info("execute()", "  Selected MET tau pt from shallow copy = %.2f GeV", ((*tauSC_itr)->pt() * 0.001));
      //Info("execute()", "  Selected MET tau pt from new Tau Container = %.2f GeV", (tau->pt() * 0.001));  
      m_MetTaus.push_back( taujet );
      //}
    }
  } // end for loop over shallow copied taus
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
  for (const auto& muon : *muonSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*muon)) {
    //if (selectAcc(*muon)) {
      //if(!overlapAcc(*muon)){
      //double muPt = (muon->pt()) * 0.001; /// GeV
      //Info("execute()", "  Selected Muon pt from shallow copy = %.2f GeV", ((*muonSC_itr)->pt() * 0.001));
      //Info("execute()", "  Selected muon pt from new Muon Container = %.2f GeV", (muon->pt() * 0.001));  
      m_MetMuons.push_back( muon );
      //}
    }
  } // end for loop over shallow copied muons
  m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
      xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
      m_met,                                  //filling this met container
      m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
      m_metMap);                              //and this association map

  met::addGhostMuonsToJets(*m_muons, *jetSC);

  // JET
  //-----------------

  //Now time to rebuild jetMet and get the soft term
  //This adds the necessary soft term for both CST and TST
  //these functions create an xAODMissingET object with the given names inside the container
  m_metMaker->rebuildJetMET("RefJet",          //name of jet met
      //"SoftClus",      //name of soft cluster term met
      softTerm,          //name of soft track term met
      m_met,             //adding to this new met container
      jetSC,             //using this jet collection to calculate jet met
      m_metCore,         //core met container
      m_metMap,          //with this association map
      true);             //apply jet jvt cut


  // MET Build
  //-----------------

  //m_metMaker->rebuildTrackMET("RefJetTrk", softTerm, m_met, jetSC, m_metCore, m_metMap, true);

  //this builds the final track or cluster met sums, using systematic varied container
  //In the future, you will be able to run both of these on the same container to easily output CST and TST
  m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());


  // Fill MET RefFinal
  float refFinal_ex = -9e9;
  float refFinal_ey = -9e9;
  float refFinal_met = -9e9;
  float refFinal_sumet = -9e9;
  float refFinal_phi = -9e9;

  refFinal_ex = ((*m_met)["Final"]->mpx()) * 0.001;
  refFinal_ey = ((*m_met)["Final"]->mpy()) * 0.001;
  //refFinal_met = sqrt(refFinal_ex*refFinal_ex+refFinal_ey*refFinal_ey);
  refFinal_met = ((*m_met)["Final"]->met()) * 0.001;
  refFinal_sumet = ((*m_met)["Final"]->sumet()) * 0.001;
  //refFinal_phi = atan2(refFinal_ey, refFinal_ex);
  refFinal_phi = ((*m_met)["Final"]->phi());

  h_refFinal_ex->Fill( refFinal_ex ); // GeV
  h_refFinal_ey->Fill( refFinal_ey ); // GeV
  h_refFinal_met->Fill( refFinal_met ); // GeV
  h_refFinal_sumet->Fill( refFinal_sumet ); // GeV
  h_refFinal_phi->Fill( refFinal_phi ); // GeV



  //------------------------
  // Define Jet Properties
  // -----------------------

  TLorentzVector jet1;
  TLorentzVector jet2;

  if (m_goodJet->size() > 1) std::partial_sort(m_goodJet->begin(), m_goodJet->begin()+2, m_goodJet->end(), DescendingPt());
  float mjj = 0;
  if (m_goodJet->size() > 1) {
    jet1 = m_goodJet->at(0)->p4();
    jet2 = m_goodJet->at(1)->p4();
    auto dijet = jet1 + jet2;
    mjj = dijet.M();

    //Info("execute()", "  jet1 = %.2f GeV, jet2 = %.2f GeV", jet1.Pt()*0.001, jet2.Pt()*0.001);
    //Info("execute()", "  mjj_test = %.2f GeV", mjjtest*0.001);
  }

  // Define HT and leading jet pT
  float lead_jet_pt = 0;
  float lead_jet_eta = 0;
  float lead_jet_phi = 0;
  float secondlead_jet_pt = 0;
  float secondlead_jet_eta = 0;
  float secondlead_jet_phi = 0;
  float signalJet_ht = 0;
  bool pass_dPhijetmet = true;
  bool pass_CJV = true; // Central Jet Veto


  // N_jet >= 1
  if (m_goodJet->size() > 0) {

    // loop over the jets in the Good Jets Container
    for (const auto& signalJets : *m_goodJet) {
      double signal_jet_pt = (signalJets->pt()) * 0.001;
      double signal_jet_eta = signalJets->eta();
      double signal_jet_phi = signalJets->phi();

      // dphijetmet
      double dPhijetmet = fabs( fabs( fabs(signal_jet_phi - refFinal_phi) - TMath::Pi() ) - TMath::Pi() );
      //Info("execute()", "  delta phi = %.2f", dPhijetmet);
      if ( dPhijetmet < 0.5 ) pass_dPhijetmet = false;

      // Central Jet Veto
      if ( m_goodJet->size() > 2 && jet1.Pt() > 55e3 && jet2.Pt() > 45e3 ){
        if (m_goodJet->at(0) != signalJets && m_goodJet->at(1) != signalJets){
          //cout << "m_goodJet->at(0) = " << m_goodJet->at(0) << " signalJets = " << signalJets << endl;
          if (signal_jet_pt > 25. && fabs(signal_jet_eta) < 4.4 ) {
            if ( (jet1.Eta() > jet2.Eta()) && (signal_jet_eta < jet1.Eta() && signal_jet_eta > jet2.Eta())){
              //Info("execute()", " Jet rapidity  = %.2f, Jet eta = %.2f", signalJets->rapidity() , signal_jet_eta);
              pass_CJV = false;
            }
            if ( (jet1.Eta() < jet2.Eta()) && (signal_jet_eta > jet1.Eta() && signal_jet_eta < jet2.Eta())){
              pass_CJV = false;
            }
          }
        }
      }

      //Info("execute()", "  Zvv Signal Jet pt = %.2f GeV, eta = %.2f", signal_pt_jet, signal_eta_jet);
      signalJet_ht += signal_jet_pt;

      // Determine leading jet
      if ( signal_jet_pt > lead_jet_pt ) {
        lead_jet_pt = signal_jet_pt;
        lead_jet_eta = signal_jet_eta;
        lead_jet_phi = signal_jet_phi;
      }
    } // Jet loop

    // Determine second leading jet
    for (const auto& signalJets : *m_goodJet) {
      double signal_jet_pt = (signalJets->pt()) * 0.001;
      double signal_jet_eta = signalJets->eta();
      double signal_jet_phi = signalJets->phi();
      if ( (signal_jet_pt < lead_jet_pt) && (signal_jet_pt > secondlead_jet_pt) ) {
        secondlead_jet_pt = signal_jet_pt;
        secondlead_jet_eta = signal_jet_eta;
        secondlead_jet_phi = signal_jet_phi;
      }
    } // Jet loop

  } //N_jet >= 1

  double dPhimetjet1;
  double dPhimetjet2;

  // ------------------
  // Get isolated track
  // ------------------

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




  //-----------
  // VBF study 
  //-----------

  //------------------------------------------------
  // Z -> vv + JET EVENT SELECTION for Zinv study
  //------------------------------------------------

  if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("MET Trigger");
    m_eventCutflow[5]+=1;
    if ( refFinal_met > 150. ) {
      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("MET > 150GeV");
      m_eventCutflow[6]+=1;
      if (m_VBFelectron->size() == 0) {
        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Electron Veto");
        m_eventCutflow[7]+=1;
        if ( m_VBFmuon->size() == 0) {
          if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Muon Veto");
          m_eventCutflow[8]+=1;
          if (m_VBFtau->size() == 0) {
            if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Tau Veto");
            m_eventCutflow[9]+=1;
            if ( m_goodJet->size() > 1 ) {
              if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("At least two Jets");
              m_eventCutflow[10]+=1;
              //if ( lead_jet_pt > 55. && secondlead_jet_pt > 45. ) {
              if ( jet1.Pt() > 55e3 && jet2.Pt() > 45e3 ) {
                if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("At least two Jets[55,45]");
                m_eventCutflow[11]+=1;
                if ( mjj > 250e3 ) {
                  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("mjj > 250GeV");
                  m_eventCutflow[12]+=1;
                  if ( pass_CJV ) {
                  //if ( dPhimetjet1 > 0.5 && dPhimetjet2 > 0.5 ) {
                  //if ( fabs(jet1.Eta() - jet2.Eta() > 3.6) ) {
                    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Central Jet Veto");
                    m_eventCutflow[13]+=1;
                    if ( passIsoTrk ) {
                      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Iso Track Veto");
                      m_eventCutflow[14]+=1;
                      h_zvv_offline_met->Fill( refFinal_met ); // GeV
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


  //////////////////////////
  // Delete copy containers
  //////////////////////////

  // Deep copies. Clearing containers deletes contents including AuxStore.
  delete m_goodJet;

  // MET study
  //delete m_signalMuon;
  //delete m_signalElectron;

  // VBF study
  delete m_VBFmuon;
  delete m_VBFelectron;
  delete m_VBFtau;
  delete m_VBFphoton;


  //////////////////////////////////
  // Delete shallow copy containers
  //////////////////////////////////

  // The containers created by the shallow copy are owned by you. Remember to delete them
  delete m_met;
  delete m_metAux;

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
    if(m_jetCleaning) {
      delete m_jetCleaning;
      m_jetCleaning = 0;
    }

    /// MET Tool
    if(m_metMaker){
      delete m_metMaker;
      m_metMaker = 0;
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
    for(int i=0; i<16 ; ++i) {
      int j = i+1;
      Info("finalize()", "Event cutflow (%i) = %i", j, m_eventCutflow[i]);
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


  EL::StatusCode ZinvxAODAnalysis :: passMuonSelection(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(mu) = false;
    selectDec(mu) = false;

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    //  if (mu.muonType()=xAOD::Muon_v1::Combined) return EL::StatusCode::SUCCESS;

    // don't bother calibrating or computing WP
    double muPt = (mu.pt()) * 0.001; /// GeV
    if ( muPt < 4. ) return EL::StatusCode::SUCCESS;

    // Muon Calibration
    if (!isData){
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
    double muPtCut = 10.0; /// GeV
    if (muPt <= muPtCut ) return EL::StatusCode::SUCCESS; /// veto muon

    // Muon eta cut
    double muEta = mu.eta();
    if (fabs(muEta) >= 2.47) return EL::StatusCode::SUCCESS;

    /*
    // eta cut
    // only for Stand-alone (SA) muons, Combined (CB) muons, Segment-tagged (ST) muons
    double muEta = mu.eta();
    if (mu.muonType()==xAOD::Muon::MuonStandAlone){
    if(fabs(mu.eta())<m_mu_minEtaSA) return EL::StatusCode::SUCCESS;
    }
    else 

    if (fabs(muEta) >= 2.5) return EL::StatusCode::SUCCESS;
    */

    // Muon eta Test after Muon selection
    //if((mu.eta()) > 2.5)
    //Info("execute()", "  selected muon eta = %.2f ", (mu.eta())); // just to print out something

    dec_baseline(mu) = true;
    selectDec(mu) = true;

    return EL::StatusCode::SUCCESS;

  }


  bool ZinvxAODAnalysis :: passMuonSignal(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo,
      xAOD::Vertex* primVertex){

    dec_signal(mu) = false;
    if (!dec_baseline(mu)) return false;

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSignal. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Muon tranverse momentum
    double muPt = (mu.pt()) * 0.001;
    double ptCut = 25.;
    if (muPt <= ptCut ) return false;
    // Muon eta cut
    double muEta = mu.eta();
    if (fabs(muEta) >= 2.4) return false;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(mu)) return false;

    /*
    // d0 / z0 cuts applied 
    // do significance 
    double d0_sig = TMath::Abs(mu.primaryTrackParticle()->d0()) /
    TMath::Sqrt(mu.primaryTrackParticle()->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>3.0) return false;

    // zo cut
    double z0_vrtPVx = mu.primaryTrackParticle()->z0() +
    mu.primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(mu.eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) return false;
    */

    dec_signal(mu) = true;
    return true;

  }


  bool ZinvxAODAnalysis :: passMuonVBF(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo,
      xAOD::Vertex* primVertex){

    dec_baseline(mu) = false;
    selectDec(mu) = false;

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSignal. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double muPt = (mu.pt()) * 0.001; /// GeV
    if ( muPt < 4. ) return false;

    // Muon Calibration
    if (!isData){
      if(m_muonCalibrationAndSmearingTool->applyCorrection(mu) == CP::CorrectionCode::OutOfValidityRange){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
        Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // MuonSelectionTool (Loose)
    if(!m_loosemuonSelection->accept(mu)) return false;

    // Muon tranverse momentum
    double ptCut = 7.;
    if (muPt <= ptCut ) return false;

    // Muon eta cut
    //double muEta = mu.eta();
    //if (fabs(muEta) >= 2.47) return false;

    // Combined (CB) or Segment-tagged (ST) muons (excluding Stand-alone (SA), Calorimeter-tagged (CaloTag) muons etc..)
    //if (!(mu.muonType() == xAOD::Muon::Combined || mu.muonType() == xAOD::Muon::SegmentTagged)) return false;
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return false;

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
      tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 3.0) return false;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return false;

    // Isolation requirement
    if ((muPt > 10. && muPt < 500.) && !m_IsoToolVBF->accept(mu)) return false;

    dec_baseline(mu) = true;
    selectDec(mu) = true;

    return true;

  }





  EL::StatusCode ZinvxAODAnalysis :: passElectronSelection(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(elec) = false;
    selectDec(elec) = false;

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_identification

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double elecPt = (elec.pt()) * 0.001; /// GeV
    if ( elecPt < 4. ) return EL::StatusCode::SUCCESS;

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return false;

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
    if (!isData){
      if(m_egammaCalibrationAndSmearingTool->applyCorrection(elec) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original electron values are taken.
        Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // Eta cut
    double Eta = elec.caloCluster()->eta();
    //double Eta = elec.caloCluster()->etaBE(2);
    if ( fabs(Eta) >= 2.47 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

    // pT cut
    double elecPtCut = 10.0; /// GeV
    if (elecPt <= elecPtCut) return EL::StatusCode::SUCCESS; /// veto electron

    dec_baseline(elec) = true;
    selectDec(elec) = true;

    return EL::StatusCode::SUCCESS;

  }


  bool ZinvxAODAnalysis :: passElectronSignal(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo,
      xAOD::Vertex* primVertex){

    dec_signal(elec) = false;
    if (!dec_baseline(elec)) return false;

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_identification:

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSignal. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    /// pT cut
    double elecPt = (elec.pt()) * 0.001;
    double ptCut = 25.0; /// GeV
    if (elecPt <= ptCut) return false; /// veto electron

    /*
    // d0 / z0 cuts applied 
    // d0 significance
    const xAOD::TrackParticle *trk = elec.trackParticle();
    float d0_sig =  TMath::Abs(trk->d0())/TMath::Sqrt(trk->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>5.0) return false;
    */

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(elec)) return false;


    dec_signal(elec) = true;
    return true;

  }



  bool ZinvxAODAnalysis :: passElectronVBF(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo,
      xAOD::Vertex* primVertex){

    dec_baseline(elec) = false;
    selectDec(elec) = false;

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double elecPt = (elec.pt()) * 0.001; /// GeV
    if ( elecPt < 4. ) return false;

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return false;

    // "Please apply the identification to uncalibrated electron object. ID scale factors are to be applied to calibrated objects."
    // LH Electron identification
    //
    // LH Electron (Loose)
    bool LHlooseSel = false;
    LHlooseSel = m_LHToolLoose2015->accept(elec);
    if (!LHlooseSel) return false;
    /*
    // LH Electron (Medium)
    bool LHmediumSel = false;
    LHmediumSel = m_LHToolMedium2015->accept(elec);
    if (!LHmediumSel) return false;
    */

    //Info("execute()", "  Selected electron pt from new Electron Container = %.2f GeV", (elec.pt() * 0.001));

    // Calibration
    if (!isData){
      if(m_egammaCalibrationAndSmearingTool->applyCorrection(elec) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original electron values are taken.
        Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // Eta cut
    //double Eta = elec.caloCluster()->eta();
    double Eta = elec.caloCluster()->etaBE(2);
    //if ( fabs(Eta) >= 2.47 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return false;
    if ( fabs(Eta) > 2.47 ) return false;

    /// pT cut
    double ptCut = 7.0; /// GeV
    if (elecPt <= ptCut) return false; /// veto electron

    // d0 / z0 cuts applied
    // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle *tp = elec.trackParticle() ; //your input track particle from the electron
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 5.0) return false;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( elec.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return false;

    // Isolation requirement
    if (!m_IsoToolVBF->accept(elec)) return false;

    dec_baseline(elec) = true;
    selectDec(elec) = true;

    return true;

  }



  EL::StatusCode ZinvxAODAnalysis :: passPhotonSelection(xAOD::Photon& phot,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(phot) = false;
    selectDec(phot) = false;

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
    if (!isData){
      if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
        Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // Eta cut
    double Eta = phot.caloCluster()->eta();
    if ( fabs(Eta) >= 2.37 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

    // pT cut
    double photPt = (phot.pt()) * 0.001; /// GeV
    double photPtCut = 10.0; /// GeV
    if (photPt <= photPtCut) return EL::StatusCode::SUCCESS; /// veto photon

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return EL::StatusCode::SUCCESS;

    // MC fudge tool
    if (!isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    //if (!m_IsolationSelectionTool->accept(phot)) return EL::StatusCode::SUCCESS;


    dec_baseline(phot) = true;
    selectDec(phot) = true;

    return EL::StatusCode::SUCCESS;

  }


  bool ZinvxAODAnalysis :: passPhotonVBF(xAOD::Photon& phot,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(phot) = false;
    selectDec(phot) = false;

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passPhotonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Photon author cuts
    if ( !(phot.author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) )
      return false;


    //Info("execute()", "  Selected photon pt from new Photon Container = %.2f GeV", (phot.pt() * 0.001));

    // Calibration
    if (!isData){
      if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
        Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // Eta cut
    double Eta = phot.caloCluster()->eta();
    if ( fabs(Eta) >= 2.47 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return false;

    // pT cut
    double photPt = (phot.pt()) * 0.001; /// GeV
    double photPtCut = 20.0; /// GeV
    if (photPt <= photPtCut) return false; /// veto photon

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return false;

    // MC fudge tool
    if (!isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return false;
    //if (!m_photonLooseIsEMSelector->accept(phot)) return false;

    // Isolation requirement
    //if (!m_IsoToolVBF->accept(phot)) return false;


    dec_baseline(phot) = true;
    selectDec(phot) = true;

    return true;

  }



  EL::StatusCode ZinvxAODAnalysis :: passTauSelection(xAOD::TauJet& tau,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(tau) = false;
    selectDec(tau) = false;

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Tau Smearing (for MC)
    if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !isData){ // it's MC!
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
    selectDec(tau) = true;

    return EL::StatusCode::SUCCESS;

  }


  bool ZinvxAODAnalysis :: passTauVBF(xAOD::TauJet& tau,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(tau) = false;
    selectDec(tau) = false;

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Tau Smearing (for MC)
    if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !isData){ // it's MC!
      if(m_tauSmearingTool->applyCorrection(tau) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original tau values are taken.
        Error("execute()", "TauSmearingTool returns Error CorrectionCode");
      }
    }

    //Info("execute()", "  original tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    // TauSelectionTool (Loose for VBF)
    if(!m_tauSelToolVBF->accept(tau)) return false;

    //Info("execute()", "  Selected tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    dec_baseline(tau) = true;
    selectDec(tau) = true;

    return true;

  }




  bool ZinvxAODAnalysis :: IsBadJet(xAOD::Jet& jet) {

    //if (overlapAcc(jet)) return false;

    double jetPt = (jet.pt()) * 0.001; /// GeV

    // Pile-up
    if ( cacc_jvt(jet) < 0.59 && fabs(jet.eta()) < 2.4 && jetPt < 50.0 ) return false;

    // pT cut
    double jetPtCut = 20.0; /// GeV
    double jetEtaCut = 4.4;
    if ( jetPt <= jetPtCut || fabs(jet.eta()) >= jetEtaCut) return false; 

    // Jet Cleaning Tool
    dec_bad(jet) = !m_jetCleaning->keep( jet );

    return dec_bad(jet);

  }


  bool ZinvxAODAnalysis :: IsSignalJet(xAOD::Jet& jet) {

    //if ( !dec_baseline(jet)  || overlapAcc(jet) ) return false;
    if ( !dec_baseline(jet) ) return false;

    double jetPt = (jet.pt()) * 0.001; /// GeV
    double jetPtCut = 20.0; /// GeV
    double jetEtaCut = 4.4;

    // pT, eta cut
    if ( jetPt <= jetPtCut || fabs(jet.eta()) >= jetEtaCut) return false;

    bool isgoodjet = !dec_bad(jet) && (cacc_jvt(jet) > 0.59 || fabs(jet.eta()) > 2.4 || jetPt > 50.0);

    dec_signal(jet) = isgoodjet;

    return isgoodjet;

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

