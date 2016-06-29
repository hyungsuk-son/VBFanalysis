#ifndef ZinvAnalysis_ZinvxAODAnalysis_H
#define ZinvAnalysis_ZinvxAODAnalysis_H

#include <EventLoop/Algorithm.h>

// Framework include(s)
#include "AsgTools/AsgTool.h"
#include "AsgTools/AsgMessaging.h"

// EDM
#include "xAODTrigger/EnergySumRoI.h"
#include "xAODTrigMissingET/TrigMissingET.h"
#include "xAODTrigMissingET/TrigMissingETContainer.h"
#include "xAODMissingET/MissingET.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonAuxContainer.h"
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/Photon.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"

// GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

// Muon
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

// Electron
#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

// Photon
#include "xAODEgamma/EgammaDefs.h"
#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"
#include "ElectronPhotonShowerShapeFudgeTool/ElectronPhotonShowerShapeFudgeTool.h"

// IsolationSelectionTool
#include "IsolationSelection/IsolationSelectionTool.h"

// Tau
#include "xAODTau/TauDefs.h"
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/TauSmearingTool.h"
#include "TauAnalysisTools/TauOverlappingElectronLLHDecorator.h"

// Jet
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetResolution/JERTool.h"
#include "JetResolution/JERSmearingTool.h"
#include "JetUncertainties/JetUncertaintiesTool.h"

// MET builder
#include "METUtilities/METMaker.h"
#include "METUtilities/CutsMETMaker.h"
#include "METUtilities/METHelpers.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETAssociationMap.h"
#include "xAODMissingET/MissingETComposition.h"

// Overlap Removal tool
#include "AssociationUtils/OverlapRemovalInit.h"
#include "AssociationUtils/DeltaROverlapTool.h"
#include "AssociationUtils/EleMuSharedTrkOverlapTool.h"
#include "AssociationUtils/EleJetOverlapTool.h"
#include "AssociationUtils/MuJetOverlapTool.h"
#include "AssociationUtils/TauLooseEleOverlapTool.h"
#include "AssociationUtils/TauLooseMuOverlapTool.h"
#include "AssociationUtils/OverlapRemovalTool.h"

// Efficiency and Scale Factor
#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"
#include "MuonEfficiencyCorrections/MuonTriggerScaleFactors.h"
#include "ElectronEfficiencyCorrection/AsgElectronEfficiencyCorrectionTool.h"
#include "JetJvtEfficiency/JetJvtEfficiency.h"
#include "TauAnalysisTools/TauEfficiencyCorrectionsTool.h"
#include "METUtilities/METSystematicsTool.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "IsolationCorrections/IsolationCorrectionTool.h"

// Event Bookkeepers
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

// Cut Flow
#include <ZinvAnalysis/BitsetCutflow.h>

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>

// include files for using the trigger tools
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"


class ZinvxAODAnalysis : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:
    // float cutValue;



    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
  public:
    // Tree *myTree; //!
    // TH1 *myHist; //!

    xAOD::TEvent *m_event; //!
    int m_eventCounter; //!
    int m_numCleanEvents; //!
    float mcEventWeight; //!
    float mcEventWeight_Zmumu; //!
    float mcEventWeight_Wmunu; //!
    float mcEventWeight_Zee; //!
    float mcEventWeight_Wenu; //!
    int mcChannelNumber; //!

    bool m_isData; //!

    // Event Channel
    bool m_isZvv; //!
    bool m_isZmumu; //!
    bool m_isWmunu; //!
    bool m_isZee; //!
    bool m_isWenu; //!

    // Cutflow
    bool m_useBitsetCutflow; //!
    int m_eventCutflow[40]; //!

    // Enable Overlap Removal tool
    bool m_doORtool; //!
    bool m_doORmanual; //!

    // Cut values
    float m_muonPtCut; //!
    float m_muonEtaCut; //!
    float m_elecPtCut; //!
    float m_elecEtaCut; //!
    float m_photPtCut; //!
    float m_photEtaCut; //!
    float m_jetPtCut; //!
    float m_jetEtaCut; //!
    float m_diJet1PtCut; //!
    float m_diJet2PtCut; //!
    float m_diJetEtaCut; //!
    float m_CJVptCut; //!
    float m_metCut; //!
    float m_mjjCut; //!
    float m_LeadLepPtCut; //!
    float m_SubLeadLepPtCut; //!
    float m_ORJETdeltaR; //!
    float m_isoMuonPtMin; //!
    float m_isoMuonPtMax; //!
    bool m_recoSF; //!
    bool m_idSF; //!
    bool m_isoSF; //!
    bool m_ttvaSF; //!

    // Some object and event counters to help roughly
    // evaluate the effects of changes in the OR tool.
    unsigned int nInputElectrons = 0; //!
    unsigned int nInputMuons = 0; //!
    unsigned int nInputJets = 0; //!
    unsigned int nInputTaus = 0; //!
    unsigned int nInputPhotons = 0; //!
    unsigned int nOverlapElectrons = 0; //!
    unsigned int nOverlapMuons = 0; //!
    unsigned int nOverlapJets = 0; //!
    unsigned int nOverlapTaus = 0; //!
    unsigned int nOverlapPhotons = 0; //!

    GoodRunsListSelectionTool *m_grl; //!


    TH1 *h_sumOfWeights; //!

    TH1 *h_met_ex; //!
    TH1 *h_met_ey; //!
    TH1 *h_met; //!
    TH1 *h_sumet; //!
    TH1 *h_met_phi; //!

    TH1 *h_emulmet_nomu_ex; //!
    TH1 *h_emulmet_nomu_ey; //!
    TH1 *h_emulmet_nomu; //!
    TH1 *h_emulsumet_nomu; //!
    TH1 *h_emulmet_nomu_phi; //!

    TH1 *h_emulmet_noelec_ex; //!
    TH1 *h_emulmet_noelec_ey; //!
    TH1 *h_emulmet_noelec; //!
    TH1 *h_emulsumet_noelec; //!
    TH1 *h_emulmet_noelec_phi; //!


    // Zinv study
    // Zvv
    TH1 *h_zvv_met; //!
    TH1 *h_zvv_njet; //!
    TH1 *h_zvv_jet1_pt; //!
    TH1 *h_zvv_jet2_pt; //!
    TH1 *h_zvv_jet3_pt; //!
    TH1 *h_zvv_jet1_phi; //!
    TH1 *h_zvv_jet2_phi; //!
    TH1 *h_zvv_jet3_phi; //!
    TH1 *h_zvv_jet1_eta; //!
    TH1 *h_zvv_jet2_eta; //!
    TH1 *h_zvv_jet3_eta; //!
    TH1 *h_zvv_jet1_rap; //!
    TH1 *h_zvv_jet2_rap; //!
    TH1 *h_zvv_jet3_rap; //!
    TH1 *h_zvv_mjj; //!
    TH1 *h_zvv_dPhijj; //!
    TH1 *h_zvv_dRjj; //!
    TH1 *h_zvv_dPhimetj1; //!
    TH1 *h_zvv_dPhimetj2; //!
    TH1 *h_zvv_dPhimetj3; //!
    TH1 *h_zvv_dPhiMinmetjet; //!


    // Zmumu
    TH1 *h_zmumu_met; //!
    TH1 *h_zmumu_emulmet_nomu; //!
    TH1 *h_zmumu_njet; //!
    TH1 *h_zmumu_jet1_pt; //!
    TH1 *h_zmumu_jet2_pt; //!
    TH1 *h_zmumu_jet3_pt; //!
    TH1 *h_zmumu_jet1_phi; //!
    TH1 *h_zmumu_jet2_phi; //!
    TH1 *h_zmumu_jet3_phi; //!
    TH1 *h_zmumu_jet1_eta; //!
    TH1 *h_zmumu_jet2_eta; //!
    TH1 *h_zmumu_jet3_eta; //!
    TH1 *h_zmumu_jet1_rap; //!
    TH1 *h_zmumu_jet2_rap; //!
    TH1 *h_zmumu_jet3_rap; //!
    TH1 *h_zmumu_mjj; //!
    TH1 *h_zmumu_dPhijj; //!
    TH1 *h_zmumu_dRjj; //!
    TH1 *h_zmumu_dPhimetj1; //!
    TH1 *h_zmumu_dPhimetj2; //!
    TH1 *h_zmumu_dPhimetj3; //!
    TH1 *h_zmumu_dPhiMinmetjet; //!
    TH1 *h_zmumu_muon1_pt; //!
    TH1 *h_zmumu_muon2_pt; //!
    TH1 *h_zmumu_muon1_phi; //!
    TH1 *h_zmumu_muon2_phi; //!
    TH1 *h_zmumu_muon1_eta; //!
    TH1 *h_zmumu_muon2_eta; //!
    TH1 *h_zmumu_mll; //!


    // Wmunu
    TH1 *h_wmunu_met; //!
    TH1 *h_wmunu_emulmet_nomu; //!
    TH1 *h_wmunu_njet; //!
    TH1 *h_wmunu_jet1_pt; //!
    TH1 *h_wmunu_jet2_pt; //!
    TH1 *h_wmunu_jet3_pt; //!
    TH1 *h_wmunu_jet1_phi; //!
    TH1 *h_wmunu_jet2_phi; //!
    TH1 *h_wmunu_jet3_phi; //!
    TH1 *h_wmunu_jet1_eta; //!
    TH1 *h_wmunu_jet2_eta; //!
    TH1 *h_wmunu_jet3_eta; //!
    TH1 *h_wmunu_jet1_rap; //!
    TH1 *h_wmunu_jet2_rap; //!
    TH1 *h_wmunu_jet3_rap; //!
    TH1 *h_wmunu_mjj; //!
    TH1 *h_wmunu_dPhijj; //!
    TH1 *h_wmunu_dRjj; //!
    TH1 *h_wmunu_dPhimetj1; //!
    TH1 *h_wmunu_dPhimetj2; //!
    TH1 *h_wmunu_dPhimetj3; //!
    TH1 *h_wmunu_dPhiMinmetjet; //!
    TH1 *h_wmunu_muon_pt; //!
    TH1 *h_wmunu_muon_phi; //!
    TH1 *h_wmunu_muon_eta; //!
    TH1 *h_wmunu_mT; //!


    // Zee
    TH1 *h_zee_met; //!
    TH1 *h_zee_emulmet_noelec; //!
    TH1 *h_zee_njet; //!
    TH1 *h_zee_jet1_pt; //!
    TH1 *h_zee_jet2_pt; //!
    TH1 *h_zee_jet3_pt; //!
    TH1 *h_zee_jet1_phi; //!
    TH1 *h_zee_jet2_phi; //!
    TH1 *h_zee_jet3_phi; //!
    TH1 *h_zee_jet1_eta; //!
    TH1 *h_zee_jet2_eta; //!
    TH1 *h_zee_jet3_eta; //!
    TH1 *h_zee_jet1_rap; //!
    TH1 *h_zee_jet2_rap; //!
    TH1 *h_zee_jet3_rap; //!
    TH1 *h_zee_mjj; //!
    TH1 *h_zee_dPhijj; //!
    TH1 *h_zee_dRjj; //!
    TH1 *h_zee_dPhimetj1; //!
    TH1 *h_zee_dPhimetj2; //!
    TH1 *h_zee_dPhimetj3; //!
    TH1 *h_zee_dPhiMinmetjet; //!
    TH1 *h_zee_electron1_pt; //!
    TH1 *h_zee_electron2_pt; //!
    TH1 *h_zee_electron1_phi; //!
    TH1 *h_zee_electron2_phi; //!
    TH1 *h_zee_electron1_eta; //!
    TH1 *h_zee_electron2_eta; //!
    TH1 *h_zee_mll; //!


    // Wenu
    TH1 *h_wenu_met; //!
    TH1 *h_wenu_emulmet_noelec; //!
    TH1 *h_wenu_njet; //!
    TH1 *h_wenu_jet1_pt; //!
    TH1 *h_wenu_jet2_pt; //!
    TH1 *h_wenu_jet3_pt; //!
    TH1 *h_wenu_jet1_phi; //!
    TH1 *h_wenu_jet2_phi; //!
    TH1 *h_wenu_jet3_phi; //!
    TH1 *h_wenu_jet1_eta; //!
    TH1 *h_wenu_jet2_eta; //!
    TH1 *h_wenu_jet3_eta; //!
    TH1 *h_wenu_jet1_rap; //!
    TH1 *h_wenu_jet2_rap; //!
    TH1 *h_wenu_jet3_rap; //!
    TH1 *h_wenu_mjj; //!
    TH1 *h_wenu_dPhijj; //!
    TH1 *h_wenu_dRjj; //!
    TH1 *h_wenu_dPhimetj1; //!
    TH1 *h_wenu_dPhimetj2; //!
    TH1 *h_wenu_dPhimetj3; //!
    TH1 *h_wenu_dPhiMinmetjet; //!
    TH1 *h_wenu_electron_pt; //!
    TH1 *h_wenu_electron_phi; //!
    TH1 *h_wenu_electron_eta; //!
    TH1 *h_wenu_mT; //!





    // trigger tools member variables
    Trig::TrigDecisionTool *m_trigDecisionTool; //!
    TrigConf::xAODConfigTool *m_trigConfigTool; //!

    // Electron and Photon
    CP::EgammaCalibrationAndSmearingTool *m_egammaCalibrationAndSmearingTool; //!

    // Muon
    CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //!
    CP::MuonSelectionTool *m_muonSelection; //!
    CP::MuonSelectionTool *m_loosemuonSelection; //!

    // Electron
    AsgElectronLikelihoodTool* m_LHToolTight2015; //!
    AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
    AsgElectronLikelihoodTool* m_LHToolLoose2015; //!

    // Photon
    AsgPhotonIsEMSelector* m_photonTightIsEMSelector; //!
    AsgPhotonIsEMSelector* m_photonMediumIsEMSelector; //!
    AsgPhotonIsEMSelector* m_photonLooseIsEMSelector; //!
    ElectronPhotonShowerShapeFudgeTool* m_electronPhotonShowerShapeFudgeTool; //!

    // IsolationSelectionTool
    CP::IsolationSelectionTool *m_IsolationSelectionTool; //!
    // IsolationSelectionTool for VBF signal
    CP::IsolationSelectionTool *m_IsoToolVBF; //!
    // Initialise Isolation Correction Tool
    CP::IsolationCorrectionTool *m_isoCorrTool; //!

    // Tau
    TauAnalysisTools::TauSelectionTool* m_tauSelTool; //!
    TauAnalysisTools::TauSmearingTool* m_tauSmearingTool; //!
    // Tau for VBF signal
    TauAnalysisTools::TauSelectionTool* m_tauSelToolVBF; //!
    // Tau
    TauAnalysisTools::TauOverlappingElectronLLHDecorator* m_tauOverlappingElectronLLHDecorator; //!

    // Jet
    JetCalibrationTool* m_jetCalibration; //!
    JetUncertaintiesTool* m_jetUncertaintiesTool; //!
    JERTool* m_jerTool; //!
    ToolHandle<IJERTool> m_jerHandle; //!
    JERSmearingTool* m_jerSmearingTool; //!
    JetVertexTaggerTool* m_jvtag; //!
    ToolHandle<IJetUpdateJvt> m_jvtagup; //!
    JetCleaningTool *m_jetCleaningTight; //!
    JetCleaningTool *m_jetCleaningLoose; //!

    // MET builder
    met::METMaker* m_metMaker; //!

    // Overlap Removal Tool
    ORUtils::ORToolBox m_toolBox; //!
    ORUtils::OverlapRemovalTool* m_orTool; //!

    // Initialise Muon Efficiency Tool
    CP::MuonEfficiencyScaleFactors* m_muonEfficiencySFTool; //!
    // Initialise Muon Isolation Tool
    CP::MuonEfficiencyScaleFactors* m_muonIsolationSFTool; //!
    // Initialise Muon TTVA Efficiency Tool
    CP::MuonEfficiencyScaleFactors* m_muonTTVAEfficiencySFTool; //!
    // Initialise Muon Trigger Scale Factor Tool
    CP::MuonTriggerScaleFactors* m_muonTriggerSFTool; //!
    // Initialise Electron Efficiency Tool
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_reco; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_id_Loose; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_id_Medium; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_id_Tight; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_iso_Loose; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_iso_Medium; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_iso_Tight; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_trigEff; //!
    AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_trigSF; //!
    CP::JetJvtEfficiency* m_jvtefficiencyTool; //!
    TauAnalysisTools::TauEfficiencyCorrectionsTool* m_tauEffTool; //!
    met::METSystematicsTool* m_metSystTool; //!
    CP::PileupReweightingTool* m_prwTool; //!


    // Cutflow
    BitsetCutflow* m_BitsetCutflow; //!


    // this is a standard constructor
    ZinvxAODAnalysis ();

    // these are the functions inherited from Algorithm
    virtual EL::StatusCode setupJob (EL::Job& job);
    virtual EL::StatusCode fileExecute ();
    virtual EL::StatusCode histInitialize ();
    virtual EL::StatusCode changeInput (bool firstFile);
    virtual EL::StatusCode initialize ();
    virtual EL::StatusCode execute ();
    virtual EL::StatusCode postExecute ();
    virtual EL::StatusCode finalize ();
    virtual EL::StatusCode histFinalize ();


    // Custom made functions

    virtual EL::StatusCode passMuonSelection(xAOD::Muon& mu,
        const xAOD::EventInfo* eventInfo,
        xAOD::Vertex* primVertex);

    virtual EL::StatusCode passMuonVBF(xAOD::Muon& mu,
        const xAOD::EventInfo* eventInfo,
        xAOD::Vertex* primVertex);

    virtual EL::StatusCode passElectronSelection(xAOD::Electron& elec,
        const xAOD::EventInfo* eventInfo,
        xAOD::Vertex* primVertex);

    virtual EL::StatusCode passElectronVBF(xAOD::Electron& elec,
        const xAOD::EventInfo* eventInfo,
        xAOD::Vertex* primVertex);

    virtual EL::StatusCode passPhotonSelection(xAOD::Photon& phot,
        const xAOD::EventInfo* eventInfo);

    virtual EL::StatusCode passTauSelection(xAOD::TauJet& tau,
        const xAOD::EventInfo* eventInfo);

    virtual EL::StatusCode passPhotonVBF(xAOD::Photon& phot,
        const xAOD::EventInfo* eventInfo);

    virtual EL::StatusCode passTauVBF(xAOD::TauJet& tau,
        const xAOD::EventInfo* eventInfo);

    bool IsBadJet(xAOD::Jet& jet);

    bool IsSignalJet(xAOD::Jet& jet);

    float GetGoodMuonSF(xAOD::Muon& mu,
        const bool recoSF, const bool isoSF, const bool ttvaSF);

    double GetTotalMuonSF(xAOD::MuonContainer& muons,
        bool recoSF, bool isoSF, bool ttvaSF);

    float GetGoodElectronSF(xAOD::Electron& elec,
        const bool recoSF, const bool idSF, const bool isoSF);

    float GetTotalElectronSF(xAOD::ElectronContainer& electrons,
        bool recoSF, bool idSF, bool isoSF);

    int NumIsoTracks(const xAOD::TrackParticleContainer* inTracks,
        xAOD::Vertex* primVertex, float Pt_Low, float Pt_High);

    int NumMuonIsoTrack(xAOD::MuonContainer* muons, const xAOD::TrackParticleContainer* inTracks,
        xAOD::Vertex* primVertex, float Pt_Low, float Pt_High);

    int NumElecIsoTrack(xAOD::ElectronContainer* electrons, const xAOD::TrackParticleContainer* inTracks,
        xAOD::Vertex* primVertex, float Pt_Low, float Pt_High);

    float DeltaPhi(float phi1, float phi2);

    float DeltaR(float eta1, float eta2, float phi1, float phi2);


    // this is needed to distribute the algorithm to the workers
    ClassDef(ZinvxAODAnalysis, 1);
};

#endif
