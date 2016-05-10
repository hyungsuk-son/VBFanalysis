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

    bool m_isData; //!

    // Event Channel
    bool m_isZvv; //!
    bool m_isZmumu; //!
    bool m_isZee; //!

    // Cutflow
    bool m_useBitsetCutflow; //!
    int m_eventCutflow[20]; //!

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


    TH1 *h_jet_selection_pt; //!

    TH1 *h_met_ex; //!
    TH1 *h_met_ey; //!
    TH1 *h_met; //!
    TH1 *h_sumet; //!
    TH1 *h_met_phi; //!

    // Zinv study
    TH1 *h_zvv_offline_met; //!


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

    // Tau
    TauAnalysisTools::TauSelectionTool* m_tauSelTool; //!
    TauAnalysisTools::TauSmearingTool* m_tauSmearingTool; //!
    // Tau for VBF signal
    TauAnalysisTools::TauSelectionTool* m_tauSelToolVBF; //!

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
        const xAOD::EventInfo* eventInfo);

    bool passMuonSignal(xAOD::Muon& mu,
        const xAOD::EventInfo* eventInfo,
        xAOD::Vertex* primVertex);

    virtual EL::StatusCode passMuonVBF(xAOD::Muon& mu,
        const xAOD::EventInfo* eventInfo,
        xAOD::Vertex* primVertex);

    virtual EL::StatusCode passElectronSelection(xAOD::Electron& elec,
        const xAOD::EventInfo* eventInfo);

    bool passElectronSignal(xAOD::Electron& elec,
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
