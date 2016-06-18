#include "xAODRootAccess/Init.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopGrid/GridDriver.h"
#include "EventLoopGrid/GridWorker.h"
#include <TSystem.h>
#include "SampleHandler/Sample.h"
#include "SampleHandler/SampleGrid.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/MetaFields.h"
#include "SampleHandler/MetaObject.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

#include "ZinvAnalysis/ZinvxAODAnalysis.h"

int main( int argc, char* argv[] ) {

  // Take the submit directory from the input if provided:
  std::string submitDir = "submitDir";
  if( argc > 1 ) submitDir = argv[ 1 ];

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // Construct the samples to run on:
  SH::SampleHandler sh;

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
//  const char* inputFilePath = gSystem->ExpandPathName ("$ALRB_TutorialData/r6630/");
//  SH::ScanDir().sampleDepth(1).samplePattern("AOD.05352803._000031.pool.root.1").scan(sh, inputFilePath);
//  const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/user/h/hson/Work/xAOD_data");
//  SH::DiskListLocal list (inputFilePath);
//  SH::scanDir (sh, list, "data15_13TeV.00271048.physics_Main.merge.AOD.f611_m1463._lb0408._0001.1"); // specifying one particular file for testing
  // If you want to use grid datasets the easiest option for discovery is scanDQ2
//  SH::scanDQ2 (sh, "data15_13TeV.00270816.physics_Main.merge.AOD.f611_m1463");
  // If you know the name of all your grid datasets you can also skip the dq2-ls step and add the datasets directly
//  SH::addGrid (sh, "data15_13TeV.00281411.physics_Main.merge.AOD.f629_m1504");
//  SH::addGrid (sh, "data15_13TeV.00282784.physics_Main.merge.AOD.f640_m1511");
  SH::addGrid (sh, "mc15_13TeV.304018.Sherpa_CT10_Znunu2JetsEW1JetQCD15GeV.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.304015.Sherpa_CT10_Wenu2JetsEW1JetQCD15GeV.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.304016.Sherpa_CT10_Wmunu2JetsEW1JetQCD15GeV.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.304017.Sherpa_CT10_Wtaunu2JetsEW1JetQCD15GeV.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.304019.Sherpa_CT10_Zee2JetsEW1JetQCD15GeVM40.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.304020.Sherpa_CT10_Zmumu2JetsEW1JetQCD15GeVM40.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.304021.Sherpa_CT10_Ztautau2JetsEW1JetQCD15GeVM40.merge.DAOD_EXOT5.e4523_s2608_r7326_r6282_p2495");
   
  SH::addGrid (sh, "mc15_13TeV.363412.Sherpa_NNPDF30NNLO_Znunu_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363413.Sherpa_NNPDF30NNLO_Znunu_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363414.Sherpa_NNPDF30NNLO_Znunu_Pt0_70_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363415.Sherpa_NNPDF30NNLO_Znunu_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363416.Sherpa_NNPDF30NNLO_Znunu_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363417.Sherpa_NNPDF30NNLO_Znunu_Pt70_140_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363418.Sherpa_NNPDF30NNLO_Znunu_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363419.Sherpa_NNPDF30NNLO_Znunu_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363420.Sherpa_NNPDF30NNLO_Znunu_Pt140_280_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363421.Sherpa_NNPDF30NNLO_Znunu_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363422.Sherpa_NNPDF30NNLO_Znunu_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363423.Sherpa_NNPDF30NNLO_Znunu_Pt280_500_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363424.Sherpa_NNPDF30NNLO_Znunu_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363425.Sherpa_NNPDF30NNLO_Znunu_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363426.Sherpa_NNPDF30NNLO_Znunu_Pt500_700_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363427.Sherpa_NNPDF30NNLO_Znunu_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363428.Sherpa_NNPDF30NNLO_Znunu_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363429.Sherpa_NNPDF30NNLO_Znunu_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363430.Sherpa_NNPDF30NNLO_Znunu_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363431.Sherpa_NNPDF30NNLO_Znunu_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363432.Sherpa_NNPDF30NNLO_Znunu_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363433.Sherpa_NNPDF30NNLO_Znunu_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363434.Sherpa_NNPDF30NNLO_Znunu_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363435.Sherpa_NNPDF30NNLO_Znunu_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
   
  SH::addGrid (sh, "mc15_13TeV.363460.Sherpa_NNPDF30NNLO_Wenu_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363461.Sherpa_NNPDF30NNLO_Wenu_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363462.Sherpa_NNPDF30NNLO_Wenu_Pt0_70_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363463.Sherpa_NNPDF30NNLO_Wenu_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363464.Sherpa_NNPDF30NNLO_Wenu_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363465.Sherpa_NNPDF30NNLO_Wenu_Pt70_140_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363466.Sherpa_NNPDF30NNLO_Wenu_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363467.Sherpa_NNPDF30NNLO_Wenu_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363468.Sherpa_NNPDF30NNLO_Wenu_Pt140_280_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363469.Sherpa_NNPDF30NNLO_Wenu_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363470.Sherpa_NNPDF30NNLO_Wenu_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363471.Sherpa_NNPDF30NNLO_Wenu_Pt280_500_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363472.Sherpa_NNPDF30NNLO_Wenu_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363473.Sherpa_NNPDF30NNLO_Wenu_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363474.Sherpa_NNPDF30NNLO_Wenu_Pt500_700_BFilter.merge.DAOD_EXOT5.e4771_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363475.Sherpa_NNPDF30NNLO_Wenu_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363476.Sherpa_NNPDF30NNLO_Wenu_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363477.Sherpa_NNPDF30NNLO_Wenu_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363478.Sherpa_NNPDF30NNLO_Wenu_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363479.Sherpa_NNPDF30NNLO_Wenu_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363480.Sherpa_NNPDF30NNLO_Wenu_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363481.Sherpa_NNPDF30NNLO_Wenu_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363482.Sherpa_NNPDF30NNLO_Wenu_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363483.Sherpa_NNPDF30NNLO_Wenu_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.363436.Sherpa_NNPDF30NNLO_Wmunu_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363437.Sherpa_NNPDF30NNLO_Wmunu_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363438.Sherpa_NNPDF30NNLO_Wmunu_Pt0_70_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363439.Sherpa_NNPDF30NNLO_Wmunu_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363440.Sherpa_NNPDF30NNLO_Wmunu_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363441.Sherpa_NNPDF30NNLO_Wmunu_Pt70_140_BFilter.merge.DAOD_EXOT5.e4771_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363442.Sherpa_NNPDF30NNLO_Wmunu_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363443.Sherpa_NNPDF30NNLO_Wmunu_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363444.Sherpa_NNPDF30NNLO_Wmunu_Pt140_280_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363445.Sherpa_NNPDF30NNLO_Wmunu_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363446.Sherpa_NNPDF30NNLO_Wmunu_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363447.Sherpa_NNPDF30NNLO_Wmunu_Pt280_500_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363448.Sherpa_NNPDF30NNLO_Wmunu_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363449.Sherpa_NNPDF30NNLO_Wmunu_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363450.Sherpa_NNPDF30NNLO_Wmunu_Pt500_700_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363451.Sherpa_NNPDF30NNLO_Wmunu_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363452.Sherpa_NNPDF30NNLO_Wmunu_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363453.Sherpa_NNPDF30NNLO_Wmunu_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363454.Sherpa_NNPDF30NNLO_Wmunu_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363455.Sherpa_NNPDF30NNLO_Wmunu_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363456.Sherpa_NNPDF30NNLO_Wmunu_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363457.Sherpa_NNPDF30NNLO_Wmunu_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363458.Sherpa_NNPDF30NNLO_Wmunu_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363459.Sherpa_NNPDF30NNLO_Wmunu_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4715_s2726_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.363331.Sherpa_NNPDF30NNLO_Wtaunu_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363332.Sherpa_NNPDF30NNLO_Wtaunu_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363333.Sherpa_NNPDF30NNLO_Wtaunu_Pt0_70_BFilter.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363334.Sherpa_NNPDF30NNLO_Wtaunu_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363335.Sherpa_NNPDF30NNLO_Wtaunu_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363336.Sherpa_NNPDF30NNLO_Wtaunu_Pt70_140_BFilter.merge.DAOD_EXOT5.e4779_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363337.Sherpa_NNPDF30NNLO_Wtaunu_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363338.Sherpa_NNPDF30NNLO_Wtaunu_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363339.Sherpa_NNPDF30NNLO_Wtaunu_Pt140_280_BFilter.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363340.Sherpa_NNPDF30NNLO_Wtaunu_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363341.Sherpa_NNPDF30NNLO_Wtaunu_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4779_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363342.Sherpa_NNPDF30NNLO_Wtaunu_Pt280_500_BFilter.merge.DAOD_EXOT5.e4779_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363343.Sherpa_NNPDF30NNLO_Wtaunu_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363344.Sherpa_NNPDF30NNLO_Wtaunu_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363345.Sherpa_NNPDF30NNLO_Wtaunu_Pt500_700_BFilter.merge.DAOD_EXOT5.e4779_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363346.Sherpa_NNPDF30NNLO_Wtaunu_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363347.Sherpa_NNPDF30NNLO_Wtaunu_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363348.Sherpa_NNPDF30NNLO_Wtaunu_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4779_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363349.Sherpa_NNPDF30NNLO_Wtaunu_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363350.Sherpa_NNPDF30NNLO_Wtaunu_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363351.Sherpa_NNPDF30NNLO_Wtaunu_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4779_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363352.Sherpa_NNPDF30NNLO_Wtaunu_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363353.Sherpa_NNPDF30NNLO_Wtaunu_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363354.Sherpa_NNPDF30NNLO_Wtaunu_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4709_s2726_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.363388.Sherpa_NNPDF30NNLO_Zee_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363389.Sherpa_NNPDF30NNLO_Zee_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363390.Sherpa_NNPDF30NNLO_Zee_Pt0_70_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363391.Sherpa_NNPDF30NNLO_Zee_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363392.Sherpa_NNPDF30NNLO_Zee_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363393.Sherpa_NNPDF30NNLO_Zee_Pt70_140_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363394.Sherpa_NNPDF30NNLO_Zee_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363395.Sherpa_NNPDF30NNLO_Zee_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363396.Sherpa_NNPDF30NNLO_Zee_Pt140_280_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363397.Sherpa_NNPDF30NNLO_Zee_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363398.Sherpa_NNPDF30NNLO_Zee_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363399.Sherpa_NNPDF30NNLO_Zee_Pt280_500_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363400.Sherpa_NNPDF30NNLO_Zee_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363401.Sherpa_NNPDF30NNLO_Zee_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363402.Sherpa_NNPDF30NNLO_Zee_Pt500_700_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363403.Sherpa_NNPDF30NNLO_Zee_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363404.Sherpa_NNPDF30NNLO_Zee_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363405.Sherpa_NNPDF30NNLO_Zee_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363406.Sherpa_NNPDF30NNLO_Zee_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363407.Sherpa_NNPDF30NNLO_Zee_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363408.Sherpa_NNPDF30NNLO_Zee_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363409.Sherpa_NNPDF30NNLO_Zee_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363410.Sherpa_NNPDF30NNLO_Zee_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363411.Sherpa_NNPDF30NNLO_Zee_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.363364.Sherpa_NNPDF30NNLO_Zmumu_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363365.Sherpa_NNPDF30NNLO_Zmumu_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363366.Sherpa_NNPDF30NNLO_Zmumu_Pt0_70_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363367.Sherpa_NNPDF30NNLO_Zmumu_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363368.Sherpa_NNPDF30NNLO_Zmumu_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363369.Sherpa_NNPDF30NNLO_Zmumu_Pt70_140_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363370.Sherpa_NNPDF30NNLO_Zmumu_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363371.Sherpa_NNPDF30NNLO_Zmumu_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363372.Sherpa_NNPDF30NNLO_Zmumu_Pt140_280_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363373.Sherpa_NNPDF30NNLO_Zmumu_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363374.Sherpa_NNPDF30NNLO_Zmumu_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363375.Sherpa_NNPDF30NNLO_Zmumu_Pt280_500_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363376.Sherpa_NNPDF30NNLO_Zmumu_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363377.Sherpa_NNPDF30NNLO_Zmumu_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363378.Sherpa_NNPDF30NNLO_Zmumu_Pt500_700_BFilter.merge.DAOD_EXOT5.e4772_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363379.Sherpa_NNPDF30NNLO_Zmumu_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363380.Sherpa_NNPDF30NNLO_Zmumu_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363381.Sherpa_NNPDF30NNLO_Zmumu_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363382.Sherpa_NNPDF30NNLO_Zmumu_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363383.Sherpa_NNPDF30NNLO_Zmumu_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363384.Sherpa_NNPDF30NNLO_Zmumu_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363385.Sherpa_NNPDF30NNLO_Zmumu_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363386.Sherpa_NNPDF30NNLO_Zmumu_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363387.Sherpa_NNPDF30NNLO_Zmumu_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4716_s2726_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.363361.Sherpa_NNPDF30NNLO_Ztautau_Pt0_70_CVetoBVeto.merge.DAOD_EXOT5.e4689_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363362.Sherpa_NNPDF30NNLO_Ztautau_Pt0_70_CFilterBVeto.merge.DAOD_EXOT5.e4689_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363363.Sherpa_NNPDF30NNLO_Ztautau_Pt0_70_BFilter.merge.DAOD_EXOT5.e4743_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363102.Sherpa_NNPDF30NNLO_Ztautau_Pt70_140_CVetoBVeto.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363103.Sherpa_NNPDF30NNLO_Ztautau_Pt70_140_CFilterBVeto.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363104.Sherpa_NNPDF30NNLO_Ztautau_Pt70_140_BFilter.merge.DAOD_EXOT5.e4792_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363105.Sherpa_NNPDF30NNLO_Ztautau_Pt140_280_CVetoBVeto.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363106.Sherpa_NNPDF30NNLO_Ztautau_Pt140_280_CFilterBVeto.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363107.Sherpa_NNPDF30NNLO_Ztautau_Pt140_280_BFilter.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363108.Sherpa_NNPDF30NNLO_Ztautau_Pt280_500_CVetoBVeto.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363109.Sherpa_NNPDF30NNLO_Ztautau_Pt280_500_CFilterBVeto.merge.DAOD_EXOT5.e4792_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363110.Sherpa_NNPDF30NNLO_Ztautau_Pt280_500_BFilter.merge.DAOD_EXOT5.e4792_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363111.Sherpa_NNPDF30NNLO_Ztautau_Pt500_700_CVetoBVeto.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363112.Sherpa_NNPDF30NNLO_Ztautau_Pt500_700_CFilterBVeto.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363113.Sherpa_NNPDF30NNLO_Ztautau_Pt500_700_BFilter.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363114.Sherpa_NNPDF30NNLO_Ztautau_Pt700_1000_CVetoBVeto.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363115.Sherpa_NNPDF30NNLO_Ztautau_Pt700_1000_CFilterBVeto.merge.DAOD_EXOT5.e4792_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363116.Sherpa_NNPDF30NNLO_Ztautau_Pt700_1000_BFilter.merge.DAOD_EXOT5.e4742_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363117.Sherpa_NNPDF30NNLO_Ztautau_Pt1000_2000_CVetoBVeto.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363118.Sherpa_NNPDF30NNLO_Ztautau_Pt1000_2000_CFilterBVeto.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363119.Sherpa_NNPDF30NNLO_Ztautau_Pt1000_2000_BFilter.merge.DAOD_EXOT5.e4666_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363120.Sherpa_NNPDF30NNLO_Ztautau_Pt2000_E_CMS_CVetoBVeto.merge.DAOD_EXOT5.e4690_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363121.Sherpa_NNPDF30NNLO_Ztautau_Pt2000_E_CMS_CFilterBVeto.merge.DAOD_EXOT5.e4690_s2726_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.363122.Sherpa_NNPDF30NNLO_Ztautau_Pt2000_E_CMS_BFilter.merge.DAOD_EXOT5.e4792_s2726_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_EXOT5.e3698_s2608_s2183_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410007.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_allhad.merge.DAOD_EXOT5.e4135_s2608_s2183_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410011.PowhegPythiaEvtGen_P2012_singletop_tchan_lept_top.merge.DAOD_EXOT5.e3824_s2608_s2183_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410012.PowhegPythiaEvtGen_P2012_singletop_tchan_lept_antitop.merge.DAOD_EXOT5.e3824_s2608_s2183_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410013.PowhegPythiaEvtGen_P2012_Wt_inclusive_top.merge.DAOD_EXOT5.e3753_s2608_s2183_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410014.PowhegPythiaEvtGen_P2012_Wt_inclusive_antitop.merge.DAOD_EXOT5.e3753_s2608_s2183_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410025.PowhegPythiaEvtGen_P2012_SingleTopSchan_noAllHad_top.merge.DAOD_EXOT5.e3998_s2608_s2183_r7326_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.410026.PowhegPythiaEvtGen_P2012_SingleTopSchan_noAllHad_antitop.merge.DAOD_EXOT5.e3998_s2608_s2183_r7326_r6282_p2495");
     
  SH::addGrid (sh, "mc15_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.merge.DAOD_EXOT5.e3569_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.merge.DAOD_EXOT5.e3569_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.merge.DAOD_EXOT5.e3668_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361023.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W.merge.DAOD_EXOT5.e3668_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361024.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4W.merge.DAOD_EXOT5.e3668_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361025.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5W.merge.DAOD_EXOT5.e3668_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361026.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ6W.merge.DAOD_EXOT5.e3569_s2608_s2183_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361027.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ7W.merge.DAOD_EXOT5.e3668_s2608_s2183_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361028.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ8W.merge.DAOD_EXOT5.e3569_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361029.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ9W.merge.DAOD_EXOT5.e3569_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361030.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ10W.merge.DAOD_EXOT5.e3569_s2576_s2132_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361031.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ11W.merge.DAOD_EXOT5.e3569_s2608_s2183_r7267_r6282_p2495");
  SH::addGrid (sh, "mc15_13TeV.361032.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ12W.merge.DAOD_EXOT5.e3668_s2608_s2183_r7267_r6282_p2495");


  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );
  sh.setMetaString("nc_grid_filter", "*AOD*");
//  sh.setMetaString("nc_grid_filter", "*");

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
//  job.options()->setDouble (EL::Job::optMaxEvents, 500);
  job.options()->setDouble (EL::Job::optRetries, 30);
  job.options()->setDouble (EL::Job::optCacheSize, 10*1024*1024);
  job.options()->setDouble (EL::Job::optCacheLearnEntries, 50);

/*  
  // For ntuple
  // define an output and an ntuple associated to that output
  EL::OutputStream output  ("myOutput");
  job.outputAdd (output);
  EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
  job.algsAdd (ntuple);
  //---------------------------------------------------------
*/

  // Add our analysis to the job:
  ZinvxAODAnalysis* alg = new ZinvxAODAnalysis();
  job.algsAdd( alg );

/*
  // For ntuple
  // Let your algorithm know the name of the output file
  alg->outputName = "myOutput"; // give the name of the output to our algorithm
  //-----------------------------------------------------------------------------
*/

  // Run the job using the local/direct driver:
//  EL::DirectDriver driver; //local
  EL::PrunDriver driver;  //grid
//  EL::GridDriver driver; //grid in the background

  driver.options()->setString("nc_outputSampleName", "user.hson.mc15_13TeV.DAOD_EXOT5.06192016.%in:name[2]%.%in:name[6]%"); //For PrunDriver
//  driver.outputSampleName = "user.hson.gridtest1.11142015.%in:name[2]%.%in:name[6]%"; //For GridDriver
//  driver.options()->setDouble("nc_nFiles", 1); // FOR TESTING!
//  driver.options()->setDouble("nc_nFilesPerJob", 1);
//  driver.options()->setDouble(EL::Job::optGridNFilesPerJob, 1);
//  driver.options()->setString("nc_excludedSite", "ANALY_SCINET,ANALY_VICTORIA,ANALY_CERN_CLOUD,ANALY_IN2P3-CC,ANALY_LAPP,ANALY_CONNECT_SHORT,ANALY_SFU,ANALY_CONNECT,ANALY_RAL_SL6,ANALY_GRIF-LPNHE,ANALY_HU_ATLAS_Tier2,ANALY_OU_OCHEP_SWT2,ANALY_IFIC,ANALY_ECDF_SL6");
//  driver.options()->setString("nc_excludedSite", "ANALY_INFN-NAPOLI-RECAS,ANALY_INFN-NAPOLI,ANALY_DESY-HH,ANALY_GRIF-IRFU,ANALY_AUSTRALIA,ANALY_SFU,ANALY_SCINET,ANALY_CPPM,ANALY_SiGNET,ANALY_LPC,ANALY_NSC,ANALY_CONNECT,ANALY_MWT2_SL6,ANALY_BU_ATLAS_Tier2_SL6,ANALY_wuppertalprod,ANALY_ARNES,ANALY_SLAC_SHORT_1HR,ANALY_SLAC,ANALY_RAL_SL6,ANALY_INFN-MILANO-ATLASC");
//  driver.options()->setString("nc_excludedSite", "ANALY_DCSC,ANALY_SiGNET");
//  driver.options()->setString("nc_site", "ANALY_CERN_SHORT,ANALY_CERN_SLC6,ANALY_PIC_SL6,ANALY_SARA"); // The Reflex dictionary build only works on a few sites
//  driver.options()->setString("nc_site", "ANALY_CERN_SLC6"); // The Reflex dictionary build only works on a few sites
//  driver.options()->setDouble(EL::Job::optGridMemory,10240); // 10 GB

//  driver.submit( job, submitDir );  // with monitoring
  driver.submitOnly(job, submitDir);  //without monitoring

  return 0;
}
