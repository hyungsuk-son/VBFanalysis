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
//period D
  SH::addGrid (sh, "data15_13TeV.00276262.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276329.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276336.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276416.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276511.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276689.physics_Main.merge.DAOD_EXOT5.f623_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276778.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276790.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276952.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
  SH::addGrid (sh, "data15_13TeV.00276954.physics_Main.merge.DAOD_EXOT5.f620_m1480_p2524");
   
//period E
  SH::addGrid (sh, "data15_13TeV.00278880.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00278912.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00278968.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279169.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279259.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279279.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279284.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279345.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279515.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279598.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279685.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279764.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279813.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279867.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
  SH::addGrid (sh, "data15_13TeV.00279928.physics_Main.merge.DAOD_EXOT5.f628_m1497_p2524");
   
//period F
  SH::addGrid (sh, "data15_13TeV.00279932.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00279984.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280231.physics_Main.merge.DAOD_EXOT5.f630_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280319.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280368.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
   
//period G
  SH::addGrid (sh, "data15_13TeV.00280423.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280464.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280500.physics_Main.merge.DAOD_EXOT5.f631_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280520.physics_Main.merge.DAOD_EXOT5.f632_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280614.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280673.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280753.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280853.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280862.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280950.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00280977.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00281070.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00281074.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00281075.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");

//period H
  SH::addGrid (sh, "data15_13TeV.00281317.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00281385.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
  SH::addGrid (sh, "data15_13TeV.00281411.physics_Main.merge.DAOD_EXOT5.f629_m1504_p2524");
   
//period J
  SH::addGrid (sh, "data15_13TeV.00282625.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00282631.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00282712.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00282784.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00282992.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00283074.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00283155.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00283270.physics_Main.merge.DAOD_EXOT5.f640_m1511_p2524");
  SH::addGrid (sh, "data15_13TeV.00283429.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00283608.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00283780.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284006.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284154.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284213.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284285.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284420.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284427.physics_Main.merge.DAOD_EXOT5.f643_m1518_p2524");
  SH::addGrid (sh, "data15_13TeV.00284484.physics_Main.merge.DAOD_EXOT5.f644_m1518_p2524");


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

  driver.options()->setString("nc_outputSampleName", "user.hson.data15_13TeV.DAOD.06192016.%in:name[2]%.%in:name[6]%"); //For PrunDriver
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
