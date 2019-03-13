// -*- mode: c++; c-basic-offset: 2; -*-
void load(std::string lib)
{
  std::cout << "." << std::flush;
  int ret = gSystem->Load(("lib"+lib).c_str());
  // In case of error, exit immediately with the error clearly showing, instead
  // of with a confusing secondary error about a page of output later.
  if(ret != 0){
    std::cout << std::endl << "gSystem->Load(lib"+lib+") failed with code " << ret << std::endl;
    exit(ret);
  }
}

void load_libs()
{
  // All the CINT exception handler does is obfuscate the stack. With this,
  // uncaught exceptions immediately show a useful backtrace under gdb.
  //  G__SetCatchException(0);

  TString qsrt = gSystem->Getenv("SRT_QUAL");
  TString qmrb = gSystem->Getenv("MRB_QUALS");
  // Mirror the optimization settings in use elsewhere
  if( qsrt.Contains("debug") || qmrb.Contains("debug") ) {
    gSystem->SetAclicMode(TSystem::kDebug);
  }
  else{
    gSystem->SetAclicMode(TSystem::kOpt);
  }

  // This magic incantation prevents ROOT doing slow cleanup work in
  // TList::RecursiveRemove() under ~TH1(). It also tends to lead to shutdown
  // crashes. HistCache now avoids almost all of that histogram destruction, so
  // now we can leave this at the default setting and get both speed and
  // stability.
  gROOT->SetMustClean(false);

  // Colorize error messages. Would be better if we could just pick up the
  // flags novasoft uses, but they don't seem to be in any env var.
  gSystem->SetFlagsDebug(TString(gSystem->GetFlagsDebug())+" -fdiagnostics-color=auto");
  gSystem->SetFlagsOpt(TString(gSystem->GetFlagsOpt())+" -fdiagnostics-color=auto -UNDEBUG"); // match gcc's maxopt behaviour of retaining assert()


  // Include path
  TString includes = "-I$SRT_PRIVATE_CONTEXT/include/ -I$SRT_PUBLIC_CONTEXT/include/ -I$ROOTSYS/include -I$NUTOOLS_INC -I$GENIE_INC/GENIE/";

  if (qmrb == "") {    // This is the SRT build
    // List of libraries to load. Dependency order.
    const std::vector<std::string> libs =
      {
        "Minuit2", // CAFReweight pulls in Genie which pulls in ROOT geometry
        // "Cintex",
        "StandardRecord",
        // "StandardRecord_dict",
        "CAFAnaCore",
        "CAFAnaVars",
        "CAFAnaCuts",
        "CAFAnaSysts",
        // "CAFAnaUnfold",
        "CAFAnaDecomp",
        "CAFAnaExtrap",
        "CAFAnaPrediction",
        //        "CAFAnaExperiment",
        // "CAFAnaFC",
        "CAFAnaAnalysis",
        // "CAFAnaXSec",
        // "CAFAna",
        //        "OscLibFunc",
        // "MCReweightFunc",
        //        "ifdh"
      };

    // Actually load the libraries
    std::cout << "Loading libraries";
    for(const std::string& lib: libs) load(lib);
    std::cout << std::endl;
  }
  else {   // This is the MRB build
    // List of libraries to load. Dependency order.
    const std::vector<std::string> libs =
      {
        "Geom", "Tree", "Minuit2", // CAFReweight pulls in Genie which pulls in ROOT geometry
        "Cintex",
        "StandardRecord", "StandardRecord_dict",
        "CAFAna",
        "OscLibFunc",
	"MCReweightFunc",
        "ifdh"
      };

    // Actually load the libraries
    std::cout << "Loading libraries (mrb)";
    for(const std::string& lib: libs) load(lib);
    std::cout << std::endl;

    // MRB requires an extra include path
    includes += " -I$NOVASOFT_INC";
  }

  // Magic incantation to get Reflex-based dictionaries converted into the
  // CINT-based information that ROOT needs
  //  Cintex::Enable();


  gSystem->SetIncludePath(includes);

  // Doesn't seem to work
  //  gSystem->Setenv("IFDH_DEBUG", "0"); // shut ifdh up

  // Pick up standard NOvA style
  gROOT->Macro("$SRT_PUBLIC_CONTEXT/Utilities/rootlogon.C");
  gROOT->ForceStyle();

  TRint* rint = dynamic_cast<TRint*>(gApplication);
  if(rint) rint->SetPrompt("cafe [%d] ");


  // Do this last so that the profiler library is unloaded before our main
  // libraries, meaning that the code in ProfilerSupport runs at the right time
  // (after the profile file is made). This does mean that we aren't able to
  // profile library loading/startup. The ideal solution would be to make a
  // mini-library for just ProfilerSupport.
  if(gSystem->Getenv("CPUPROFILE")){
    // We were passed the --prof option
    const std::string cpuprof = gSystem->Getenv("CPUPROFILE");

    std::cout << "Profiling enabled." << std::endl;

    if(!qsrt.Contains("debug") && !qmrb.Contains("debug")){
      std::cout << "Note: profiling works much better in debug mode." << std::endl;
    }

    const char* pd = getenv("GPERFTOOLS_DIR");
    if(pd){
      gSystem->Load((std::string(pd)+"/lib/libprofiler.so").c_str());
      // Apparently the library load manages to corrupt the env var? Put it
      // back.
      gSystem->Setenv("CPUPROFILE", cpuprof.c_str());
    }
    else{
      std::cout << "Couldn't find gperftools library" << std::endl;
    }
  }
}
