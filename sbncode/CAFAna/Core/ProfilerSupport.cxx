#include <cstdlib>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <vector>

namespace ana
{
  /// Support for the --prof commandline option
  class ProfilerSupport
  {
  public:
    ~ProfilerSupport()
    {
      if(getenv("CPUPROFILE")){
        char tmp[] = "/tmp/XXXXXX.pdf";
        mkstemps(tmp, 4);

	const char* pd = getenv("GPERFTOOLS_DIR");
        if(!pd){
          std::cout << "Couldn't find pprof executable" << std::endl;
          return;
        }
        const std::string perfdir = pd;

        const std::string pprof = perfdir + "/bin/pprof";

	// pprof needs both of these set
	setenv("PATH", (perfdir+"/bin/:"+std::string(getenv("PATH"))).c_str(), 1);
	setenv("LD_LIBRARY_PATH", (perfdir+"/lib/:"+std::string(getenv("LD_LIBRARY_PATH"))).c_str(), 1);

        std::cout << "Creating profile. This can take some time..." << std::endl;
        const std::string cmd1 = pprof + " --pdf `which root.exe` " + std::string(getenv("CPUPROFILE")) + " > " + tmp;
        std::cout << cmd1 << std::endl;
        const int ret = system(cmd1.c_str());
        if(ret != 0) return;

        std::cout << "Displaying profile..." << std::endl;
        const std::string cmd2 = "evince "+std::string(tmp);
        std::cout << cmd2 << std::endl;
        system(cmd2.c_str());
      }
    }
  };

  namespace{
    // Run destructor at end of process
    static ProfilerSupport ps;
  }
}
