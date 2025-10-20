////////////////////////////////////////////////////////////////////////
/// \file  CORSIKAGenSBN_module.cc
/// \brief Generator for cosmic-ray secondaries based on pre-generated CORSIKA shower databases.
///        Overrides the larsim CORSIKAGen module to allow to be inherited by a base class.
///
/// \author  gputnam@fnal.gov, taken from module by: Matthew.Bass@physics.ox.ac.uk
////////////////////////////////////////////////////////////////////////

#include "sbncode/EventGenerator/WireModGen/CORSIKAGenSBN.h"

namespace evgen {

  CORSIKAGenSBN::CORSIKAGenSBN(fhicl::ParameterSet const& p)
    : EDProducer{p}
    , fProjectToHeight(p.get<double>("ProjectToHeight", 0.))
    , fShowerInputFiles(p.get<std::vector<std::string>>("ShowerInputFiles"))
    , fShowerCopyType(p.get<std::string>("ShowerCopyType", "IFDH"))
    , fShowerFluxConstants(p.get<std::vector<double>>("ShowerFluxConstants"))
    , fSampleTime(p.get<double>("SampleTime", 0.))
    , fToffset(p.get<double>("TimeOffset", 0.))
    , fBuffBox(p.get<std::vector<double>>("BufferBox", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}))
    , fShowerAreaExtension(p.get<double>("ShowerAreaExtension", 0.))
    , fRandomXZShift(p.get<double>("RandomXZShift", 0.))
    , fGenEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
        createEngine(0, "HepJamesRandom", "gen"),
        "HepJamesRandom",
        "gen",
        p,
        {"Seed", "SeedGenerator"}))
    , fPoisEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
        createEngine(0, "HepJamesRandom", "pois"),
        "HepJamesRandom",
        "pois",
        p,
        "SeedPoisson"))
  {
    if (fShowerInputFiles.size() != fShowerFluxConstants.size() || fShowerInputFiles.size() == 0 ||
        fShowerFluxConstants.size() == 0)
      throw cet::exception("CORSIKAGenSBN")
        << "ShowerInputFiles and ShowerFluxConstants have different or invalid sizes!"
        << "\n";
    fShowerInputs = fShowerInputFiles.size();

    if (fSampleTime == 0.) throw cet::exception("CORSIKAGenSBN") << "SampleTime not set!";

    if (fProjectToHeight == 0.) mf::LogInfo("CORSIKAGenSBN") << "Using 0. for fProjectToHeight!";
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"

    this->openDBs(p.get<std::string>("module_label"));
    this->populateNShowers();
    this->populateTOffset();

    produces<std::vector<simb::MCTruth>>();
    produces<sumdata::RunData, art::InRun>();
  }

  CORSIKAGenSBN::~CORSIKAGenSBN()
  {
    for (int i = 0; i < fShowerInputs; i++) {
      sqlite3_close(fdb[i]);
    }
    //cleanup temp files
    fIFDH->cleanup();
  }

  void CORSIKAGenSBN::ProjectToBoxEdge(const double xyz[],
                                    const double indxyz[],
                                    const double xlo,
                                    const double xhi,
                                    const double ylo,
                                    const double yhi,
                                    const double zlo,
                                    const double zhi,
                                    double xyzout[])
  {

    //we want to project backwards, so take mirror of momentum
    const double dxyz[3] = {-indxyz[0], -indxyz[1], -indxyz[2]};

    // Compute the distances to the x/y/z walls
    double dx = 99.E99;
    double dy = 99.E99;
    double dz = 99.E99;
    if (dxyz[0] > 0.0) { dx = (xhi - xyz[0]) / dxyz[0]; }
    else if (dxyz[0] < 0.0) {
      dx = (xlo - xyz[0]) / dxyz[0];
    }
    if (dxyz[1] > 0.0) { dy = (yhi - xyz[1]) / dxyz[1]; }
    else if (dxyz[1] < 0.0) {
      dy = (ylo - xyz[1]) / dxyz[1];
    }
    if (dxyz[2] > 0.0) { dz = (zhi - xyz[2]) / dxyz[2]; }
    else if (dxyz[2] < 0.0) {
      dz = (zlo - xyz[2]) / dxyz[2];
    }

    // Choose the shortest distance
    double d = 0.0;
    if (dx < dy && dx < dz)
      d = dx;
    else if (dy < dz && dy < dx)
      d = dy;
    else if (dz < dx && dz < dy)
      d = dz;

    // Make the step
    for (int i = 0; i < 3; ++i) {
      xyzout[i] = xyz[i] + dxyz[i] * d;
    }
  }

  void CORSIKAGenSBN::openDBs(std::string const& module_label)
  {
    //choose files based on fShowerInputFiles, copy them with ifdh, open them
    // for c2: statement is unused
    //sqlite3_stmt *statement;
    CLHEP::RandFlat flat(fGenEngine);

    //setup ifdh object
    if (!fIFDH) fIFDH = new ifdh_ns::ifdh;
    const char* ifdh_debug_env = std::getenv("IFDH_DEBUG_LEVEL");
    if (ifdh_debug_env) {
      mf::LogInfo("CORSIKAGenSBN") << "IFDH_DEBUG_LEVEL: " << ifdh_debug_env << "\n";
      fIFDH->set_debug(ifdh_debug_env);
    }

    //get ifdh path for each file in fShowerInputFiles, put into selectedflist
    //if 1 file returned, use that file
    //if >1 file returned, randomly select one file
    //if 0 returned, make exeption for missing files
    std::vector<std::pair<std::string, long>> selectedflist;
    for (int i = 0; i < fShowerInputs; i++) {
      if (fShowerInputFiles[i].find("*") == std::string::npos) {
        //if there are no wildcards, don't call findMatchingFiles
        selectedflist.push_back(std::make_pair(fShowerInputFiles[i], 0));
        mf::LogInfo("CorsikaGen") << "Selected" << selectedflist.back().first << "\n";
      }
      else {
        //use findMatchingFiles
        //default to IFDH approach when wildcards in path
        std::vector<std::pair<std::string, long>> flist;
        std::string path(gSystem->DirName(fShowerInputFiles[i].c_str()));
        std::string pattern(gSystem->BaseName(fShowerInputFiles[i].c_str()));
        flist = fIFDH->findMatchingFiles(path, pattern);
        unsigned int selIndex = -1;
        if (flist.size() == 1) { //0th element is the search path:pattern
          selIndex = 0;
        }
        else if (flist.size() > 1) {
          selIndex = (unsigned int)(flat() * (flist.size() - 1) +
                                    0.5); //rnd with rounding, dont allow picking the 0th element
        }
        else {
          throw cet::exception("CORSIKAGenSBN")
            << "No files returned for path:pattern: " << path << ":" << pattern << std::endl;
        }
        selectedflist.push_back(flist[selIndex]);
        mf::LogInfo("CorsikaGen") << "For " << fShowerInputFiles[i] << ":" << pattern << "\nFound "
                                  << flist.size() << " candidate files"
                                  << "\nChoosing file number " << selIndex << "\n"
                                  << "\nSelected " << selectedflist.back().first << "\n";
      }
    }

    //do the fetching, store local filepaths in locallist
    std::vector<std::string> locallist;
    for (unsigned int i = 0; i < selectedflist.size(); i++) {
      mf::LogInfo("CorsikaGen") << "Fetching: " << selectedflist[i].first << " "
                                << selectedflist[i].second << "\n";
      if (fShowerCopyType == "IFDH") {
        std::string fetchedfile(fIFDH->fetchInput(selectedflist[i].first));
        MF_LOG_DEBUG("CorsikaGen") << "Fetched; local path: " << fetchedfile << "; copy type: IFDH";
        locallist.push_back(fetchedfile);
      }
      else if (fShowerCopyType == "DIRECT") {
        std::string fetchedfile(selectedflist[i].first);
        MF_LOG_DEBUG("CorsikaGen")
          << "Fetched; local path: " << fetchedfile << "; copy type: DIRECT";
        locallist.push_back(fetchedfile);
      }
      else
        throw cet::exception("CORSIKAGenSBN")
          << "Error copying shower db: ShowerCopyType must be \"IFDH\" or \"DIRECT\"\n";
    }

    //open the files in fShowerInputFilesLocalPaths with sqlite3
    for (unsigned int i = 0; i < locallist.size(); i++) {
      //prepare and execute statement to attach db file
      int res = sqlite3_open(locallist[i].c_str(), &fdb[i]);
      if (res != SQLITE_OK)
        throw cet::exception("CORSIKAGenSBN")
          << "Error opening db: (" << locallist[i].c_str() << ") (" << res
          << "): " << sqlite3_errmsg(fdb[i]) << "; memory used:<<" << sqlite3_memory_used() << "/"
          << sqlite3_memory_highwater(0) << "\n";
      else
        mf::LogInfo("CORSIKAGenSBN") << "Attached db " << locallist[i] << "\n";
    }
  }

  double CORSIKAGenSBN::wrapvar(const double var, const double low, const double high)
  {
    //wrap variable so that it's always between low and high
    return (var - (high - low) * floor(var / (high - low))) + low;
  }

  double CORSIKAGenSBN::wrapvarBoxNo(const double var, const double low, const double high, int& boxno)
  {
    //wrap variable so that it's always between low and high
    boxno = int(floor(var / (high - low)));
    return (var - (high - low) * floor(var / (high - low))) + low;
  }

  void CORSIKAGenSBN::populateTOffset()
  {
    //populate TOffset_corsika by finding minimum ParticleTime from db file

    sqlite3_stmt* statement;
    const std::string kStatement("select min(t) from particles");
    double t = 0.;

    for (int i = 0; i < fShowerInputs; i++) {
      //build and do query to get run min(t) from each db
      if (sqlite3_prepare(fdb[i], kStatement.c_str(), -1, &statement, 0) == SQLITE_OK) {
        int res = 0;
        res = sqlite3_step(statement);
        if (res == SQLITE_ROW) {
          t = sqlite3_column_double(statement, 0);
          mf::LogInfo("CORSIKAGenSBN")
            << "For showers input " << i << " found particles.min(t)=" << t << "\n";
          if (i == 0 || t < fToffset_corsika) fToffset_corsika = t;
        }
        else {
          throw cet::exception("CORSIKAGenSBN")
            << "Unexpected sqlite3_step return value: (" << res << "); "
            << "ERROR:" << sqlite3_errmsg(fdb[i]) << "\n";
        }
      }
      else {
        throw cet::exception("CORSIKAGenSBN") << "Error preparing statement: (" << kStatement << "); "
                                           << "ERROR:" << sqlite3_errmsg(fdb[i]) << "\n";
      }
    }

    mf::LogInfo("CORSIKAGenSBN") << "Found corsika timeoffset [ns]: " << fToffset_corsika << "\n";
  }

  void CORSIKAGenSBN::populateNShowers()
  {
    //populate vector of the number of showers per event based on:
    //AREA the showers are being distributed over
    //TIME of the event (fSampleTime)
    //flux constants that determine the overall normalizations (fShowerFluxConstants)
    //Energy range over which the sample was generated (ERANGE_*)
    //power spectrum over which the sample was generated (ESLOPE)

    //compute shower area based on the maximal x,z dimensions of cryostat boundaries + fShowerAreaExtension
    art::ServiceHandle<geo::Geometry const> geom;
    for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {
      double bounds[6] = {0.};
      cryostat.Boundaries(bounds);
      for (unsigned int bnd = 0; bnd < 6; bnd++) {
        mf::LogVerbatim("CORSIKAGenSBN") << "Cryo Boundary: " << bnd << "=" << bounds[bnd]
                                      << " ( + Buffer=" << fBuffBox[bnd] << ")\n";
        if (fabs(bounds[bnd]) > fabs(fShowerBounds[bnd])) { fShowerBounds[bnd] = bounds[bnd]; }
      }
    }
    //add on fShowerAreaExtension without being clever
    fShowerBounds[0] = fShowerBounds[0] - fShowerAreaExtension;
    fShowerBounds[1] = fShowerBounds[1] + fShowerAreaExtension;
    fShowerBounds[4] = fShowerBounds[4] - fShowerAreaExtension;
    fShowerBounds[5] = fShowerBounds[5] + fShowerAreaExtension;

    double showersArea = (fShowerBounds[1] / 100 - fShowerBounds[0] / 100) *
                         (fShowerBounds[5] / 100 - fShowerBounds[4] / 100);

    mf::LogInfo("CORSIKAGenSBN") << "Area extended by : " << fShowerAreaExtension
                              << "\nShowers to be distributed betweeen: x=" << fShowerBounds[0]
                              << "," << fShowerBounds[1] << " & z=" << fShowerBounds[4] << ","
                              << fShowerBounds[5]
                              << "\nShowers to be distributed with random XZ shift: "
                              << fRandomXZShift << " cm"
                              << "\nShowers to be distributed over area: " << showersArea << " m^2"
                              << "\nShowers to be distributed over time: " << fSampleTime << " s"
                              << "\nShowers to be distributed with time offset: " << fToffset
                              << " s"
                              << "\nShowers to be distributed at y: " << fShowerBounds[3] << " cm";

    //db variables
    sqlite3_stmt* statement;
    const std::string kStatement("select erange_high,erange_low,eslope,nshow from input");
    double upperLimitOfEnergyRange = 0., lowerLimitOfEnergyRange = 0., energySlope = 0.,
           oneMinusGamma = 0., EiToOneMinusGamma = 0., EfToOneMinusGamma = 0.;

    for (int i = 0; i < fShowerInputs; i++) {
      //build and do query to get run info from databases
      //  double thisrnd=flat();//need a new random number for each query
      if (sqlite3_prepare(fdb[i], kStatement.c_str(), -1, &statement, 0) == SQLITE_OK) {
        int res = 0;
        res = sqlite3_step(statement);
        if (res == SQLITE_ROW) {
          upperLimitOfEnergyRange = sqlite3_column_double(statement, 0);
          lowerLimitOfEnergyRange = sqlite3_column_double(statement, 1);
          energySlope = sqlite3_column_double(statement, 2);
          fMaxShowers.push_back(sqlite3_column_int(statement, 3));
          oneMinusGamma = 1 + energySlope;
          EiToOneMinusGamma = pow(lowerLimitOfEnergyRange, oneMinusGamma);
          EfToOneMinusGamma = pow(upperLimitOfEnergyRange, oneMinusGamma);
          mf::LogVerbatim("CORSIKAGenSBN")
            << "For showers input " << i << " found e_hi=" << upperLimitOfEnergyRange
            << ", e_lo=" << lowerLimitOfEnergyRange << ", slope=" << energySlope
            << ", k=" << fShowerFluxConstants[i] << "\n";
        }
        else {
          throw cet::exception("CORSIKAGenSBN")
            << "Unexpected sqlite3_step return value: (" << res << "); "
            << "ERROR:" << sqlite3_errmsg(fdb[i]) << "\n";
        }
      }
      else {
        throw cet::exception("CORSIKAGenSBN") << "Error preparing statement: (" << kStatement << "); "
                                           << "ERROR:" << sqlite3_errmsg(fdb[i]) << "\n";
      }

      //this is computed, how?
      double NShowers = (M_PI * showersArea * fShowerFluxConstants[i] *
                         (EfToOneMinusGamma - EiToOneMinusGamma) / oneMinusGamma) *
                        fSampleTime;
      fNShowersPerEvent.push_back(NShowers);
      mf::LogVerbatim("CORSIKAGenSBN")
        << "For showers input " << i << " the number of showers per event is " << (int)NShowers
        << "\n";
    }
  }

  void CORSIKAGenSBN::GetSample(simb::MCTruth& mctruth)
  {
    //for each input, randomly pull fNShowersPerEvent[i] showers from the Particles table
    //and randomly place them in time (between -fSampleTime/2 and fSampleTime/2)
    //wrap their positions based on the size of the area under consideration
    //based on http://nusoft.fnal.gov/larsoft/doxsvn/html/CRYHelper_8cxx_source.html (Sample)

    //query from sqlite db with select * from particles where shower in (select id from showers ORDER BY substr(id*0.51123124141,length(id)+2) limit 100000) ORDER BY substr(shower*0.51123124141,length(shower)+2);
    //where 0.51123124141 is a random seed to allow randomly selecting rows and should be randomly generated for each query
    //the inner order by is to select randomly from the possible shower id's
    //the outer order by is to make sure the shower numbers are ordered randomly (without this, the showers always come out ordered by shower number
    //and 100000 is the number of showers to be selected at random and needs to be less than the number of showers in the showers table

    //TDatabasePDG is for looking up particle masses
    static TDatabasePDG* pdgt = TDatabasePDG::Instance();

    //db variables
    sqlite3_stmt* statement;
    const TString kStatement(
      "select shower,pdg,px,py,pz,x,z,t,e from particles where shower in (select id from showers "
      "ORDER BY substr(id*%f,length(id)+2) limit %d) ORDER BY substr(shower*%f,length(shower)+2)");

    CLHEP::RandFlat flat(fGenEngine);
    CLHEP::RandPoissonQ randpois(fPoisEngine);

    // get geometry and figure where to project particles to, based on CRYHelper
    art::ServiceHandle<geo::Geometry const> geom;
    double x1, x2;
    double y1, y2;
    double z1, z2;
    geom->WorldBox(&x1, &x2, &y1, &y2, &z1, &z2);

    // make the world box slightly smaller so that the projection to
    // the edge avoids possible rounding errors later on with Geant4
    double fBoxDelta = 1.e-5;
    x1 += fBoxDelta;
    x2 -= fBoxDelta;
    y1 += fBoxDelta;
    y2 = fProjectToHeight;
    z1 += fBoxDelta;
    z2 -= fBoxDelta;

    //populate mctruth
    int ntotalCtr = 0; //count number of particles added to mctruth
    int lastShower =
      0; //keep track of last shower id so that t can be randomized on every new shower
    int nShowerCntr = 0; //keep track of how many showers are left to be added to mctruth
    int nShowerQry = 0;  //number of showers to query from db
    int shower, pdg;
    double px, py, pz, x, z, tParticleTime, etot,
      showerTime = 0., showerTimex = 0., showerTimez = 0., showerXOffset = 0., showerZOffset = 0.,
      t;
    for (int i = 0; i < fShowerInputs; i++) {
      nShowerCntr = randpois.fire(fNShowersPerEvent[i]);
      mf::LogInfo("CORSIKAGEN") << " Shower input " << i << " with mean " << fNShowersPerEvent[i]
                                << " generating " << nShowerCntr;

      while (nShowerCntr > 0) {
        //how many showers should we query?
        if (nShowerCntr > fMaxShowers[i]) {
          nShowerQry = fMaxShowers[i]; //take the group size
        }
        else {
          nShowerQry = nShowerCntr; //take the rest that are needed
        }
        //build and do query to get nshowers
        double thisrnd = flat(); //need a new random number for each query
        TString kthisStatement = TString::Format(kStatement.Data(), thisrnd, nShowerQry, thisrnd);
        MF_LOG_DEBUG("CORSIKAGenSBN") << "Executing: " << kthisStatement;
        if (sqlite3_prepare(fdb[i], kthisStatement.Data(), -1, &statement, 0) == SQLITE_OK) {
          int res = 0;
          //loop over database rows, pushing particles into mctruth object
          while (1) {
            res = sqlite3_step(statement);
            if (res == SQLITE_ROW) {
              /*
               * Memo columns:
               * [0] shower
               * [1] particle ID (PDG)
               * [2] momentum: x component [GeV/c]
               * [3] momentum: y component [GeV/c]
               * [4] momentum: z component [GeV/c]
               * [5] position: x component [cm]
               * [6] position: z component [cm]
               * [7] time [ns]
               * [8] energy [GeV]
               */
              shower = sqlite3_column_int(statement, 0);
              if (shower != lastShower) {
                //each new shower gets its own random time and position offsets
                showerTime = 1e9 * (flat() * fSampleTime);  //converting from s to ns
                showerTimex = 1e9 * (flat() * fSampleTime); //converting from s to ns
                showerTimez = 1e9 * (flat() * fSampleTime); //converting from s to ns
                //and a random offset in both z and x controlled by the fRandomXZShift parameter
                showerXOffset = flat() * fRandomXZShift - (fRandomXZShift / 2);
                showerZOffset = flat() * fRandomXZShift - (fRandomXZShift / 2);
              }
              pdg = sqlite3_column_int(statement, 1);
              //get mass for this particle
              double m = 0.; // in GeV
              TParticlePDG* pdgp = pdgt->GetParticle(pdg);
              if (pdgp) m = pdgp->Mass();

              //Note: position/momentum in db have north=-x and west=+z, rotate so that +z is north and +x is west
              //get momentum components
              px = sqlite3_column_double(statement, 4); //uboone x=Particlez
              py = sqlite3_column_double(statement, 3);
              pz = -sqlite3_column_double(statement, 2); //uboone z=-Particlex
              etot = sqlite3_column_double(statement, 8);

              //get/calculate position components
              int boxnoX = 0, boxnoZ = 0;
              x = wrapvarBoxNo(sqlite3_column_double(statement, 6) + showerXOffset,
                               fShowerBounds[0],
                               fShowerBounds[1],
                               boxnoX);
              z = wrapvarBoxNo(-sqlite3_column_double(statement, 5) + showerZOffset,
                               fShowerBounds[4],
                               fShowerBounds[5],
                               boxnoZ);
              tParticleTime = sqlite3_column_double(
                statement, 7); //time offset, includes propagation time from top of atmosphere
              //actual particle time is particle surface arrival time
              //+ shower start time
              //+ global offset (fcl parameter, in s)
              //- propagation time through atmosphere
              //+ boxNo{X,Z} time offset to make grid boxes have different shower times
              t = tParticleTime + showerTime + (1e9 * fToffset) - fToffset_corsika +
                  showerTimex * boxnoX + showerTimez * boxnoZ;
              //wrap surface arrival so that it's in the desired time window
              t = wrapvar(t, (1e9 * fToffset), 1e9 * (fToffset + fSampleTime));

              simb::MCParticle p(ntotalCtr, pdg, "primary", -200, m, 1);

              //project back to wordvol/fProjectToHeight
              /*
               * This back propagation goes from a point on the upper surface of
               * the cryostat back to the edge of the world, except that that
               * world is cut short by `fProjectToHeight` (`y2`) ceiling.
               * The projection will most often lie on that ceiling, but it may
               * end up instead on one of the side edges of the world, or even
               * outside it.
               */
              double xyzo[3];
              double x0[3] = {x, fShowerBounds[3], z};
              double dx[3] = {px, py, pz};
              this->ProjectToBoxEdge(x0, dx, x1, x2, y1, y2, z1, z2, xyzo);

              TLorentzVector pos(
                xyzo[0], xyzo[1], xyzo[2], t); // time needs to be in ns to match GENIE, etc
              TLorentzVector mom(px, py, pz, etot);
              p.AddTrajectoryPoint(pos, mom);
              mctruth.Add(p);
              ntotalCtr++;
              lastShower = shower;
            }
            else if (res == SQLITE_DONE) {
              break;
            }
            else {
              throw cet::exception("CORSIKAGenSBN")
                << "Unexpected sqlite3_step return value: (" << res << "); "
                << "ERROR:" << sqlite3_errmsg(fdb[i]) << "\n";
            }
          }
        }
        else {
          throw cet::exception("CORSIKAGenSBN")
            << "Error preparing statement: (" << kthisStatement << "); "
            << "ERROR:" << sqlite3_errmsg(fdb[i]) << "\n";
        }
        nShowerCntr = nShowerCntr - nShowerQry;
      }
    }
  }

  void CORSIKAGenSBN::beginRun(art::Run& run)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()), art::fullRun());
  }

  void CORSIKAGenSBN::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<simb::MCTruth>> truthcol(new std::vector<simb::MCTruth>);

    art::ServiceHandle<geo::Geometry const> geom;

    simb::MCTruth truth;
    truth.SetOrigin(simb::kCosmicRay);

    simb::MCTruth pretruth;
    GetSample(pretruth);
    mf::LogInfo("CORSIKAGenSBN") << "GetSample number of particles returned: " << pretruth.NParticles()
                              << "\n";
    // loop over particles in the truth object
    for (int i = 0; i < pretruth.NParticles(); ++i) {
      simb::MCParticle particle = pretruth.GetParticle(i);
      const TLorentzVector& v4 = particle.Position();
      const TLorentzVector& p4 = particle.Momentum();
      double x0[3] = {v4.X(), v4.Y(), v4.Z()};
      double dx[3] = {p4.Px(), p4.Py(), p4.Pz()};

      // now check if the particle goes through any cryostat in the detector
      // if so, add it to the truth object.
      for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {
        double bounds[6] = {0.};
        cryostat.Boundaries(bounds);

        //add a buffer box around the cryostat bounds to increase the acceptance and account for scattering
        //By default, the buffer box has zero size
        for (unsigned int cb = 0; cb < 6; cb++)
          bounds[cb] = bounds[cb] + fBuffBox[cb];

        //calculate the intersection point with each cryostat surface
        bool intersects_cryo = false;
        for (int bnd = 0; bnd != 6; ++bnd) {
          if (bnd < 2) {
            double p2[3] = {bounds[bnd],
                            x0[1] + (dx[1] / dx[0]) * (bounds[bnd] - x0[0]),
                            x0[2] + (dx[2] / dx[0]) * (bounds[bnd] - x0[0])};
            if (p2[1] >= bounds[2] && p2[1] <= bounds[3] && p2[2] >= bounds[4] &&
                p2[2] <= bounds[5]) {
              intersects_cryo = true;
              break;
            }
          }
          else if (bnd >= 2 && bnd < 4) {
            double p2[3] = {x0[0] + (dx[0] / dx[1]) * (bounds[bnd] - x0[1]),
                            bounds[bnd],
                            x0[2] + (dx[2] / dx[1]) * (bounds[bnd] - x0[1])};
            if (p2[0] >= bounds[0] && p2[0] <= bounds[1] && p2[2] >= bounds[4] &&
                p2[2] <= bounds[5]) {
              intersects_cryo = true;
              break;
            }
          }
          else if (bnd >= 4) {
            double p2[3] = {x0[0] + (dx[0] / dx[2]) * (bounds[bnd] - x0[2]),
                            x0[1] + (dx[1] / dx[2]) * (bounds[bnd] - x0[2]),
                            bounds[bnd]};
            if (p2[0] >= bounds[0] && p2[0] <= bounds[1] && p2[1] >= bounds[2] &&
                p2[1] <= bounds[3]) {
              intersects_cryo = true;
              break;
            }
          }
        }

        if (intersects_cryo) {
          truth.Add(particle);
          break; //leave loop over cryostats to avoid adding particle multiple times
        }        // end if particle goes into a cryostat
      }          // end loop over cryostats in the detector

    } // loop on particles

    mf::LogInfo("CORSIKAGenSBN")
      << "Number of particles from getsample crossing cryostat + bounding box: "
      << truth.NParticles() << "\n";

    truthcol->push_back(truth);
    evt.put(std::move(truthcol));

    return;
  } // end produce

} // end namespace
