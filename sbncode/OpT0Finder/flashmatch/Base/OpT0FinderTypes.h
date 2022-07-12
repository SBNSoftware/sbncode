#ifndef OPT0FINDER_OPT0FINDERTYPES_H
#define OPT0FINDER_OPT0FINDERTYPES_H

#include <vector>
#include <numeric>
#include "OpT0FinderConstants.h"
#include <string>
#include <cmath>

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/GeoAlgo/GeoVector.h"
#else
#include "sbncode/OpT0Finder/flashmatch/GeoAlgo/GeoVector.h"
#endif

namespace flashmatch {

  /// Index used to identify Flash_t/QPointCollection_t uniquely in an event
  typedef size_t ID_t;
  /// Invalid ID
  const ID_t kINVALID_ID = kINVALID_SIZE;

  /// Enumerator for different types of algorithm
  enum Algorithm_t {
    kTPCFilter,       ///< Algorithm type to filter out TPC objects from matching candidate list
    kFlashFilter,     ///< Algorithm type to filter out flash from matching candidate list
    kFlashMatch,      ///< Algorithm type to match flash hypothesis and reconstructed flash
    kMatchProhibit,   ///< Algorithm type to prohibit a match between a flash and a cluster
    kFlashHypothesis, ///< Algorithm type to make QCluster_t => Flash_t hypothesis
    kCustomAlgo,      ///< Algorithm type that does not play a role in the framework execution but inherits from BaseAlgorithm
    kTouchMatch,      ///< Algorithm type to match flash time and TPC object based on geometry
    kMatchSelect,     ///< Algorithm type to select the final set of matched pairs
    kAlgorithmTypeMax ///< enum flag for algorithm type count & invalid type
  };

  /// Struct to represent an optical flash
  struct Flash_t {
  public:

    std::vector<double> pe_v; ///< PE distribution over photo-detectors
    std::vector<double> pe_true_v; ///< PE distribution over photo-detectors of MCFlash
    std::vector<double> pe_err_v; ///< PE value error
    double x,y,z;             ///< Flash position
    double x_err,y_err,z_err; ///< Flash position error
    double time;              ///< Flash timing, a candidate T0
    double time_true;         ///< MCFlash timing (if it was matched to a MCFlash)
    double time_width;        ///< Flash time integration window
    double dt_next, dt_prev;
    ID_t idx;                 ///< index from original larlite vector
    //ID_t ROOT_idx;            ///< index in root file
    /// Default ctor assigns invalid values
    Flash_t() : pe_v(), pe_true_v() , pe_err_v() {
      x = y = z = kINVALID_DOUBLE;
      x_err = y_err = z_err = kINVALID_DOUBLE;
      time = kINVALID_DOUBLE;
      time_true = kINVALID_DOUBLE;
      time_width = kINVALID_DOUBLE;
      dt_next = kINVALID_DOUBLE;
      dt_prev = kINVALID_DOUBLE;
      idx = kINVALID_ID;
      //ROOT_idx = kINVALID_ID;
    }

    /// maximum x
    inline double max_pe() const
    { double x=-flashmatch::kINVALID_DOUBLE; for(auto const& v : pe_v) x = std::max(x, v); return x; }

    /// maximum x in pe_true_v
    inline double max_pe_true() const
    { double x=-flashmatch::kINVALID_DOUBLE; for(auto const& v : pe_true_v) x = std::max(x, v); return x; }

    /// minimum x
    inline double min_pe() const
    { double x=flashmatch::kINVALID_DOUBLE; for(auto const& v : pe_v) x = std::min(x, v); return x; }

    /// minimum x in pe_true_v
    inline double min_pe_true() const
    { double x=flashmatch::kINVALID_DOUBLE; for(auto const& v : pe_true_v) x = std::min(x, v); return x; }

    /// Total PE calculation
    double TotalPE() const {
      double res=0.;
      for(auto const& v : pe_v) if(v>=0.) res+=v;
      return res;
    }
    /// Total true PE calculation
    double TotalTruePE() const {
        double res=0.;
        for (auto const& v : pe_true_v) if (v>=0.) res+=v;
        return res;
    }
    /// Check validity
    bool Valid(size_t nopdet=0) const {
      return (nopdet ? (pe_v.size() == nopdet && pe_err_v.size() == nopdet) : (pe_v.size() == pe_err_v.size()));
    }
    //double TotalPE() const{ return std::accumulate(pe_v.begin(),pe_v.end(),0.0);}
  };

  /// Struct to represent an energy deposition point in 3D space
  struct QPoint_t{

    double x,y,z; ///< Spatial position in [cm]
    double q;     ///< Charge in an arbitrary unit
    /// Default ctor assigns invalid values
    QPoint_t()
      : x(kINVALID_DOUBLE)
      , y(kINVALID_DOUBLE)
      , z(kINVALID_DOUBLE)
      , q(kINVALID_DOUBLE)
    {}
    /// Alternative ctor
    QPoint_t(double xvalue,
	     double yvalue,
	     double zvalue,
	     double qvalue)
      : x(xvalue)
      , y(yvalue)
      , z(zvalue)
      , q(qvalue)
    {}
    /// distance
    inline double dist(const QPoint_t& pt) const
    { return sqrt(pow(pt.x-x,2)+pow(pt.y-y,2)+pow(pt.z-z,2)); }
  };

  /// Collection of charge deposition 3D point (cluster)
  class QCluster_t : public std::vector<QPoint_t>{
  public:
    ID_t idx;     ///< index from original container
    //ID_t ROOT_idx;///< index in original root file
    double time;  ///< assumed time w.r.t. trigger for reconstruction
    double time_true; ///< Time from MCTrack information
    double min_x_true; ///< True x-minimum value

    /// Default constructor
    QCluster_t() : idx(kINVALID_ID), time(kINVALID_DOUBLE), time_true(kINVALID_DOUBLE), min_x_true(kINVALID_DOUBLE) {}
    ~QCluster_t() {}

    /// returns the sum of "q" from QPoint_t
    double sum() const;

    /// returns the total trajectory length
    double length() const;

    /// drop points outside the x range specified
    void drop(double xmin, double xmax);

    /// drop points outside xyz range specified
    void drop(double xmin, double ymin, double zmin,
              double xmax, double ymax, double zmax);

    /// minimum x point
    void min_x(flashmatch::QPoint_t& qpt) const;

    /// minimum x point
    void min_x(geoalgo::Point_t& pt) const;

    /// maximum x point
    void max_x(flashmatch::QPoint_t& qpt) const;

    /// maximum x point
    void max_x(geoalgo::Point_t& pt) const;

    /// minimum x
    inline double min_x() const
    { double x=flashmatch::kINVALID_DOUBLE; for(auto const& pt : (*this)) x = std::min(x,pt.x); return x; }

    /// maximum x
    inline double max_x() const
    { double x=-flashmatch::kINVALID_DOUBLE; for(auto const& pt : (*this)) x = std::max(x,pt.x); return x; }

    inline QCluster_t& operator+=(const double shift)
    { for(auto& pt : (*this)) pt.x += shift; return (*this); }

    inline QCluster_t operator+(const double shift)
    { auto result = (*this); result += shift; return result; }

    inline QCluster_t& operator-=(const double shift)
    { for(auto& pt : (*this)) pt.x -= shift; return (*this); }

    inline QCluster_t operator-(const double shift)
    { auto result = (*this); result -= shift; return result; }

    inline QCluster_t& operator+=(const QCluster_t& rhs) {
      this->reserve(rhs.size() + this->size());
      for(auto const& pt : rhs) this->push_back(pt);
      return (*this);
    }

    inline QCluster_t operator+(const QCluster_t& rhs) const {
      QCluster_t res((*this));
      res += rhs;
      return res;
    }

  };
  std::ostream& operator << (std::ostream& out, const flashmatch::QCluster_t& obj);

  /// Collection of 3D point clusters (one use case is TPC object representation for track(s) and shower(s))
  typedef std::vector<flashmatch::QCluster_t> QClusterArray_t;
  /// Collection of Flash objects
  typedef std::vector<flashmatch::Flash_t> FlashArray_t;

  /// Index collection
  typedef std::vector<flashmatch::ID_t> IDArray_t;

  /// Enum to define a type of touch match
  enum TouchMatch_t {
    kNoTouchMatch,
    kAnodeCrossing,
    kCathodeCrossing,
    kAnodeCathodeCrossing
  };

  /// Flash-TPC match info
  struct FlashMatch_t {
    ID_t tpc_id;   ///< matched TPC object ID
    ID_t flash_id; ///< matched Flash ID
    std::vector<double> hypothesis; ///< Hypothesis flash object
    double score;             ///< floating point representing the "goodness" (algorithm dependent), should be 0.0=>1.0
    QPoint_t tpc_point;       ///< estimated & matched 3D flash hypothesis point from TPC information
    QPoint_t tpc_point_err;   ///< error on the estimated point
    TouchMatch_t touch_match; ///< see TouchMatch_t enum
    double touch_score;       ///< score for touch match, should be 0.0=>1.0
    QPoint_t touch_point;     ///< estimated & matched 3D point based on BaseTouchMatch inherited algorithm

    // FIXME: attributes below are meant for a particular algorithm QLLMatch, not meant to be here for long... 
	  unsigned int duration;  ///< Computation time of the match algorithm on this match (ns)
    unsigned int num_steps; ///< Number of MIGRAD steps
    double minimizer_min_x; ///< the minimum X value MIGRAD tried out
    double minimizer_max_x; ///< the maximum X value MIGRAD tried out

    /// Default ctor assigns invalid values
    FlashMatch_t() : tpc_id(kINVALID_ID), flash_id(kINVALID_ID), hypothesis(),
    score(-1), touch_match(kNoTouchMatch), touch_score(-1), duration(0)
    {}

  };

  /// Enum to define MC source type (MCTrack or MCShower) for a given QCluster
  enum MCAncestor_t {
    kMCShowerAncestor,
    kMCTrackAncestor,
    kUnknownAncestor
  };

  /// Struct to represent the ancestor information for a specific interaction (QCluster)
  struct MCSource_t {
    int     index_id; ///< MCTrac/MCShower collection index ID of the ancestor
    double  g4_time;  ///< Interaction G4 time in micro-seconds
    double  energy_deposit;   ///< Deposited energy total
    MCAncestor_t source_type; ///< Ancestor source type
    MCSource_t()
    {
      index_id = -1;
      g4_time  = kINVALID_DOUBLE;
      source_type = kUnknownAncestor;
    }
  };

  namespace msg {
    /// Verbosity message level
    enum Level_t {
      kDEBUG,
      kINFO,
      kNORMAL,
      kWARNING,
      kERROR,
      kCRITICAL,
      kMSG_TYPE_MAX
    };

    const std::string kStringPrefix[kMSG_TYPE_MAX] =
      {
	"\033[94m     [DEBUG]  \033[00m", ///< DEBUG message prefix
	"\033[92m      [INFO]  \033[00m", ///< INFO message prefix
	"\033[95m    [NORMAL]  \033[00m", ///< NORMAL message prefix
	"\033[93m   [WARNING]  \033[00m", ///< WARNING message prefix
	"\033[91m     [ERROR]  \033[00m", ///< ERROR message prefix
	"\033[5;1;33;41m [EXCEPTION]  \033[00m"  ///< CRITICAL message prefix
      };
    ///< Prefix of message
  }
}
#endif
