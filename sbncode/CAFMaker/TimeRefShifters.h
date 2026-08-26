/**
 * @file   sbncode/CAFMaker/TimeRefShifters.h
 * @brief  Simple utilities to apply a reference time shift.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 21, 2026
 * 
 * The purpose of these utilities is to ensure that the time reference shifts
 * are applied in a consistent way.
 */

#ifndef SBNCODE_CAFMAKER_TIMEREFSHIFTERS_H
#define SBNCODE_CAFMAKER_TIMEREFSHIFTERS_H


// C++ standard libraries
#include <cstdint> // std::int64_t


// -----------------------------------------------------------------------------
namespace caf {
  template <typename RefTime = double> struct TimeRefShifter;
  struct CRTtimeRefShifter;
}
/**
 * @brief Applies a standard time reference change to a time.
 * @tparam RefTime type used to internally store the reference time shift
 * 
 * The shifter object is constructed with a time shift value: the reference time
 * is delayed in time by that value. As a consequence, all time points will
 * acquire an "earlier" value.
 * 
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * double T1 = 0.5, T2 = 1.0; // [us]
 * int T3 = 50; // [ns]
 * caf::TimeRefShifter const shifter{ 0.235 };
 * double shiftedT1 = shifter.shiftedTime(T1);
 * shifter.shift(T2);
 * int shiftedT3 = shifter.shiftedTimeNS(T3);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will shift the time reference forward by 235 ns; the value of all shifted
 * times will become 235 ns smaller (`0.265` for `shiftedT1`, `0.765` for `T2`,
 * `-185` for `T3`; `T1` was not shifted).
 * 
 */
template <typename RefTime>
struct caf::TimeRefShifter {
  
  /// Forward shift of the reference time [us]
  RefTime const refTimeShift = RefTime{ 0 };
  /// Returns a time with reference shifted compared to `t` [us]
  
  /// Returns `t` [us], shifted.
  template <typename T>
  constexpr T shiftedTime(T t) const { return t - shiftAs<T>(); }
  
  /// Returns `t` [us], shifted only if `doShift` is `true`.
  template <typename T>
  constexpr T shiftedTimeIf(T t, bool doShift) const
    { return doShift? shiftedTime(t): t; }
  
  /// Returns `t` [us], shifted only if it's not exactly equal to `defVal`.
  template <typename T>
  constexpr T shiftedTimeIfNot(T t, T defVal) const
    { return shiftedTimeIf(t, t != defVal); }
  
  /// Returns a time with reference shifted compared to `t` [ns]
  template <typename T>
  constexpr T shiftedTimeNS(T t) const { return t - shiftNSas<T>(); }
  
  /// Changes the value of `t` shifting its reference [us]
  template <typename T>
  T& shift(T& t) const { return t -= shiftAs<T>(); }
  
  /// Changes the value of `t` shifting its reference [ns]
  template <typename T>
  T& shiftNS(T& t) const { return t -= shiftNSas<T>(); }
  
  
    private:
  
  template <typename T>
  constexpr T shiftAs() const { return static_cast<T>(refTimeShift); }
  template <typename T>
  constexpr T shiftNSas() const { return static_cast<T>(refTimeShift * 1000); }
  
}; // caf::TimeRefShifter

// template deduction guide: if presented with an argument T...
namespace caf { template <typename T> TimeRefShifter(T) -> TimeRefShifter<T>; }

// -----------------------------------------------------------------------------
/**
 * @brief Applies time reference shifts to CRT times and timestamps.
 * 
 * This shifter works similarly to `caf::TimeRefShifter` but it handles relative
 * times and absolute timestamps independent shifts.
 * 
 * This object largely adopts the terminology from `sbn::crt::CRTHit`, calling
 * time stamps with "TS0" and relative times as "TS1".
 * In addition, it supports a "merged" time which is either from TS0 or TS1
 * depending on a fixes parameter.
 * 
 * Note that the convention of `sbn::crt::CRTHit` are:
 *  * timestamps are expressed as signed 64-bit integers;
 *  * relative times are expressed as `double`.
 * 
 * They are both expressed in nanoseconds. However, across the CRT-related data
 * products conventions do change.
 * 
 * The output follows CAF conventions instead, so it is always in `double`
 * precision and in microseconds. That may not play well with actual
 * 64-bit timestamps (TS0-style), which may lose precision.
 */
struct caf::CRTtimeRefShifter {

  bool const useTS0 = false; ///< Use TS0 for time?
  std::int64_t const refTimeShiftTS = 0; ///< Reference shift for timestamps [ns]
  double const refTimeShift = 0.0; ///< Reference shifts for relative times [us]
  
  
  constexpr bool usingTS0() const { return useTS0; }
  
  /// Returns the timestamp `ts_ns` with shifted reference [us]
  constexpr double shiftTimestamp(std::int64_t ts_ns) const
    { return (ts_ns - refTimeShiftTS) / 1000.; }
    
  /// Returns the time `t` [ns] with shifted reference [us]
  constexpr double shiftTime(double t_ns) const
    { return t_ns / 1000. - refTimeShift; }
  
  /// Definitions here follow `sbn::crt::CRTHit` (i.e. TS1 is not really a timestamp)
  template <typename T>
  void fillCRTtimes(
    std::int64_t TS0, std::int64_t TS1,
    T* t0 = nullptr, T* t1 = nullptr, T* time = nullptr
    ) const
    {
      double const shT0 = shiftTimestamp(TS0);
      double const shT1 = shiftTime(static_cast<double>(TS1));
      if (t0) *t0 = shT0;
      if (t1) *t1 = shT1;
      if (time) *time = useTS0? shT0: shT1;
    }
  
  
  /// Creates a timestamp (of type `TS`) from a second and a nanosecond part.
  template <typename S, typename NS, typename TS = std::int64_t>
  static constexpr TS makeTS(S s, NS ns)
    {
      return static_cast<TS>
        (static_cast<std::int64_t>(s) * 1'000'000'000 + static_cast<std::int64_t>(ns));
    }
  
}; // caf::CRTtimeRefShifter


// -----------------------------------------------------------------------------


#endif // SBNCODE_CAFMAKER_TIMEREFSHIFTERS_H
