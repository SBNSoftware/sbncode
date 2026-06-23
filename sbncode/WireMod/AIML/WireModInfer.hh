/**
 * WireModInfer.hh - WaireMLod interface for inference at SBN
 *
 * Loads weights for the Julia WaireMLod network from binary file
 * and runs model inference on hits to modify them to be more ``data-like``
 *
 * ── Caller-facing input (11 primitive features, 0-indexed) ──────────────────
 *
 * The external ABI keeps the legacy 11-feature primitive layout so the SBN
 * caller code (`produceNew` in the .cc) does not have to change. Feature
 * augmentation (sin/cos θXW, sin/cos 2·θ_plane, block_is_odd) happens
 * *inside* this header (see `derive_features`), keeping all wire-geometry
 * knowledge embedded in the inference library.
 *
 *   Index  Symbol      Unit       Description
 *   0      x           cm         Drift position
 *   1      y           cm         Vertical position
 *   2      z           cm         Beam-direction position
 *   3      dir_x       —          Track unit vector x-component
 *   4      dir_y       —          Track unit vector y-component
 *   5      dir_z       —          Track unit vector z-component
 *   6      dir_y_rel   —          Plane-relative dir_y
 *   7      dir_z_rel   —          Plane-relative dir_z
 *   8      channel     —          Global wire channel number [0, 55295]
 *   9      theta_xw    degrees    Wire-crossing angle (atan(dir_x/dir_z))
 *   10     dqdx        ADC/cm     Charge deposition per unit path length
 *
 * ── Network-facing input (15 production features, 0-indexed) ────────────────
 *
 * Built by `derive_features` from the 11-primitive caller input. Must match
 * the row order in src/DataPipeline.jl (the Julia training pipeline).
 *
 *   Index  Symbol         Unit      Source
 *   0      x              cm        raw[0]
 *   1      y              cm        raw[1]
 *   2      z              cm        raw[2]
 *   3      dir_x          —         raw[3]
 *   4      dir_y          —         raw[4]
 *   5      dir_z          —         raw[5]
 *   6      dir_y_rel      —         raw[6]
 *   7      dir_z_rel      —         raw[7]
 *   8      channel        —         raw[8]
 *   9      sin_θXW        —         sin(raw[9] · π/180)
 *   10     cos_θXW        —         cos(raw[9] · π/180)
 *   11     sin_2θ_plane   —         sin(2 · plane_theta(raw[8]))
 *   12     cos_2θ_plane   —         cos(2 · plane_theta(raw[8]))
 *   13     block_is_odd   {0,1}     ((int)raw[8] / 13824) is odd ? 1 : 0
 *   14     dqdx           ADC/cm    raw[10]
 *
 * ── Weight-file format version ──────────────────────────────────────────────
 *
 * Accepts format versions 2, 3, 4, 5, and 6 (all headers carry `n_in`).
 *
 *   v2 — single-routing model. One Dense(128→4)+softmax `cls` head produces the
 *        class routing weights; the MoE runs once; output heads use softplus
 *        and the (Â, ŵ) → (Î, ŵ) combine.
 *   v3 — split-routing model. TWO independent Dense(128→4)+softmax routers
 *        (`routing_a`, `routing_w`) replace the single `cls` head. The backbone
 *        (InputLayer → Dense → ResBlock1 → ResBlock2) emits 128-d features; each
 *        router prepends its softmax to those features and the MoE runs once per
 *        router. The integral Î is taken from the routing_a run and the width ŵ
 *        from the routing_w run. The output heads use identity activation, and
 *        the combine is log-space:
 *           Â = exp(clamp(ahead, -50, 30));  ŵ = exp(clamp(whead, -50, 30))
 *           Î = √(2π)·ŵ·Â;  ŵ = ŵ
 *   v4 — same as v3, but the output heads are 2-layer MLPs:
 *           Dense(132 → HEAD_HIDDEN_DIM, gelu) → Dense(HEAD_HIDDEN_DIM → 1, id)
 *        Tensors are exposed as `ahead.l1.*`, `ahead.l2.*`, `whead.l1.*`,
 *        `whead.l2.*` instead of single `ahead.*` / `whead.*`. Combine
 *        semantics (exp + clamp + √(2π)) are identical to v3.
 *   v5 — asymmetric "lean" model. ONE router + a 3-CWL A-expert bank → logÂ
 *        (2-layer log head `ahead.l1/l2`), and a 3-layer DENSE w-head
 *        (`whead.l1/l2/l3`) → logŵ. Expert count K is read at runtime from
 *        `routing_a` (not the compile-time N_CLASSES). Combine as v3/v4.
 *   v6 — decoupled model. Shared backbone + TWO independent 3-CWL expert banks:
 *        routing_a→bank_A→logÂ (`cwl_a1/2/3`, `ahead.l1/l2`) and
 *        routing_w→bank_W→logŵ (`cwl_w1/2/3`, `whead.l1/l2`), same K for both.
 *        Single pass (each bank run once); combine as v3/v4 (exp+clamp+√(2π)).
 *
 * Version-5/6 readers auto-detect K from `routing_a.bias`, so 4- or 9-expert
 * files both load.
 *
 * Version-1 weight files were produced before the 15-feature augmentation
 * landed and are rejected with a clear error; retrain with the current
 * src/Training.jl to upgrade.
 */

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

namespace sys
{
  // ─────────────────────────────────────────────────────────────────────────────
  // WMHit Class
  // ─────────────────────────────────────────────────────────────────────────────
  class WMHit : public recob::Hit
  {
    public:
      WMHit(raw::ChannelID_t channel,
            raw::TDCtick_t   start_tick,
            raw::TDCtick_t   end_tick,
            float            peak_time,
            float            sigma_peak_time,
            float            rms,
            float            peak_amplitude,
            float            sigma_peak_amplitude,
            float            ROISummedADC,
            float            HitSummedADC,
            float            hit_integral,
            float            hit_sigma_integral,
            short int        multiplicity,
            short int        local_index,
            float            goodness_of_fit,
            int              dof,
            geo::View_t      view,
            geo::SigType_t   signal_type,
            geo::WireID      wireID,
            float            integralRatio,
            float            widthRatio) : recob::Hit(channel,
                                                      start_tick,
                                                      end_tick,
                                                      peak_time,
                                                      sigma_peak_time,
                                                      rms,
                                                      peak_amplitude,
                                                      sigma_peak_amplitude,
                                                      ROISummedADC,
                                                      HitSummedADC,
                                                      hit_integral,
                                                      hit_sigma_integral,
                                                      multiplicity,
                                                      local_index,
                                                      goodness_of_fit,
                                                      dof,
                                                      view,
                                                      signal_type,
                                                      wireID),
                                           fIntegralRatio(integralRatio),
                                           fWidthRatio(widthRatio)
                                           {}
      WMHit(recob::Hit hit, float integralRatio = 1.0, float widthRatio = 1.0)
        : recob::Hit(hit.Channel(),
                     hit.StartTick(),
                     hit.EndTick(),
                     hit.PeakTime(),
                     hit.SigmaPeakTime(),
                     hit.RMS() * widthRatio,
                     hit.PeakAmplitude() * integralRatio / widthRatio,
                     hit.SigmaPeakAmplitude() * integralRatio / widthRatio,
                     hit.ROISummedADC() * integralRatio,
                     hit.HitSummedADC() * integralRatio,
                     hit.Integral() * integralRatio,
                     hit.SigmaIntegral() * integralRatio,
                     hit.Multiplicity(),
                     hit.LocalIndex(),
                     hit.GoodnessOfFit(),
                     hit.DegreesOfFreedom(),
                     hit.View(),
                     hit.SignalType(),
                     hit.WireID()),
          fIntegralRatio(integralRatio),
          fWidthRatio(widthRatio)
          {}
      float IntegralRatio() { return fIntegralRatio; }
      float WidthRatio() { return fWidthRatio; }
    private:
      float fIntegralRatio;
      float fWidthRatio;
  };

  // ─────────────────────────────────────────────────────────────────────────────
  // Model constants  (must match Constants.jl + DataPipeline.jl row order)
  // ─────────────────────────────────────────────────────────────────────────────
  //
  // Two distinct widths matter here:
  //   • N_RAW       — caller-facing primitive input (legacy 11-feature ABI).
  //                   This is the layout `produceNew` and external code must
  //                   hand to `infer()`.
  //   • N_FEATURES  — network-facing input width AFTER `derive_features`
  //                   builds the augmented (sin/cos θXW, sin/cos 2θ_plane,
  //                   block_is_odd) columns. Must equal the `n_in` field
  //                   stored in the v2 weight-file header; mismatch is a
  //                   hard error.
  //
  // `N_FIXED` is the count of network-facing features that use the
  // InputLayer's fixed-scale path (everything except dQ/dx, which uses the
  // trainable log-compress scale).
  static constexpr int N_RAW        =  11;                    // caller ABI
  static constexpr int N_FEATURES   =  15;                    // network input
  static constexpr int N_CLASSES    =   4;                    // classes
  static constexpr int NEURON_N     = 128;                    // neuron_n
  static constexpr int NORM_GROUPS  =   8;                    // normgroups
  static constexpr int N_FIXED      =  14;                    // fixed input scales
  static constexpr int HEAD_HIDDEN  =  64;                    // v4 2-layer head hidden dim
  static constexpr float SQRT2PI    =   2.5066282746310002f;  // sqrt(2*pi)
  static constexpr float GN_EPS     = 1e-5f;                  // GroupNorm epsilon
  static constexpr float FLT_EPS    =   1.1754944e-38f;       // eps(Float32)
  static constexpr float DEG2RAD    = 0.017453292519943295f;  // π / 180

  // ICARUS wire-plane geometry (mirrors PhysicsHelpers.jl _Planeθ)
  static constexpr int   CHANNELS_PER_BLOCK = 13824;          // 4 blocks × 13824 = 55296
  static constexpr int   MAX_CHANNEL        = 55295;
  static constexpr float PLANE_THETA_1      = 0.0f;           // U plane
  static constexpr float PLANE_THETA_2      = -1.0471975511965976f;  // -π/3
  static constexpr float PLANE_THETA_3      =  1.0471975511965976f;  // +π/3

  // ICARUS Gain per plane
  // Seems not to be in the detector property service???
  static constexpr std::array<float, 3> GAIN = {0.016751, 0.012755, 0.012513}; // ADC per e

  // ═════════════════════════════════════════════════════════════════════════════
  // Weight file reader
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * TensorStore - loads all named tensors from a WAIRE1 binary weight file
   * and provides access by name.
   *
   * File format: see inference/README.md and Training.jl:export_weights_for_cpp.
   * All values are little-endian.
   *
   * v2 format adds an `n_in` UInt32 to the header between `version` and
   * `n_tensors`; this loader validates `n_in == N_FEATURES`.
   */
  class TensorStore {
    public:
      explicit TensorStore(const std::string& path) { load(path); }

      /** Return a const pointer to the raw Float32 data of tensor `name` */
      const float* data(const std::string& name) const
      {
        auto it = tensors_.find(name);
        if (it == tensors_.end())
          throw std::runtime_error("WaireMLod --- Tensor not found: " + name);
        return it->second.data();
      }

      /** Number of elements in tensor `name` */
      std::size_t size(const std::string& name) const
      {
        auto it = tensors_.find(name);
        if (it == tensors_.end())
          throw std::runtime_error("WaireMLod --- Tensor not found: " + name);
        return it->second.size();
      }

      /** WAIRE1 format version of the loaded file (2 = single-route, 3 = split-route) */
      uint32_t version() const { return version_; }

    private:
      std::map<std::string, std::vector<float>> tensors_;
      uint32_t version_ = 0;

      template<typename T>
        static T read_le(std::ifstream& f)
        {
          T val;
          f.read(reinterpret_cast<char*>(&val), sizeof(T));
          return val;
        }

      void load(const std::string& path)
      {
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("WaireMLod --- Cannot open weight file: " + path);

        // Magic
        char label[6];
        f.read(label, 6);
        if (std::strncmp(label, "WAIRE1", 6) != 0)
          throw std::runtime_error("WaireMLod --- Not a WAIRE1 file: " + path);

        uint32_t version = read_le<uint32_t>(f);
        if (version != 2 && version != 3 && version != 4 && version != 5 && version != 6)
        {
          std::string msg = "WaireMLod --- Unsupported WAIRE format version: " +
                            std::to_string(version) +
                            " (this reader expects version 2, 3, 4, 5, or 6). ";
          if (version == 1)
            msg += "Version-1 weight files were produced before the 15-feature "
                   "augmentation landed. Retrain with the current "
                   "src/Training.jl, or pin to a pre-augmentation build of "
                   "WireModInfer to consume legacy weights.";
          throw std::runtime_error(msg);
        }
        version_ = version;

        // v2/v3 headers both carry the network-facing input width n_in.
        // Validate against this build's compile-time N_FEATURES so a stale
        // binary never silently consumes a model with a different feature count.
        uint32_t n_in = read_le<uint32_t>(f);
        if (n_in != static_cast<uint32_t>(N_FEATURES))
          throw std::runtime_error(
            "WaireMLod --- Weight-file n_in (" + std::to_string(n_in) +
            ") does not match this binary's N_FEATURES (" +
            std::to_string(N_FEATURES) +
            "). Rebuild WireModInfer with matching N_FEATURES, or retrain "
            "with matching n_in in src/Training.jl.");

        uint32_t n_tensors = read_le<uint32_t>(f);

        for (uint32_t t = 0; t < n_tensors; ++t)
        {
          // Name
          uint32_t name_len = read_le<uint32_t>(f);
          std::string name(name_len, '\0');
          f.read(name.data(), name_len);

          // Shape
          uint8_t ndims = read_le<uint8_t>(f);
          std::size_t total = 1;
          for (uint8_t d = 0; d < ndims; ++d)
          {
            uint32_t dim = read_le<uint32_t>(f);
            total *= dim;
          }

          // Data
          std::vector<float> buf(total);
          f.read(reinterpret_cast<char*>(buf.data()), total * sizeof(float));
          tensors_[name] = std::move(buf);
        }

        if (!f) throw std::runtime_error("WaireMLod --- Unexpected end of weight file");
        std::cout << "Loaded " << tensors_.size() << " tensors from " << path << std::endl;
      }
  };

  // ═════════════════════════════════════════════════════════════════════════════
  // Activations
  // ═════════════════════════════════════════════════════════════════════════════
  inline float gelu(float x)
  {
    // tanh-approx GELU — MUST match the model's activation. Current Flux/NNlib
    // `gelu` is `gelu_tanh` (Dense layers show `typeof(gelu_tanh)`), so we use
    // the tanh form here, NOT the exact erf form (they differ ~0.1%, which
    // compounds through the network):
    //   x·0.5·(1 + tanh(√(2/π)·(x + 0.044715·x³)))
    const float l = 0.7978845608028654f;   // √(2/π)
    return 0.5f * x * (1.0f + std::tanh(l * (x + 0.044715f * x * x * x)));
  }

  inline Eigen::VectorXf gelu(const Eigen::VectorXf& v)
  {
    return v.unaryExpr([](float x){ return gelu(x); });
  }

  inline float softplus(float x)
  {
    // Numerically stable for large x
    return (x > 20.0f) ? x : std::log1p(std::exp(x));
  }

  inline Eigen::VectorXf softplus(const Eigen::VectorXf& v)
  {
    return v.unaryExpr([](float x){ return softplus(x); });
  }

  inline Eigen::VectorXf softmax(const Eigen::VectorXf& v)
  {
    Eigen::VectorXf shifted = v.array() - v.maxCoeff();
    Eigen::VectorXf ex      = shifted.array().exp();
    return ex / ex.sum();
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // GroupNorm (single sample, C channels, G groups)
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * Apply GroupNorm to a single vector of `channels` divided into `groups`
   * groups of size `channels / groups`.
   *
   * For each group g:
   *   mean_g = mean(x[g*gs ... (g+1)*gs - 1])
   *   var_g  = variance(x[g*gs ... (g+1)*gs - 1])
   *   x[c]   = (x[c] - mean_g) / sqrt(var_g + eps) for c in group g
   *   x[c]   = gamma[c] * x[c] + beta[c]
   *
   * Matches Flux.GroupNorm applied to a (channels, 1) tensor.
   *
   * gamma and beta are vectors of length `channels`.
   */
  inline Eigen::VectorXf group_norm(
    const Eigen::VectorXf& x,
    int channels,
    int groups,
    const float* gamma,
    const float* beta)
  {
    assert(channels % groups == 0);
    int group_size = channels / groups;
    Eigen::VectorXf out(channels);

    for (int g = 0; g < groups; ++g)
    {
      int start = g * group_size;
      // Mean
      float mean = 0.0f;
      for (int i = 0; i < group_size; ++i) mean += x[start + i];
      mean /= group_size;
      // Variance
      float var = 0.0f;
      for (int i = 0; i < group_size; ++i)
      {
        float d = x[start + i] - mean;
        var += d * d;
      }
      var /= group_size;
      float inv_std = 1.0f / std::sqrt(var + GN_EPS);
      // Normalize and affine
      for (int i = 0; i < group_size; ++i)
      {
        int c = start + i;
        out[c] = gamma[c] * (x[c] - mean) * inv_std + beta[c];
      }
    }

    return out;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // Dense
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * Dense layer forward pass: y = W * x + b, optionally followed by an activation.
   *
   * W is stored as a Julia (out, in) column-major matrix.
   * Eigen::Map<const MatrixXf>(ptr, out, in) reads it directly.
   */
  inline Eigen::VectorXf dense_linear(
    const Eigen::VectorXf& x,
    const float* W_ptr,
    const float* b_ptr,
    int out_dim,
    int in_dim)
  {
    Eigen::Map<const Eigen::MatrixXf> W(W_ptr, out_dim, in_dim);
    Eigen::Map<const Eigen::VectorXf> b(b_ptr, out_dim);
    return W * x + b;
  }

  inline Eigen::VectorXf dense_gelu(
    const Eigen::VectorXf& x,
    const float* W_ptr,
    const float* b_ptr,
    int out_dim,
    int in_dim)
  {
    return gelu(dense_linear(x, W_ptr, b_ptr, out_dim, in_dim));
  }

  inline Eigen::VectorXf dense_softplus(
    const Eigen::VectorXf& x,
    const float* W_ptr,
    const float* b_ptr,
    int out_dim,
    int in_dim)
  {
    return softplus(dense_linear(x, W_ptr, b_ptr, out_dim, in_dim));
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // InputLayer
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * InputLayer forward pass.
   *
   * Features 0..N_FIXED-1 are scaled by fixed scales.
   * Features N_FIXED..N_FEATURES-1 are compressed with log1p(|x * s|).
   *
   * In the production layout (N_FIXED=14, N_FEATURES=15) only the dQ/dx
   * feature at index 14 takes the trainable log-compress path.
   */
  inline Eigen::VectorXf input_layer(
    const std::array<float, N_FEATURES>& features,
    const float* scales_fixed,   // length N_FIXED
    const float* scales_train)   // length (N_FEATURES - N_FIXED)
  {
    Eigen::VectorXf out(N_FEATURES);
    for (int i = 0; i < N_FIXED; ++i)
      out[i] = features[i] * scales_fixed[i];
    for (int i = N_FIXED; i < N_FEATURES; ++i)
      out[i] = std::log1p(std::fabs(features[i] * scales_train[i - N_FIXED]));
    return out;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // Wire-plane geometry helper (mirrors PhysicsHelpers.jl _Planeθ)
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * plane_theta(channel) — wire-plane angle θ (radians) for the given global
   * wire channel. Each block of CHANNELS_PER_BLOCK = 13824 channels is split
   * into:
   *   [0,    2304)  → U plane (θ = 0)
   *   [2304, 8064)  → V plane (θ = ±π/3, sign alternates by block parity)
   *   [8064, 13824) → Y plane (θ = ∓π/3, opposite sign of V)
   *
   * NOTE: the .cc currently rotates dir_y/dir_z to plane-relative coordinates
   * via `wire_geom->WireAngleToVertical(...) - 0.5*M_PI`. That path predates
   * the production augmentation and may differ subtly from `plane_theta(ch)`
   * here. If validated discrepancies appear in test inferences, unify on
   * `plane_theta(channel)` for the .cc rotation as well. Out of scope for
   * this update.
   */
  inline float plane_theta(int channel)
  {
    if (channel < 0 || channel > MAX_CHANNEL)
      throw std::runtime_error(
        "WaireMLod --- plane_theta: channel " + std::to_string(channel) +
        " out of range [0, " + std::to_string(MAX_CHANNEL) + "]");

    int block            = channel / CHANNELS_PER_BLOCK;
    int channel_in_block = channel % CHANNELS_PER_BLOCK;
    bool block_is_even   = (block % 2 == 0);
    if (channel_in_block < 2304)
      return PLANE_THETA_1;
    if (channel_in_block < 8064)
      return block_is_even ? PLANE_THETA_2 : PLANE_THETA_3;
    return block_is_even ? PLANE_THETA_3 : PLANE_THETA_2;
  }

  /**
   * derive_features — build the 15-feature network input from the 11
   * caller-facing primitives. Embeds the augmentation that the Julia
   * training pipeline applies in src/DataPipeline.jl, so the external ABI
   * stays at the legacy 11 primitives.
   */
  inline std::array<float, N_FEATURES> derive_features(
    const std::array<float, N_RAW>& raw)
  {
    std::array<float, N_FEATURES> f{};

    // Pass-through primitives (positions, directions, plane-relative dirs,
    // channel as Float32). Indices 0..8 unchanged from caller layout.
    for (int i = 0; i < 9; ++i)
      f[i] = raw[i];

    // (sin, cos) θXW — caller-supplied θXW (degrees) directly. Matches
    // DataPipeline.jl which derives it once from (dir_x, pitch).
    float thxw_rad = raw[9] * DEG2RAD;
    f[9]  = std::sin(thxw_rad);
    f[10] = std::cos(thxw_rad);

    // (sin 2·θ_plane, cos 2·θ_plane) and block parity from channel.
    int channel_int = static_cast<int>(raw[8]);
    float theta_p   = plane_theta(channel_int);
    f[11] = std::sin(2.0f * theta_p);
    f[12] = std::cos(2.0f * theta_p);
    int block = channel_int / CHANNELS_PER_BLOCK;
    f[13] = (block % 2 != 0) ? 1.0f : 0.0f;

    // dQ/dx — last feature, takes the trainable log-compress path inside
    // InputLayer.
    f[14] = raw[10];

    return f;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // Slicing operations (mirror SlicingOps.jl, 0-indexed)
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * _astream_select: x[0 : N_CLASSES + NEURON_N - 1]
   * Drops the ŵ-stream rows. Input (260,) → output (132,).
   */
  inline Eigen::VectorXf astream_select(const Eigen::VectorXf& x)
  { return x.head(N_CLASSES + NEURON_N); }

  /**
   * _wstream_select: vcat(x[0:N_CLASSES-1], x[N_CLASSES+NEURON_N : end])
   * Drops the Â-stream rows. Input (260,) → output (132,).
   */
  inline Eigen::VectorXf wstream_select(const Eigen::VectorXf& x)
  {
    int total_out = N_CLASSES + NEURON_N;
    Eigen::VectorXf out(total_out);
    out.head(N_CLASSES) = x.head(N_CLASSES);
    out.tail(NEURON_N)  = x.tail(NEURON_N);
    return out;
  }

  /**
   * _prepend_classes: vcat(x[0:N_CLASSES-1], cwl_out)
   * Used as the SkipConnection combiner after the first CWL.
   * x:       (132,) classifier output
   * cwl_out: (256,) CWL output
   * result:  (260,)
   */
  inline Eigen::VectorXf prepend_classes(const Eigen::VectorXf& cwl_out, const Eigen::VectorXf& x)
  {
    int out_size = N_CLASSES + cwl_out.size();
    Eigen::VectorXf out(out_size);
    out.head(N_CLASSES) = x.head(N_CLASSES);
    out.tail(cwl_out.size()) = cwl_out;
    return out;
  }

  /**
   * _residual_combine: vcat(x[0:N_CLASSES-1], x[N_CLASSES:end] + residual)
   * SkipConnection combiner in the two parallel residual blocks.
   * residual: (256,)  parallel block output
   * x:        (260,)  working tensor
   * result:   (260,)
   */
  inline Eigen::VectorXf residual_combine(const Eigen::VectorXf& residual, const Eigen::VectorXf& x)
  {
    Eigen::VectorXf out = x;
    out.tail(out.size() - N_CLASSES) += residual;
    return out;
  }

  /**
   * _output_combine: vcat(sqrt2pi * w_hat * A_hat, w_hat)
   * Converts the (Â, ŵ) pair to the (Î, ŵ) output pair.
   * A_hat, w_hat are scalars (size-1 vectors). Result (2,): [Î, ŵ].
   */
  inline Eigen::VectorXf output_combine(const Eigen::VectorXf& A_hat, const Eigen::VectorXf& w_hat)
  {
    Eigen::VectorXf out(2);
    out[0] = SQRT2PI * w_hat[0] * A_hat[0];
    out[1] = w_hat[0];
    return out;
  }

  /**
   * _output_combine_log: vcat(sqrt2pi * ŵ * Â, ŵ) where
   *   Â = exp(clamp(logA, -50, 30)),  ŵ = exp(clamp(logW, -50, 30)).
   * v3 log-output variant: the heads emit log-space scalars (identity
   * activation), so the exp + clamp lives here. Matches the Julia
   * `_output_combine_log` (clamp bounds and √(2π) constant inlined there).
   * logA, logW are scalars (size-1 vectors). Result (2,): [Î, ŵ].
   */
  inline Eigen::VectorXf output_combine_log(const Eigen::VectorXf& logA, const Eigen::VectorXf& logW)
  {
    float A_hat = std::exp(std::max(-50.0f, std::min(30.0f, logA[0])));
    float w_hat = std::exp(std::max(-50.0f, std::min(30.0f, logW[0])));
    Eigen::VectorXf out(2);
    out[0] = SQRT2PI * w_hat * A_hat;
    out[1] = w_hat;
    return out;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // ClassWeightedLayer (single-sample forward pass, flat weight layout)
  // ═════════════════════════════════════════════════════════════════════════════
  inline Eigen::VectorXf cwl_forward(
    const Eigen::VectorXf& x,
    int n_classes,
    int dense_in,
    int dense_out,
    int norm_groups,
    const float* weights,
    const float* biases,
    const float* norm_gamma,
    const float* norm_beta,
    const float* expert_scale,
    const float* expert_bias)
  {
    int flat_out = dense_out * n_classes;

    // Split input
    Eigen::Map<const Eigen::VectorXf> cw(x.data(), n_classes);
    Eigen::Map<const Eigen::VectorXf> neurons(x.data() + n_classes, dense_in);

    // Single GEMM: weights * neurons + biases
    Eigen::Map<const Eigen::MatrixXf> W(weights, flat_out, dense_in);
    Eigen::Map<const Eigen::VectorXf> b(biases, flat_out);
    Eigen::VectorXf act_flat = gelu(W * neurons + b);

    // GroupNorm + per-expert affine + weighted mixture
    Eigen::VectorXf mixed = Eigen::VectorXf::Zero(dense_out);

    for (int k = 0; k < n_classes; ++k)
    {
      int offset = k * dense_out;

      Eigen::VectorXf act_k = act_flat.segment(offset, dense_out);

      Eigen::VectorXf normed_k = group_norm(act_k, dense_out, norm_groups, norm_gamma, norm_beta);

      Eigen::Map<const Eigen::VectorXf> es_k(expert_scale + offset, dense_out);
      Eigen::Map<const Eigen::VectorXf> eb_k(expert_bias  + offset, dense_out);
      normed_k = normed_k.cwiseProduct(es_k) + eb_k;

      mixed += normed_k * cw[k];
    }

    // L2-normalise by class weight magnitude
    float cw_norm_sq = cw.squaredNorm();
    float scale = 1.0f / std::sqrt(cw_norm_sq + FLT_EPS);
    return scale * mixed;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // WaireMLod  —  full forward pass
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * WaireMLod
   *
   * Loads all weights from a v2, v3, or v4 WAIRE1 binary and infers data-like
   * hit integrals and widths. `infer` dispatches on the loaded format version:
   * v2 runs the single-routing classifier + one MoE pass; v3/v4 run the shared
   * backbone, two routers, and one MoE pass per router with the log-output
   * combine (v4 uses 2-layer output heads; see the file header for the version
   * semantics). External callers pass the legacy 11-feature primitive layout
   * (see file header); the augmentation to the 15-feature network input happens
   * inside `infer`.
   */
  class WaireMLod
  {
    public:
      explicit WaireMLod(const std::string& weights_path)
        : store_(weights_path) { validate_shapes(); cache_pointers(); }

      /**
       * Run inference on a single hit.
       *
       * @param raw  Array of N_RAW = 11 caller-facing primitive features in
       *             the order documented at the top of this header. Internally
       *             augmented to the 15-feature network input via
       *             `derive_features`.
       * @return     {predicted_integral, predicted_width}
       */
      std::array<float, 2> infer(const std::array<float, N_RAW>& raw) const
      {
        std::array<float, N_FEATURES> features = derive_features(raw);

        if (store_.version() == 6)
        {
          // ── v6 decoupled: shared backbone → TWO independent expert banks
          //    (routing_a→bank_A→logÂ, routing_w→bank_W→logŵ); single-pass
          //    combine I=√(2π)·Â·ŵ, w=ŵ.
          return infer_v6(features);
        }

        if (store_.version() == 5)
        {
          // ── v5 asymmetric: one router → A-expert bank → logÂ; dense w-head
          //    → logŵ; single-pass combine I=√(2π)·Â·ŵ, w=ŵ.
          return infer_v5(features);
        }

        if (store_.version() >= 3)
        {
          // ── v3/v4 split-routing ───────────────────────────────────────────
          // Shared backbone emits 128-d features; two routers each prepend
          // their softmax to form a 132-d MoE input. The MoE runs once per
          // router (log-output combine). Take Î from the A-routed run and ŵ
          // from the W-routed run. v4 differs only in the output-head depth,
          // handled inside run_mixture_of_experts.
          Eigen::VectorXf feats = run_backbone(features);

          Eigen::VectorXf rA = softmax(dense_linear(
            feats, wp_.routing_a_w, wp_.routing_a_b, N_CLASSES, NEURON_N));
          Eigen::VectorXf rW = softmax(dense_linear(
            feats, wp_.routing_w_w, wp_.routing_w_b, N_CLASSES, NEURON_N));

          Eigen::VectorXf inA(N_CLASSES + NEURON_N);
          inA.head(N_CLASSES) = rA;
          inA.tail(NEURON_N)  = feats;
          Eigen::VectorXf inW(N_CLASSES + NEURON_N);
          inW.head(N_CLASSES) = rW;
          inW.tail(NEURON_N)  = feats;

          Eigen::VectorXf outA = run_mixture_of_experts(inA, /*log_output=*/true);
          Eigen::VectorXf outW = run_mixture_of_experts(inW, /*log_output=*/true);
          return {outA[0], outW[1]};
        }

        // ── v2 single-routing ───────────────────────────────────────────────
        Eigen::VectorXf x = run_classifier(features);
        Eigen::VectorXf y = run_mixture_of_experts(x);
        return {y[0], y[1]};
      }

      /**
       * Run inference on a batch of hits.
       *
       * @param raw_vec  Vector of N_RAW = 11 primitive-feature arrays.
       * @return         vector<{predicted_integral, predicted_width}>
       */
      std::vector<std::array<float, 2>>
        infer_batch(const std::vector<std::array<float, N_RAW>>& raw_vec) const
        {
          std::vector<std::array<float, 2>> out;
          out.reserve(raw_vec.size());
          for (auto const& raw : raw_vec) out.push_back(infer(raw));
          return out;
        }

      /**
       * Construct a vector of data-like hits from simulated hits.
       */
      std::vector<recob::Hit> produceNew(
        const std::vector<art::Ptr<recob::Hit>> old_hits,
        const cheat::BackTrackerService* back_tracker,
        const cheat::ParticleInventoryService* particles,
        const detinfo::DetectorClocksData* det_clock,
        const detinfo::DetectorPropertiesData* det_prop,
        const geo::WireReadoutGeom* wire_geom) const;

    private:
      TensorStore store_;

      /**
       * Cache pointers to the weight tensors in the store.
       */
      struct WeightPointers
      {
        const float *il_scales_fixed, *il_scales_train;
        const float *d1_w, *d1_b;
        const float *rb1_gamma, *rb1_beta, *rb1_w, *rb1_b;
        const float *rb2_gamma, *rb2_beta, *rb2_w, *rb2_b;
        const float *cls_w = nullptr, *cls_b = nullptr;                 // v2 only
        const float *routing_a_w = nullptr, *routing_a_b = nullptr;     // v3/v4 only
        const float *routing_w_w = nullptr, *routing_w_b = nullptr;     // v3/v4 only
        struct CWLPtrs
        {
          const float *weights, *biases;
          const float *gamma, *beta;
          const float *escale, *ebias;
        };
        CWLPtrs cwl1, cwl2a, cwl2w, cwl3a, cwl3w;
        // Output heads. v2/v3: single Dense, use *_w / *_b only.
        // v4: 2-layer MLP, use *_l1_* and *_l2_*.
        const float *ahead_w = nullptr, *ahead_b = nullptr;             // v2/v3
        const float *whead_w = nullptr, *whead_b = nullptr;             // v2/v3
        const float *ahead_l1_w = nullptr, *ahead_l1_b = nullptr;       // v4
        const float *ahead_l2_w = nullptr, *ahead_l2_b = nullptr;       // v4
        const float *whead_l1_w = nullptr, *whead_l1_b = nullptr;       // v4/v5
        const float *whead_l2_w = nullptr, *whead_l2_b = nullptr;       // v4/v5
        const float *whead_l3_w = nullptr, *whead_l3_b = nullptr;       // v5 (3-layer dense w-head)
        CWLPtrs cwl_a1, cwl_a2, cwl_a3;                                 // v5/v6 A-expert bank
        CWLPtrs cwl_w1, cwl_w2, cwl_w3;                                 // v6 W-expert bank
        int n_exp = N_CLASSES;                                          // v5/v6: runtime expert count
      };
      WeightPointers wp_;

      void cache_pointers()
      {
        wp_.il_scales_fixed = store_.data("il.scales_fixed");
        wp_.il_scales_train = store_.data("il.scales_train");
        wp_.d1_w  = store_.data("d1.weight");
        wp_.d1_b  = store_.data("d1.bias");
        wp_.rb1_gamma = store_.data("rb1.gn.gamma");
        wp_.rb1_beta  = store_.data("rb1.gn.beta");
        wp_.rb1_w     = store_.data("rb1.d.weight");
        wp_.rb1_b     = store_.data("rb1.d.bias");
        wp_.rb2_gamma = store_.data("rb2.gn.gamma");
        wp_.rb2_beta  = store_.data("rb2.gn.beta");
        wp_.rb2_w     = store_.data("rb2.d.weight");
        wp_.rb2_b     = store_.data("rb2.d.bias");

        // ── v6 (decoupled): shared backbone + TWO independent expert banks. Same
        //    A-bank as v5, plus a W-bank (routing_w + cwl_w1/2/3 + 2-layer log-ŵ
        //    head). K read from routing_a (same K for both banks). Early-return —
        //    the v2/3/4/5 tensors below do NOT all exist in a v6 file. ─────────
        if (store_.version() == 6)
        {
          auto load_cwl6 = [&](const std::string& p) {
            return WeightPointers::CWLPtrs{
              store_.data(p + ".weights"),      store_.data(p + ".biases"),
              store_.data(p + ".norm.gamma"),   store_.data(p + ".norm.beta"),
              store_.data(p + ".expert_scale"), store_.data(p + ".expert_bias") };
          };
          wp_.n_exp       = static_cast<int>(store_.size("routing_a.bias"));
          wp_.routing_a_w = store_.data("routing_a.weight");
          wp_.routing_a_b = store_.data("routing_a.bias");
          wp_.cwl_a1 = load_cwl6("cwl_a1"); wp_.cwl_a2 = load_cwl6("cwl_a2"); wp_.cwl_a3 = load_cwl6("cwl_a3");
          wp_.ahead_l1_w = store_.data("ahead.l1.weight"); wp_.ahead_l1_b = store_.data("ahead.l1.bias");
          wp_.ahead_l2_w = store_.data("ahead.l2.weight"); wp_.ahead_l2_b = store_.data("ahead.l2.bias");
          wp_.routing_w_w = store_.data("routing_w.weight");
          wp_.routing_w_b = store_.data("routing_w.bias");
          wp_.cwl_w1 = load_cwl6("cwl_w1"); wp_.cwl_w2 = load_cwl6("cwl_w2"); wp_.cwl_w3 = load_cwl6("cwl_w3");
          wp_.whead_l1_w = store_.data("whead.l1.weight"); wp_.whead_l1_b = store_.data("whead.l1.bias");
          wp_.whead_l2_w = store_.data("whead.l2.weight"); wp_.whead_l2_b = store_.data("whead.l2.bias");
          return;
        }

        // ── v5 (asymmetric): one router, one A-expert bank (3 CWLs), a 2-layer
        //    log-Â head, and a 3-layer dense log-ŵ head. Expert count is read
        //    from routing_a (no compile-time class assumption). Early-return —
        //    the v2/3/4 tensors below (cls / routing_w / cwl1.. / single heads)
        //    do NOT exist in a v5 file. ──────────────────────────────────────
        if (store_.version() == 5)
        {
          auto load_cwl5 = [&](const std::string& p) {
            return WeightPointers::CWLPtrs{
              store_.data(p + ".weights"),      store_.data(p + ".biases"),
              store_.data(p + ".norm.gamma"),   store_.data(p + ".norm.beta"),
              store_.data(p + ".expert_scale"), store_.data(p + ".expert_bias") };
          };
          wp_.n_exp       = static_cast<int>(store_.size("routing_a.bias"));
          wp_.routing_a_w = store_.data("routing_a.weight");
          wp_.routing_a_b = store_.data("routing_a.bias");
          wp_.cwl_a1 = load_cwl5("cwl_a1");
          wp_.cwl_a2 = load_cwl5("cwl_a2");
          wp_.cwl_a3 = load_cwl5("cwl_a3");
          wp_.ahead_l1_w = store_.data("ahead.l1.weight"); wp_.ahead_l1_b = store_.data("ahead.l1.bias");
          wp_.ahead_l2_w = store_.data("ahead.l2.weight"); wp_.ahead_l2_b = store_.data("ahead.l2.bias");
          wp_.whead_l1_w = store_.data("whead.l1.weight"); wp_.whead_l1_b = store_.data("whead.l1.bias");
          wp_.whead_l2_w = store_.data("whead.l2.weight"); wp_.whead_l2_b = store_.data("whead.l2.bias");
          wp_.whead_l3_w = store_.data("whead.l3.weight"); wp_.whead_l3_b = store_.data("whead.l3.bias");
          return;
        }

        // Routing head(s): v2 has a single `cls` head; v3/v4 split it into two
        // independent routers (`routing_a` for the integral, `routing_w` for
        // the width).
        if (store_.version() >= 3)
        {
          wp_.routing_a_w = store_.data("routing_a.weight");
          wp_.routing_a_b = store_.data("routing_a.bias");
          wp_.routing_w_w = store_.data("routing_w.weight");
          wp_.routing_w_b = store_.data("routing_w.bias");
        }
        else
        {
          wp_.cls_w = store_.data("cls.weight");
          wp_.cls_b = store_.data("cls.bias");
        }

        auto load_cwl = [&](const std::string& p)
        {
          return WeightPointers::CWLPtrs
          {
            store_.data(p + ".weights"),
            store_.data(p + ".biases"),
            store_.data(p + ".norm.gamma"),
            store_.data(p + ".norm.beta"),
            store_.data(p + ".expert_scale"),
            store_.data(p + ".expert_bias"),
          };
        };
        wp_.cwl1  = load_cwl("cwl1");
        wp_.cwl2a = load_cwl("cwl2a");
        wp_.cwl2w = load_cwl("cwl2w");
        wp_.cwl3a = load_cwl("cwl3a");
        wp_.cwl3w = load_cwl("cwl3w");

        // Output heads. v2/v3: single Dense per head. v4: 2-layer MLP per head.
        if (store_.version() == 4)
        {
          wp_.ahead_l1_w = store_.data("ahead.l1.weight");
          wp_.ahead_l1_b = store_.data("ahead.l1.bias");
          wp_.ahead_l2_w = store_.data("ahead.l2.weight");
          wp_.ahead_l2_b = store_.data("ahead.l2.bias");
          wp_.whead_l1_w = store_.data("whead.l1.weight");
          wp_.whead_l1_b = store_.data("whead.l1.bias");
          wp_.whead_l2_w = store_.data("whead.l2.weight");
          wp_.whead_l2_b = store_.data("whead.l2.bias");
        }
        else
        {
          wp_.ahead_w = store_.data("ahead.weight");
          wp_.ahead_b = store_.data("ahead.bias");
          wp_.whead_w = store_.data("whead.weight");
          wp_.whead_b = store_.data("whead.bias");
        }
      }

      /**
       * Validate that all loaded tensors have the expected number of elements.
       * Uses N_FEATURES (15) for the InputLayer and first Dense-layer shapes,
       * matching the v2 weight format.
       */
      void validate_shapes() const
      {
        auto check = [&](const std::string& name, std::size_t expected)
        {
          std::size_t actual = store_.size(name);
          if (actual != expected)
            throw std::runtime_error(
              "WaireMLod --- Shape mismatch for tensor '" + name + "' expected " +
              std::to_string(expected) + " elements, got " + std::to_string(actual));
        };

        // Classifier
        check("il.scales_fixed",  N_FIXED);
        check("il.scales_train",  N_FEATURES - N_FIXED);
        check("d1.weight",        NEURON_N * N_FEATURES);
        check("d1.bias",          NEURON_N);
        check("rb1.gn.gamma",     NEURON_N);
        check("rb1.gn.beta",      NEURON_N);
        check("rb1.d.weight",     NEURON_N * NEURON_N);
        check("rb1.d.bias",       NEURON_N);
        check("rb2.gn.gamma",     NEURON_N);
        check("rb2.gn.beta",      NEURON_N);
        check("rb2.d.weight",     NEURON_N * NEURON_N);
        check("rb2.d.bias",       NEURON_N);

        // ── v6 (decoupled): backbone + A-bank + W-bank. Each bank: a router, 3
        //    CWLs (dense_out=NEURON_N), and a 2-layer log head. K = expert count
        //    (same for both banks, read from routing_a). ──────────────────────
        if (store_.version() == 6)
        {
          int K = static_cast<int>(store_.size("routing_a.bias"));
          for (const auto& rt : {"routing_a", "routing_w"})
          {
            std::string p(rt);
            check(p + ".weight", K * NEURON_N);
            check(p + ".bias",   K);
          }
          for (const auto& prefix : {"cwl_a1", "cwl_a2", "cwl_a3", "cwl_w1", "cwl_w2", "cwl_w3"})
          {
            std::string p(prefix);
            check(p + ".weights",      NEURON_N * K * NEURON_N);
            check(p + ".biases",       NEURON_N * K);
            check(p + ".norm.gamma",   NEURON_N);
            check(p + ".norm.beta",    NEURON_N);
            check(p + ".expert_scale", NEURON_N * K);
            check(p + ".expert_bias",  NEURON_N * K);
          }
          int head_in = K + NEURON_N;
          for (const auto& hd : {"ahead", "whead"})
          {
            std::string p(hd);
            check(p + ".l1.weight", HEAD_HIDDEN * head_in); check(p + ".l1.bias", HEAD_HIDDEN);
            check(p + ".l2.weight", HEAD_HIDDEN);           check(p + ".l2.bias", 1);
          }
          return;
        }

        // ── v5 (asymmetric): one router + A-expert bank (3 CWLs, dense_out=NEURON_N)
        //    + 2-layer log-Â head + 3-layer dense log-ŵ head. K = expert count. ──
        if (store_.version() == 5)
        {
          int K = static_cast<int>(store_.size("routing_a.bias"));
          check("routing_a.weight", K * NEURON_N);
          check("routing_a.bias",   K);
          for (const auto& prefix : {"cwl_a1", "cwl_a2", "cwl_a3"})
          {
            std::string p(prefix);
            check(p + ".weights",      NEURON_N * K * NEURON_N);
            check(p + ".biases",       NEURON_N * K);
            check(p + ".norm.gamma",   NEURON_N);
            check(p + ".norm.beta",    NEURON_N);
            check(p + ".expert_scale", NEURON_N * K);
            check(p + ".expert_bias",  NEURON_N * K);
          }
          int head_in = K + NEURON_N;
          check("ahead.l1.weight", HEAD_HIDDEN * head_in);  check("ahead.l1.bias", HEAD_HIDDEN);
          check("ahead.l2.weight", HEAD_HIDDEN);            check("ahead.l2.bias", 1);
          check("whead.l1.weight", HEAD_HIDDEN * NEURON_N); check("whead.l1.bias", HEAD_HIDDEN);
          check("whead.l2.weight", HEAD_HIDDEN * HEAD_HIDDEN); check("whead.l2.bias", HEAD_HIDDEN);
          check("whead.l3.weight", HEAD_HIDDEN);            check("whead.l3.bias", 1);
          return;
        }

        // Routing head(s): single `cls` for v2, dual `routing_a`/`routing_w`
        // for v3/v4. Each head is a Dense(NEURON_N → N_CLASSES).
        if (store_.version() >= 3)
        {
          check("routing_a.weight", N_CLASSES * NEURON_N);
          check("routing_a.bias",   N_CLASSES);
          check("routing_w.weight", N_CLASSES * NEURON_N);
          check("routing_w.bias",   N_CLASSES);
        }
        else
        {
          check("cls.weight",       N_CLASSES * NEURON_N);
          check("cls.bias",         N_CLASSES);
        }

        // MoE: CWL1 (dense_out = 2*NEURON_N)
        int cwl1_out = 2 * NEURON_N;
        check("cwl1.weights",      cwl1_out * N_CLASSES * NEURON_N);
        check("cwl1.biases",       cwl1_out * N_CLASSES);
        check("cwl1.norm.gamma",   cwl1_out);
        check("cwl1.norm.beta",    cwl1_out);
        check("cwl1.expert_scale", cwl1_out * N_CLASSES);
        check("cwl1.expert_bias",  cwl1_out * N_CLASSES);

        // MoE: CWL2a / CWL2w / CWL3a / CWL3w (dense_out = NEURON_N)
        for (const auto& prefix : {"cwl2a", "cwl2w", "cwl3a", "cwl3w"})
        {
          std::string p(prefix);
          check(p + ".weights",      NEURON_N * N_CLASSES * NEURON_N);
          check(p + ".biases",       NEURON_N * N_CLASSES);
          check(p + ".norm.gamma",   NEURON_N);
          check(p + ".norm.beta",    NEURON_N);
          check(p + ".expert_scale", NEURON_N * N_CLASSES);
          check(p + ".expert_bias",  NEURON_N * N_CLASSES);
        }

        // Output heads. v2/v3: single Dense(head_in → 1). v4: 2-layer MLP
        // Dense(head_in → HEAD_HIDDEN, gelu) → Dense(HEAD_HIDDEN → 1, id).
        int head_in = N_CLASSES + NEURON_N;
        if (store_.version() == 4)
        {
          check("ahead.l1.weight", HEAD_HIDDEN * head_in);
          check("ahead.l1.bias",   HEAD_HIDDEN);
          check("ahead.l2.weight", HEAD_HIDDEN);
          check("ahead.l2.bias",   1);
          check("whead.l1.weight", HEAD_HIDDEN * head_in);
          check("whead.l1.bias",   HEAD_HIDDEN);
          check("whead.l2.weight", HEAD_HIDDEN);
          check("whead.l2.bias",   1);
        }
        else
        {
          check("ahead.weight", head_in);
          check("ahead.bias",   1);
          check("whead.weight", head_in);
          check("whead.bias",   1);
        }
      }

      /**
       * Shared classifier backbone (routing-head agnostic).
       *   InputLayer → Dense(N_FEATURES → 128, gelu) → ResBlock1 → ResBlock2
       *   Output: (128,) activations.
       *
       * v2 wraps this with the single `cls` softmax head (see run_classifier);
       * v3 feeds it to the two independent routers in `infer`.
       */
      Eigen::VectorXf run_backbone(const std::array<float, N_FEATURES>& features) const
      {
        // InputLayer
        Eigen::VectorXf x = input_layer(features,
          wp_.il_scales_fixed,
          wp_.il_scales_train);

        // Dense N_FEATURES → NEURON_N, gelu
        x = dense_gelu(x, wp_.d1_w, wp_.d1_b, NEURON_N, N_FEATURES);

        // ResBlock 1: GroupNorm → Dense(128→128, gelu) → residual
        {
          Eigen::VectorXf normed = group_norm(x, NEURON_N, NORM_GROUPS,
            wp_.rb1_gamma, wp_.rb1_beta);
          Eigen::VectorXf inner  = dense_gelu(normed,
            wp_.rb1_w, wp_.rb1_b,
            NEURON_N, NEURON_N);
          x += inner;
        }

        // ResBlock 2: same structure
        {
          Eigen::VectorXf normed = group_norm(x, NEURON_N, NORM_GROUPS,
            wp_.rb2_gamma, wp_.rb2_beta);
          Eigen::VectorXf inner  = dense_gelu(normed,
            wp_.rb2_w, wp_.rb2_b,
            NEURON_N, NEURON_N);
          x += inner;
        }

        return x;
      }

      /**
       * v5 asymmetric forward. Backbone (128-d) → A-router → A-expert bank
       * (first CWL via _prepend_classes, then two residual CWLs) → 2-layer log-Â
       * head; a 3-layer dense log-ŵ head reads the backbone features directly.
       * Combine: Â=exp(clamp(logÂ,-50,30)), ŵ=exp(clamp(logŵ,-50,30));
       * I=√(2π)·Â·ŵ, w=ŵ. Expert count K = wp_.n_exp (read from routing_a).
       */
      std::array<float, 2> infer_v5(const std::array<float, N_FEATURES>& features) const
      {
        const int K = wp_.n_exp;
        Eigen::VectorXf feats = run_backbone(features);                 // (128,)

        Eigen::VectorXf rA = softmax(dense_linear(
          feats, wp_.routing_a_w, wp_.routing_a_b, K, NEURON_N));        // (K,)

        // A-expert bank. Working tensor h = [rA ; neurons]; neurons start = feats.
        Eigen::VectorXf h(K + NEURON_N);
        h.head(K) = rA;
        h.tail(NEURON_N) = feats;
        auto run_cwl = [&](const WeightPointers::CWLPtrs& c) {
          return cwl_forward(h, K, NEURON_N, NEURON_N, NORM_GROUPS,
            c.weights, c.biases, c.gamma, c.beta, c.escale, c.ebias);
        };
        Eigen::VectorXf c1 = run_cwl(wp_.cwl_a1); h.tail(NEURON_N)  = c1;  // _prepend_classes
        Eigen::VectorXf c2 = run_cwl(wp_.cwl_a2); h.tail(NEURON_N) += c2;  // _residual_combine
        Eigen::VectorXf c3 = run_cwl(wp_.cwl_a3); h.tail(NEURON_N) += c3;  // _residual_combine

        // log-Â head: Dense(K+128 → 64, gelu) → Dense(64 → 1, identity).
        Eigen::VectorXf ah   = dense_gelu(h, wp_.ahead_l1_w, wp_.ahead_l1_b, HEAD_HIDDEN, K + NEURON_N);
        Eigen::VectorXf logA = dense_linear(ah, wp_.ahead_l2_w, wp_.ahead_l2_b, 1, HEAD_HIDDEN);

        // log-ŵ head (dense, off the backbone): Dense(128→64,gelu)→Dense(64→64,gelu)→Dense(64→1,id).
        Eigen::VectorXf w1   = dense_gelu(feats, wp_.whead_l1_w, wp_.whead_l1_b, HEAD_HIDDEN, NEURON_N);
        Eigen::VectorXf w2   = dense_gelu(w1,    wp_.whead_l2_w, wp_.whead_l2_b, HEAD_HIDDEN, HEAD_HIDDEN);
        Eigen::VectorXf logW = dense_linear(w2,  wp_.whead_l3_w, wp_.whead_l3_b, 1, HEAD_HIDDEN);

        Eigen::VectorXf out = output_combine_log(logA, logW);
        return {out[0], out[1]};
      }

      // Routed expert bank: routing softmax → 3 residual CWLs (prepend, then two
      // residual-combines) → 2-layer log head. Identical structure to v5's A-bank;
      // v6 runs it twice (A→logÂ, W→logŵ). Returns the (1,) log output.
      Eigen::VectorXf run_expert_bank(
        const Eigen::VectorXf& feats, int K,
        const float* route_w, const float* route_b,
        const WeightPointers::CWLPtrs& c1,
        const WeightPointers::CWLPtrs& c2,
        const WeightPointers::CWLPtrs& c3,
        const float* h1_w, const float* h1_b,
        const float* h2_w, const float* h2_b) const
      {
        Eigen::VectorXf r = softmax(dense_linear(feats, route_w, route_b, K, NEURON_N));  // (K,)
        Eigen::VectorXf h(K + NEURON_N);
        h.head(K) = r;
        h.tail(NEURON_N) = feats;
        auto run_cwl = [&](const WeightPointers::CWLPtrs& c) {
          return cwl_forward(h, K, NEURON_N, NEURON_N, NORM_GROUPS,
            c.weights, c.biases, c.gamma, c.beta, c.escale, c.ebias);
        };
        Eigen::VectorXf x1 = run_cwl(c1); h.tail(NEURON_N)  = x1;   // _prepend_classes
        Eigen::VectorXf x2 = run_cwl(c2); h.tail(NEURON_N) += x2;   // _residual_combine
        Eigen::VectorXf x3 = run_cwl(c3); h.tail(NEURON_N) += x3;   // _residual_combine
        Eigen::VectorXf hh = dense_gelu(h, h1_w, h1_b, HEAD_HIDDEN, K + NEURON_N);
        return dense_linear(hh, h2_w, h2_b, 1, HEAD_HIDDEN);        // (1,) log output
      }

      std::array<float, 2> infer_v6(const std::array<float, N_FEATURES>& features) const
      {
        const int K = wp_.n_exp;
        Eigen::VectorXf feats = run_backbone(features);             // (128,)
        Eigen::VectorXf logA = run_expert_bank(feats, K,
          wp_.routing_a_w, wp_.routing_a_b, wp_.cwl_a1, wp_.cwl_a2, wp_.cwl_a3,
          wp_.ahead_l1_w, wp_.ahead_l1_b, wp_.ahead_l2_w, wp_.ahead_l2_b);  // logÂ
        Eigen::VectorXf logW = run_expert_bank(feats, K,
          wp_.routing_w_w, wp_.routing_w_b, wp_.cwl_w1, wp_.cwl_w2, wp_.cwl_w3,
          wp_.whead_l1_w, wp_.whead_l1_b, wp_.whead_l2_w, wp_.whead_l2_b);  // logŵ
        Eigen::VectorXf out = output_combine_log(logA, logW);
        return {out[0], out[1]};
      }

      /**
       * Classifier subnetwork (v2 single-routing).
       *   run_backbone → SkipConnection(Dense(128→4)→softmax, vcat)
       *   Output: (132,)  rows 0..3 = softmax, rows 4..131 = activations.
       */
      Eigen::VectorXf run_classifier(const std::array<float, N_FEATURES>& features) const
      {
        Eigen::VectorXf x = run_backbone(features);

        // Classifier head
        Eigen::VectorXf cls_scores = dense_linear(x,
          wp_.cls_w, wp_.cls_b,
          N_CLASSES, NEURON_N);
        cls_scores = softmax(cls_scores);
        Eigen::VectorXf out(N_CLASSES + NEURON_N);
        out.head(N_CLASSES) = cls_scores;
        out.tail(NEURON_N)  = x;

        return out;
      }

      /**
       * MixtureOfExperts subnetwork.
       *   Input (132,) → CWL1 → (260,)
       *   → ResBlock3 (Parallel Â/ŵ chains) → (260,)
       *   → ResBlock4 (Parallel Â/ŵ chains) → (260,)
       *   → Parallel(combine, Â-head, ŵ-head) → (2,)
       *
       * @param log_output  false (v2): heads use softplus, combine is
       *                     output_combine. true (v3/v4): heads use identity
       *                     activation and combine is output_combine_log
       *                     (exp + clamp lives in the combine). v4 heads are
       *                     2-layer MLPs; v3 heads are a single Dense.
       */
      Eigen::VectorXf run_mixture_of_experts(const Eigen::VectorXf& clf_out,
                                             bool log_output = false) const
      {
        // CWL1: (132,) → (256,), then prepend the 4 class scores → (260,)
        Eigen::VectorXf cwl1_out = cwl_forward(
          clf_out, N_CLASSES, NEURON_N, 2*NEURON_N, NORM_GROUPS,
          wp_.cwl1.weights, wp_.cwl1.biases,
          wp_.cwl1.gamma,   wp_.cwl1.beta,
          wp_.cwl1.escale,  wp_.cwl1.ebias);
        Eigen::VectorXf x = prepend_classes(cwl1_out, clf_out);

        // ResBlock 3: Parallel Â and ŵ chains
        {
          Eigen::VectorXf a_out = cwl_forward(
            astream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            wp_.cwl2a.weights, wp_.cwl2a.biases,
            wp_.cwl2a.gamma,   wp_.cwl2a.beta,
            wp_.cwl2a.escale,  wp_.cwl2a.ebias);

          Eigen::VectorXf w_out = cwl_forward(
            wstream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            wp_.cwl2w.weights, wp_.cwl2w.biases,
            wp_.cwl2w.gamma,   wp_.cwl2w.beta,
            wp_.cwl2w.escale,  wp_.cwl2w.ebias);

          Eigen::VectorXf parallel_out(2 * NEURON_N);
          parallel_out.head(NEURON_N) = a_out;
          parallel_out.tail(NEURON_N) = w_out;

          x = residual_combine(parallel_out, x);
        }

        // ResBlock 4: identical structure
        {
          Eigen::VectorXf a_out = cwl_forward(
            astream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            wp_.cwl3a.weights, wp_.cwl3a.biases,
            wp_.cwl3a.gamma,   wp_.cwl3a.beta,
            wp_.cwl3a.escale,  wp_.cwl3a.ebias);

          Eigen::VectorXf w_out = cwl_forward(
            wstream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            wp_.cwl3w.weights, wp_.cwl3w.biases,
            wp_.cwl3w.gamma,   wp_.cwl3w.beta,
            wp_.cwl3w.escale,  wp_.cwl3w.ebias);

          Eigen::VectorXf parallel_out(2 * NEURON_N);
          parallel_out.head(NEURON_N) = a_out;
          parallel_out.tail(NEURON_N) = w_out;

          x = residual_combine(parallel_out, x);
        }

        // Output heads
        if (log_output)
        {
          // v3/v4: identity-activation heads emit log-space scalars; exp + clamp
          // and the (Î, ŵ) combine happen in output_combine_log.
          int head_in = N_CLASSES + NEURON_N;
          Eigen::VectorXf logA, logW;
          if (store_.version() == 4)
          {
            // 2-layer MLP head: Dense(head_in → HEAD_HIDDEN, gelu) → Dense(HEAD_HIDDEN → 1, id)
            logA = dense_linear(
              dense_gelu(astream_select(x), wp_.ahead_l1_w, wp_.ahead_l1_b, HEAD_HIDDEN, head_in),
              wp_.ahead_l2_w, wp_.ahead_l2_b, 1, HEAD_HIDDEN);
            logW = dense_linear(
              dense_gelu(wstream_select(x), wp_.whead_l1_w, wp_.whead_l1_b, HEAD_HIDDEN, head_in),
              wp_.whead_l2_w, wp_.whead_l2_b, 1, HEAD_HIDDEN);
          }
          else
          {
            // v3: single Dense head.
            logA = dense_linear(astream_select(x), wp_.ahead_w, wp_.ahead_b, 1, head_in);
            logW = dense_linear(wstream_select(x), wp_.whead_w, wp_.whead_b, 1, head_in);
          }

          return output_combine_log(logA, logW);
        }

        Eigen::VectorXf A_hat = dense_softplus(
          astream_select(x),
          wp_.ahead_w, wp_.ahead_b,
          1, N_CLASSES + NEURON_N);
        Eigen::VectorXf w_hat = dense_softplus(
          wstream_select(x),
          wp_.whead_w, wp_.whead_b,
          1, N_CLASSES + NEURON_N);

        return output_combine(A_hat, w_hat);
      }
  };
}
