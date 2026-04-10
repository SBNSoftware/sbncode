/**
 * WireModInfer.hh - WaireMLod interface for inference at SBN
 *
 * Loads weights for the Julia WaireMLod network from binary file
 * and runs model inference on hits to modify them to be more ``data-like``
 */

#include <Eigen/Dense>
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
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"

namespace sys
{
  // ─────────────────────────────────────────────────────────────────────────────
  // Model constants  (must match Constants.jl used in training)
  // ─────────────────────────────────────────────────────────────────────────────
  static constexpr int N_FEATURES   =  11;
  static constexpr int N_CLASSES    =   4;                    // classes
  static constexpr int NEURON_N     = 128;                    // neuron_n
  static constexpr int NORM_GROUPS  =   8;                    // normgroups
  static constexpr int N_FIXED      =  10;                    // fixed input scales
  static constexpr float SQRT2PI    =   2.5066282746310002f;  // sqrt(2*pi)
  static constexpr float GN_EPS     = 1e-5f;                  // GroupNorm epsilon
  static constexpr float FLT_EPS    =   1.1754944e-38f;       // eps(Float32)
  
  // ═════════════════════════════════════════════════════════════════════════════
  // Weight file reader
  // ═════════════════════════════════════════════════════════════════════════════
  
  /**
   * TensorStore - loads all named tensors from a WAIRE1 binary weight file
   * and provides access by name
   *
   * File format: see inference/README.md and Training.jl:export_weights_for_cpp from
   * WaireMLod julia project
   *
   * All values are little-endian
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
  
    private:
      std::map<std::string, std::vector<float>> tensors_;
  
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
  
        // check file is correct
        char label[6];
        f.read(label, 6);
        if (std::strncmp(label, "WAIRE1", 6) != 0)
          throw std::runtime_error("WaireMLod --- Not a WAIRE1 file: " + path);
  
        uint32_t version = read_le<uint32_t>(f);
        if (version != 1)
          throw std::runtime_error("WaireMLod --- Unsupported version: " + std::to_string(version));
  
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
    return x * 0.5f * (1.0f + std::erf(x * std::sqrt(2.0f)));
  }

  inline Eigen::VectorXf gelu(const Eigen::VectorXf& v)
  {
    return v.unaryExpr([](float x){ return gelu(x); });
  }
  
  inline float softplus(float x)
  {
    return std::log1p(std::exp(x));
  }
  
  inline Eigen::VectorXf softplus(const Eigen::VectorXf& v)
  {
    return v.unaryExpr([](float x){ return softplus(x); });
  }
  
  Eigen::VectorXf softmax(const Eigen::VectorXf& v)
  {
    Eigen::VectorXf shifted = v.array() - v.maxCoeff();
    Eigen::VectorXf ex      = shifted.array().exp();
    return ex / ex.sum();
  }
  
  // ═════════════════════════════════════════════════════════════════════════════
  // GroupNorm (single sample, C channels, G groups)
  // ═════════════════════════════════════════════════════════════════════════════
  
  /**
   * Apply GroupNorm to a single vector of `channels` divided into `groups` of size `channels / groups`
   *
   * For each group g:
   *   mean_g = mean(x[g*gs ... (g+1)*gs - 1])
   *   var_g  = variance(x[g*gs ... (g+1)*gs - 1])
   *   x[c]   = (x[c] - mean_g) / sqrt(var_g + eps) for c in group g
   *   x[c]   = gamma[c] * x[c] + beta[c]
   *
   * Matches Flux.GroupNorm applied to a (channels, 1) tensor
   *
   * gamma and beta are vectors of length `channels`
   */
  Eigen::VectorXf group_norm(
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
      // Normaize and affine
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
   * Dense layer forward pass: y = W * x + b
   *
   * W is stored as a Julia (out, in) column-major matrix
   * Eigen::Map<const MatrixXf>(ptr, out, in) reads it directly
   */
  Eigen::VectorXf dense_linear(
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

  Eigen::VectorXf dense_gelu(
    const Eigen::VectorXf& x,
    const float* W_ptr,
    const float* b_ptr,
    int out_dim,
    int in_dim)
  {
    return gelu(dense_linear(x, W_ptr, b_ptr, out_dim, in_dim));
  }

  Eigen::VectorXf dense_softplus(
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
   * Input scaling for forward pass
   *
   * Features 0...N_FIXED-1 are scaled by fixed amounts
   * Features N_FIXED.. are log1p(|s*x|) compressed
   */
  Eigen::VectorXf input_layer(
    const std::array<float, N_FEATURES>& raw,
    const float* scales_fixed, // length N_FIXED
    const float* scales_train) // length (N_FEATURES - N_FIXED)
  {
    Eigen::VectorXf out(N_FEATURES);
    for (int i = 0; i < N_FIXED; ++i)
      out[i] = raw[i] * scales_fixed[i];
    for (int i = N_FIXED; i < N_FEATURES; ++i)
      out[i] = std::log1p(std::fabs(raw[i] * scales_train[i - N_FIXED]));
    return out;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // Slicing operations
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * _astream_select: x[0 : N_CLASSES + NEURON_N - 1]
   * Drops the w-stream rows (indices N_CLASSES+NEURON_N ... end)
   * Input:  (N_CLASSES + 2*NEURON_N,) = (260,)
   * Output: (N_CLASSES + NEURON_N,)   = (132,)
   */
  Eigen::VectorXf astream_select(const Eigen::VectorXf& x) { return x.head(N_CLASSES + NEURON_N); }

  /**
   * _wstream_select: vcat(x[0:N_CLASSES-1], x[N_CLASSES+NEURON_N : end])
   * Drops the A-stream rows (indices N_CLASSES ... N_CLASSES+NEURON_N-1).
   * Input:  (N_CLASSES + 2*NEURON_N,) = (260,)
   * Output: (N_CLASSES + NEURON_N,)   = (132,)
   */
  Eigen::VectorXf wstream_select(const Eigen::VectorXf& x) {
    int total_out = N_CLASSES + NEURON_N;
    Eigen::VectorXf out(total_out);
    out.head(N_CLASSES)         = x.head(N_CLASSES);
    out.tail(NEURON_N)          = x.tail(NEURON_N);
    return out;
  }
 
  /**
   * _prepend_classes: vcat(x[0:N_CLASSES-1], cwl_out)
   * Used as the SkipConnection combiner after the first CWL.
   * x:       (N_CLASSES + NEURON_N,)   = (132,) classifier output
   * cwl_out: (2*NEURON_N,)             = (256,) CWL output
   * result:  (N_CLASSES + 2*NEURON_N,) = (260,)
   */
  Eigen::VectorXf prepend_classes(const Eigen::VectorXf& cwl_out, const Eigen::VectorXf& x) {
    int out_size = N_CLASSES + cwl_out.size();
    Eigen::VectorXf out(out_size);
    out.head(N_CLASSES) = x.head(N_CLASSES);
    out.tail(cwl_out.size()) = cwl_out;
    return out;
  }
 
  /**
   * _residual_combine: vcat(x[0:N_CLASSES-1], x[N_CLASSES:end] + residual)
   * Used as the SkipConnection combiner in the two residual blocks.
   * residual: (2*NEURON_N,)             = (256,) parallel block output
   * x:        (N_CLASSES + 2*NEURON_N,) = (260,) working tensor
   * result:   (260,)
   */
  Eigen::VectorXf residual_combine(const Eigen::VectorXf& residual, const Eigen::VectorXf& x) {
    Eigen::VectorXf out = x;
    out.tail(out.size() - N_CLASSES) += residual;
    return out;
  }
 
  /**
   * _output_combine: vcat(sqrt2pi * w_hat * A_hat, w_hat)
   * Converts the (A, w) pair to the (I, w) output pair.
   * A_hat, w_hat are scalars (size-1 vectors).
   * result: (2,) — row 0 = I, row 1 = w
   */
  Eigen::VectorXf output_combine(const Eigen::VectorXf& A_hat, const Eigen::VectorXf& w_hat) {
    Eigen::VectorXf out(2);
    out[0] = SQRT2PI * w_hat[0] * A_hat[0];
    out[1] = w_hat[0];
    return out;
  }

  // ═════════════════════════════════════════════════════════════════════════════
  // ClassWeightedLayer
  // ═════════════════════════════════════════════════════════════════════════════

  /**
   * ClassWeightedLayer forward pass for batch size = 1, flat weight layout
   */
  Eigen::VectorXf cwl_forward(
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

    // Inital pass
    Eigen::Map<const Eigen::MatrixXf> W(weights, flat_out, dense_in);
    Eigen::Map<const Eigen::VectorXf> b(biases, flat_out);
    Eigen::VectorXf act_flat = gelu(W * neurons + b);

    // Normize and mix
    Eigen::VectorXf mixed = Eigen::VectorXf::Zero(dense_out);

    for (int k = 0; k < n_classes; ++k)
    {
      int offset = k * dense_out;

      // Slice expert k's activations
      Eigen::VectorXf act_k = act_flat.segment(offset, dense_out);

      // Norm with shared gamma/beta
      Eigen::VectorXf normed_k = group_norm(act_k, dense_out, norm_groups, norm_gamma, norm_beta);

      // Per-expert affine
      Eigen::Map<const Eigen::VectorXf> es_k(expert_scale + offset, dense_out);
      Eigen::Map<const Eigen::VectorXf> eb_k(expert_bias  + offset, dense_out);
      normed_k = normed_k.cwiseProduct(es_k) + eb_k;

      // Weighted accumulation
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
   * Loads all weights from a WAIRE1 binary and infers data-like hit integrals and widths
   */
  class WaireMLod
  {
    public:
      explicit WaireMLod(const std::string& weights_path) : store_(weights_path) { validate_shapes(); }

      /**
       * Run inference on a single hit
       *
       * @param raw Array of 11 floats in the input feature order
       * @return {predicted_integral, predicted_width}
       */
      std::array<float, 2> infer(const std::array<float, N_FEATURES>& raw) const
      {
        Eigen::VectorXf x = run_classifier(raw);
        Eigen::VectorXf y = run_mixture_of_experts(x);
        return {y[0], y[1]};
      }

      /**
       * Construct vector of data-like hits from simulated hits
       */
      std::vector<recob::Hit> produceNew(
        const std::vector<recob::Hit> old_hits,
        cheat::BackTrackerService* back_tracker,
        detinfo::DetectorClocksData* det_clock,
        geo::WireReadoutGeom const* wire_geom) const;

    private:
      TensorStore store_;

      /**
       * Validate that all loaded tensors have the expected number of elements
       * Catches corrupted or mismatched weight files before they cause silent errors
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
        check("cls.weight",       N_CLASSES * NEURON_N);
        check("cls.bias",         N_CLASSES);

        // MoE: CWL1 (dense_out = 2*NEURON_N)
        int cwl1_out = 2 * NEURON_N;
        check("cwl1.weights",      cwl1_out * N_CLASSES * NEURON_N);
        check("cwl1.biases",       cwl1_out * N_CLASSES);
        check("cwl1.norm.gamma",   cwl1_out);
        check("cwl1.norm.beta",    cwl1_out);
        check("cwl1.expert_scale", cwl1_out * N_CLASSES);
        check("cwl1.expert_bias",  cwl1_out * N_CLASSES);

        // MoE: CWL2a/2w/3a/3w (dense_out = NEURON_N)
        for (const auto& prefix : {"cwl2a", "cwl2w", "cwl3a", "cwl3w"}) {
          std::string p(prefix);
          check(p + ".weights",      NEURON_N * N_CLASSES * NEURON_N);
          check(p + ".biases",       NEURON_N * N_CLASSES);
          check(p + ".norm.gamma",   NEURON_N);
          check(p + ".norm.beta",    NEURON_N);
          check(p + ".expert_scale", NEURON_N * N_CLASSES);
          check(p + ".expert_bias",  NEURON_N * N_CLASSES);
        }

        // Output heads
        int head_in = N_CLASSES + NEURON_N;
        check("ahead.weight", head_in);
        check("ahead.bias",   1);
        check("whead.weight", head_in);
        check("whead.bias",   1);
      }
  
      /**
       * Classifier subnetwork
       */
      Eigen::VectorXf run_classifier(const std::array<float, N_FEATURES>& raw) const
      {
        // Inputs
        Eigen::VectorXf x = input_layer(raw, 
          store_.data("il.scales_fixed"),
          store_.data("il.scales_train"));

        // Dense N_FEATURES -> NEURON_N, gelu
        x = dense_gelu(x, store_.data("d1.weight"), store_.data("d1.bias"), NEURON_N, N_FEATURES);

        // ResBlock 1: GroupNorm -> Dense, gelu -> residual
        {
          Eigen::VectorXf normed = group_norm(x, NEURON_N, NORM_GROUPS,
            store_.data("rb1.gn.gamma"), store_.data("rb1.gn.beta"));
          Eigen::VectorXf inner  = dense_gelu(normed,
            store_.data("rb1.d.weight"), store_.data("rb1.d.bias"),
            NEURON_N, NEURON_N);
          x += inner;
        }

        // ResBlock 2: GroupNorm -> Dense, gelu -> residual
        {
          Eigen::VectorXf normed = group_norm(x, NEURON_N, NORM_GROUPS,
            store_.data("rb2.gn.gamma"), store_.data("rb2.gn.beta"));
          Eigen::VectorXf inner  = dense_gelu(normed,
            store_.data("rb2.d.weight"), store_.data("rb2.d.bias"),
            NEURON_N, NEURON_N);
          x += inner;
        }

        // Classifier
        Eigen::VectorXf cls_scores = dense_linear(x,
          store_.data("cls.weight"), store_.data("cls.bias"),
          N_CLASSES, NEURON_N);
        cls_scores = softmax(cls_scores);
        Eigen::VectorXf out(N_CLASSES + NEURON_N);
        out.head(N_CLASSES) = cls_scores;
        out.tail(NEURON_N)  = x;

        return out;
      }

      /**
       * Mixure of Expert subnetwork
       */
      Eigen::VectorXf run_mixture_of_experts(const Eigen::VectorXf& clf_out) const
      {
        // Prep classifier output for A/w attention heads
        Eigen::VectorXf cwl1_out = cwl_forward(
          clf_out, N_CLASSES, NEURON_N, 2*NEURON_N, NORM_GROUPS,
          store_.data("cwl1.weights"),
          store_.data("cwl1.biases"),
          store_.data("cwl1.norm.gamma"),
          store_.data("cwl1.norm.beta"),
          store_.data("cwl1.expert_scale"),
          store_.data("cwl1.expert_bias"));
        Eigen::VectorXf x = prepend_classes(cwl1_out, clf_out);
        
        // ResBlock 3: Parallel A and w chains
        {
          Eigen::VectorXf a_out = cwl_forward(
            astream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            store_.data("cwl2a.weights"),
            store_.data("cwl2a.biases"),
            store_.data("cwl2a.norm.gamma"),
            store_.data("cwl2a.norm.beta"),
            store_.data("cwl2a.expert_scale"),
            store_.data("cwl2a.expert_bias"));

          Eigen::VectorXf w_out = cwl_forward(
            wstream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            store_.data("cwl2w.weights"),
            store_.data("cwl2w.biases"),
            store_.data("cwl2w.norm.gamma"),
            store_.data("cwl2w.norm.beta"),
            store_.data("cwl2w.expert_scale"),
            store_.data("cwl2w.expert_bias"));

          Eigen::VectorXf parallel_out(2 * NEURON_N);
          parallel_out.head(NEURON_N) = a_out;
          parallel_out.tail(NEURON_N) = w_out;

          x = residual_combine(parallel_out, x);
        }
        
        // ResBlock 4: Parallel A and w chains
        {
          Eigen::VectorXf a_out = cwl_forward(
            astream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            store_.data("cwl3a.weights"),
            store_.data("cwl3a.biases"),
            store_.data("cwl3a.norm.gamma"),
            store_.data("cwl3a.norm.beta"),
            store_.data("cwl3a.expert_scale"),
            store_.data("cwl3a.expert_bias"));

          Eigen::VectorXf w_out = cwl_forward(
            wstream_select(x), N_CLASSES, NEURON_N, NEURON_N, NORM_GROUPS,
            store_.data("cwl3w.weights"),
            store_.data("cwl3w.biases"),
            store_.data("cwl3w.norm.gamma"),
            store_.data("cwl3w.norm.beta"),
            store_.data("cwl3w.expert_scale"),
            store_.data("cwl3w.expert_bias"));

          Eigen::VectorXf parallel_out(2 * NEURON_N);
          parallel_out.head(NEURON_N) = a_out;
          parallel_out.tail(NEURON_N) = w_out;

          x = residual_combine(parallel_out, x);
        }

        // Output
        Eigen::VectorXf A_hat = dense_softplus(
          astream_select(x),
          store_.data("ahead.weight"), store_.data("ahead.bias"),
          1, N_CLASSES + NEURON_N);
        Eigen::VectorXf w_hat = dense_softplus(
          wstream_select(x),
          store_.data("whead.weight"), store_.data("whead.bias"),
          1, N_CLASSES + NEURON_N);

        return output_combine(A_hat, w_hat);
      }
  };
} 
