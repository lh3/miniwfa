//
//  wfa_lm.hpp
//
// MIT License
//
// Copyright (c) 2022 Jordan Eizenga
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef wfa_lm_hpp_included
#define wfa_lm_hpp_included

#include <cstdint>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>
#include <deque>
#include <sstream>
#include <limits>
#include <tuple>
#include <array>

#if defined(__SSE4_1__) || defined(__SSE4_2__)
#include <x86intrin.h>
#else
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse4.1.h"
#endif

namespace wfalm {

/*
 * Operation in a CIGAR string (as defined in SAM format specification)
 */
struct CIGAROp {
    CIGAROp(char op, uint32_t len) : op(op), len(len) {}
    CIGAROp() = default;
    
    char op;
    uint32_t len;
};

template<int NumPW> class WFAligner; // forward declaration

/*
 * A convenience specialization that uses a linear gap penalty.
 */
typedef WFAligner<0> LinearWFAligner;
/// Construct with WFA-style scoring parameters
inline LinearWFAligner make_linear_wfaligner(uint32_t mismatch, uint32_t gap);
/// Construct with Smith-Waterman-Gotoh-style scoring parameters
inline LinearWFAligner make_linear_wfaligner(uint32_t match, uint32_t mismatch, uint32_t gap);

/*
 * A convenience specialization that uses an affine gap penalty.
 */
typedef WFAligner<1> AffineWFAligner;
/// Construct with WFA-style scoring parameters
inline AffineWFAligner make_affine_wfaligner(uint32_t mismatch, uint32_t gap_extend, uint32_t gap_open);
/// Construct with Smith-Waterman-Gotoh-style scoring parameters
inline AffineWFAligner make_affine_wfaligner(uint32_t match, uint32_t mismatch, uint32_t gap_extend, uint32_t gap_open);

/*
 * A convenience specialization that uses a convex gap penalty.
 */
typedef WFAligner<2> ConvexWFAligner;
/// Construct with WFA-style scoring parameters
inline ConvexWFAligner make_convex_wfaligner(uint32_t mismatch, uint32_t gap_extend_1, uint32_t gap_open_1,
                                             uint32_t gap_extend_2, uint32_t gap_open_2);
/// Construct with Smith-Waterman-Gotoh-style scoring parameters
inline ConvexWFAligner make_convex_wfaligner(uint32_t match, uint32_t mismatch,
                                             uint32_t gap_extend_1, uint32_t gap_open_1,
                                             uint32_t gap_extend_2, uint32_t gap_open_2);

/*
 * Class that performs WFA with the standard, low-memory, or recusive variants. The aligner objects are
 * relatively lightweight and can be constructed for single-use with limited overhead. They are also
 * threadsafe, as long as the configurable parameters are not being altered.
 * The template integer controls how many  piecewise-affine segments there are in the gap penalty. As a
 * special case, a template parameter of 0 performs linear gap alignment.
 */
template<int NumPW>
class WFAligner {
    
public:
    
    /// Initialize with WFA-style score parameters. On opening an insertion or
    /// deletion, *both* the gap open and gap extend penalties are applied.
    /// Scores returned by alignment methods will also be WFA-style.
    /// One gap extend and gap open penalty are required for each of the piecewise
    /// affine segments. If the template parameter is set to 0, then one gap extend
    /// is still required, ang the aligner then performs linear gap alignment.
    WFAligner(uint32_t mismatch,
              std::array<uint32_t, std::max(1, NumPW)> gap_extend,
              std::array<uint32_t, NumPW> gap_open);
    
    /// Initialize with Smith-Waterman-Gotoh-style score parameters. On opening an
    /// insertion or deletion, *both* the gap open and gap extend penalties are applied.
    /// Scores returned by alignment methods will also be Smith-Waterman-Gotoh-style.
    /// One gap extend and gap open penalty are required for each of the piecewise
    /// affine segments. If the template parameter is set to 0, then one gap extend
    /// is still required, ang the aligner then performs linear gap alignment.
    WFAligner(uint32_t match, uint32_t mismatch,
              std::array<uint32_t, std::max(1, NumPW)> gap_extend,
              std::array<uint32_t, NumPW> gap_open);
    
    /// Default constructor
    WFAligner();

    /*****************************
     *  Configurable parameters  *
     *****************************/
    
    /// If set to a number >= 0, prunes diagonals that fall this many anti-diagonals
    /// behind the furthest-reaching diagonal. This can lead to suboptimal alignments,
    /// but it also can increase speed and reduce memory use. If set to a number < 0,
    /// no pruning is performed.
    int32_t lagging_diagonal_prune = -1;
    
    /***********************
     *  Alignment methods  *
     ***********************/
    
    /// Globally align two sequences using the fastest algorithm that will remain
    /// constrained to a given maximum memory.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///   max_mem      The target maximum memory use in bytes
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_adaptive(const char* seq1, size_t len1,
                             const char* seq2, size_t len2,
                             uint64_t max_mem) const;
    
    /// Globally align two sequences in O(s log s) memory and O(sN log s) time using the recursive WFA algorithm.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_recursive(const char* seq1, size_t len1,
                              const char* seq2, size_t len2) const;
    
    /// Globally align two sequences in O(s^3/2) memory and O(sN) time using the low-memory WFA algorithm.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_low_mem(const char* seq1, size_t len1,
                            const char* seq2, size_t len2) const;
    
    /// Globally align two sequences in O(s^2) memory and O(sN) time using the standard WFA algorithm.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align(const char* seq1, size_t len1,
                    const char* seq2, size_t len2) const;
    
    /// Locally align two sequences from a seed using the adaptive WFA as an
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   max_mem          The target maximum memory use in bytes
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_adaptive(const char* seq1, size_t len1,
                                   const char* seq2, size_t len2,
                                   uint64_t max_mem,
                                   size_t anchor_begin_1, size_t anchor_end_1,
                                   size_t anchor_begin_2, size_t anchor_end_2,
                                   bool anchor_is_match = true) const;
    
    /// Locally align two sequences from a seed using the recursive WFA as an
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_recursive(const char* seq1, size_t len1,
                                    const char* seq2, size_t len2,
                                    size_t anchor_begin_1, size_t anchor_end_1,
                                    size_t anchor_begin_2, size_t anchor_end_2,
                                    bool anchor_is_match = true) const;
    
    /// Locally align two sequences from a seed using the low memory WFA as an
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_low_mem(const char* seq1, size_t len1,
                                  const char* seq2, size_t len2,
                                  size_t anchor_begin_1, size_t anchor_end_1,
                                  size_t anchor_begin_2, size_t anchor_end_2,
                                  bool anchor_is_match = true) const;
    
    /// Locally align two sequences from a seed using the standard WFA as an
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local(const char* seq1, size_t len1,
                          const char* seq2, size_t len2,
                          size_t anchor_begin_1, size_t anchor_end_1,
                          size_t anchor_begin_2, size_t anchor_end_2,
                          bool anchor_is_match = true) const;
    
    
    
    
    
    
    
    
    /*******************************************************************
     *******************************************************************
     ***             Only internal functions below here              ***
     *******************************************************************
     *******************************************************************/
    
    
    
    
    
    
    
    
    
protected:

    static const int StdMem = 0;
    static const int LowMem = 1;
    static const int Recursive = 2;
    static const int AdaptiveMem = 3;
    
    // forward declaration
    template<typename IntType> struct Wavefront;
    
    /*
     * Allows the WF bank to be used in wf_next
     */
    template<typename IntType>
    struct WFBankAdapter {
        WFBankAdapter(const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank) : wf_bank(wf_bank) {
            
        }
        
        inline const Wavefront<IntType>& operator[](size_t i) const {
            return wf_bank[i].second;
        };
        
        inline size_t size() const {
            return wf_bank.size();
        }
        
        const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank;
    };
    
    /*
     * wrappers to perform WFA alignment algorithms for dependency injection
     */
    template<typename MatchFunc, typename StringType>
    struct StandardGlobalWFA {
        StandardGlobalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner<NumPW>::StdMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                       std::numeric_limits<uint64_t>::max());
        }
        const WFAligner<NumPW>* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct StandardSemilocalWFA {
        StandardSemilocalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner<NumPW>::StdMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                      std::numeric_limits<uint64_t>::max());
        }
        const WFAligner<NumPW>* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct LowMemGlobalWFA {
        LowMemGlobalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner<NumPW>::LowMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                       std::numeric_limits<uint64_t>::max());
        }
        const WFAligner<NumPW>* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct LowMemSemilocalWFA {
        LowMemSemilocalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner<NumPW>::LowMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                      std::numeric_limits<uint64_t>::max());
        }
        const WFAligner<NumPW>* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct RecursiveGlobalWFA {
        RecursiveGlobalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner<NumPW>::Recursive, MatchFunc, StringType>(seq1, seq2,
                                                                                                          std::numeric_limits<uint64_t>::max());
        }
        const WFAligner<NumPW>* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct RecursiveSemilocalWFA {
        RecursiveSemilocalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner<NumPW>::Recursive, MatchFunc, StringType>(seq1, seq2,
                                                                                                         std::numeric_limits<uint64_t>::max());
        }
        const WFAligner<NumPW>* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct AdaptiveGlobalWFA {
        AdaptiveGlobalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner), max_mem(max_mem) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner<NumPW>::AdaptiveMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                            max_mem);
        }
        const WFAligner<NumPW>* aligner = nullptr;
        uint64_t max_mem = 0;
    };
    template<typename MatchFunc, typename StringType>
    struct AdaptiveSemilocalWFA {
        AdaptiveSemilocalWFA(const WFAligner<NumPW>* aligner, uint64_t max_mem) : aligner(aligner), max_mem(max_mem) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner<NumPW>::AdaptiveMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                           max_mem);
        }
        const WFAligner<NumPW>* aligner = nullptr;
        uint64_t max_mem = 0;
    };
    
    // give the wrapper access to the dispatch function
    template<class MatchFunc, class StringType> friend class StandardGlobalWFA;
    template<class MatchFunc, class StringType> friend class StandardSemilocalWFA;
    template<class MatchFunc, class StringType> friend class LowMemGlobalWFA;
    template<class MatchFunc, class StringType> friend class LowMemSemilocalWFA;
    template<class MatchFunc, class StringType> friend class RecursiveGlobalWFA;
    template<class MatchFunc, class StringType> friend class RecursiveSemilocalWFA;
    template<class MatchFunc, class StringType> friend class AdaptiveGlobalWFA;
    template<class MatchFunc, class StringType> friend class AdaptiveSemilocalWFA;
    
    inline void init(uint32_t mismatch,
                     std::array<uint32_t, std::max(1, NumPW)> gap_extend,
                     std::array<uint32_t, NumPW> gap_open);
    
    // central routine that back-ends the user-facing interface
    template<bool Local, int Mem, typename MatchFunc, typename StringType>
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_dispatch(const StringType& seq1, const StringType& seq2, uint64_t max_mem) const;
    
    // greedy take matches in WFA
    template<typename StringType, typename MatchFunc, typename IntType>
    inline void wavefront_extend(const StringType& seq1, const StringType& seq2,
                                 Wavefront<IntType>& wf, const MatchFunc& match_func) const;
    
    // increase score and initialize row in WFA
    template<typename IntType, typename WFVector, typename StringType>
    inline Wavefront<IntType> wavefront_next(const StringType& seq1, const StringType& seq2,
                                             const WFVector& wfs) const;
    
    // prune lagging diagonals in WFA
    template<bool Local, typename IntType>
    inline void wavefront_prune(Wavefront<IntType>& wf) const;
    
    // check if reached complete alignment
    template<typename IntType>
    inline bool wavefront_reached(const Wavefront<IntType>& wf, int32_t diag, int32_t anti_diag) const;
    
    // traceback through a complete DP structure
    template <typename StringType, typename IntType>
    inline std::vector<CIGAROp> wavefront_traceback(const StringType& seq1, const StringType& seq2,
                                                    const std::vector<Wavefront<IntType>>& wfs,
                                                    int64_t s, int64_t d) const;
    
    // traceback as far as possible through a possibly incomplete DP structure
    // creates the CIGAR in reverse order,
    // appends new CIGAR operations to the CIGAR string that is passed in
    template<typename WFVector, typename StringType>
    void wavefront_traceback_internal(const StringType& seq1, const StringType& seq2,
                                      const WFVector& wfs, int64_t& d, int64_t& lead_matches,
                                      int& mat, int64_t& s, std::vector<CIGAROp>& cigar) const;
    
    // traceback through the low memory algorithms DP structure, with recomputation
    template <bool Local, typename StringType, typename MatchFunc, typename IntType>
    std::vector<CIGAROp> wavefront_traceback_low_mem(const StringType& seq1, const StringType& seq2,
                                                     std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                                     int64_t s, int64_t d, const MatchFunc& match_func) const;
    
    // do classic WFA
    template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_core(const StringType& seq1, const StringType& seq2,
                         const MatchFunc& match_func, uint64_t max_mem) const;
    
    // do low memory WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_low_mem_core(const StringType& seq1, const StringType& seq2,
                                 const MatchFunc& match_func) const;
    
    // internal routine for low memory WFA after initialization or fall-back
    template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_low_mem_core_internal(const StringType& seq1, const StringType& seq2,
                                          const MatchFunc& match_func, uint64_t max_mem,
                                          int64_t epoch_len, int64_t epoch_end, int64_t sample_rate,
                                          int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                          std::deque<Wavefront<IntType>>& wf_buffer, size_t s,
                                          std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                          uint64_t curr_mem_buffer, uint64_t curr_mem_bank,
                                          uint64_t curr_mem_block, uint64_t max_mem_block) const;
    
    // do recursive WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_recursive_core(const StringType& seq1, const StringType& seq2,
                                   const MatchFunc& match_func) const;
    
    // internal routine for recursive WFA after initialization or fall-back
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_recursive_core_internal(const StringType& seq1, const StringType& seq2,
                                            const MatchFunc& match_func,
                                            int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                            std::vector<Wavefront<IntType>>& opt_stripe,
                                            std::vector<Wavefront<IntType>>& first_stripe,
                                            std::deque<Wavefront<IntType>>& wf_buffer, size_t s) const;
    
    // recursive calls of recursive WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc, typename WFVector1, typename WFVector2>
    void wavefront_align_recursive_internal(const StringType& seq1, const StringType& seq2,
                                            int64_t lower_stripe_num, WFVector1& lower_stripe,
                                            int64_t upper_stripe_num, WFVector2& upper_stripe,
                                            int64_t& traceback_s, int64_t& traceback_d,
                                            int& traceback_mat, int64_t& lead_matches,
                                            std::vector<CIGAROp>& cigar,
                                            const MatchFunc& match_func) const;
    
    // give up on a standard WFA and switch to doing low memory WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    fall_back_to_low_mem(const StringType& seq1, const StringType& seq2,
                         const MatchFunc& match_func, uint64_t max_mem,
                         int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                         std::vector<Wavefront<IntType>>& wfs) const;
    
    // give up on low memory WFA and switch to doing recursive WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    fall_back_to_recursive(const StringType& seq1, const StringType& seq2,
                           const MatchFunc& match_func,
                           int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                           std::deque<Wavefront<IntType>>& wf_buffer, size_t s,
                           std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank) const;
    
    // check score that are SWG-locally optimal
    template<typename StringType, typename IntType>
    inline void find_local_opt(const StringType& seq1, const StringType& seq2,
                               int64_t s, const Wavefront<IntType>& wf,
                               int64_t& opt, int64_t& opt_diag, int64_t& opt_s, size_t& max_s) const;
    
    // do anchored local alignment
    template<typename PrefSemilocalWFA, typename AnchorGlobalWFA, typename SuffSemilocalWFA>
    std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_core(const char* seq1, size_t len1,
                               const char* seq2, size_t len2,
                               uint64_t max_mem,
                               size_t anchor_begin_1, size_t anchor_end_1,
                               size_t anchor_begin_2, size_t anchor_end_2,
                               bool anchor_is_match) const;
    
    // greatest commmon denominator
    uint32_t gcd(uint32_t a, uint32_t b) const;
    
    // convert a WFA score into a SWG score
    inline int32_t convert_score(size_t len1, size_t len2, int32_t score) const;
    
    // merge adjacent CIGAR ops that are the same type
    inline void coalesce_cigar(std::vector<CIGAROp>& cigar) const;
    
    // get the length of the aligned sequence on both sequences
    inline std::pair<size_t, size_t> cigar_base_length(const std::vector<CIGAROp>& cigar) const;
    
    // WFA-style parameters
    uint32_t mismatch = 1;
    std::array<uint32_t, std::max(1, NumPW)> gap_extend;
    std::array<uint32_t, NumPW> gap_open;
    
    // scores are reduced by a constant factor if they have a non-trivial common
    // denominator to reduce run time
    uint32_t factor = 1;
    
    // SWG match score, or 0 if not using SWG scoring
    uint32_t match = 0;
    
    // how many subsequent wavefronts are required to do DP with these scoring parameters
    int64_t stripe_width = 1;
    
private:
    
    /* debug functions*/
    template<typename IntType>
    inline void print_wf(const Wavefront<IntType>& wf, int diag_begin, int diag_end) const;
    template<typename WFVector>
    void print_wfs(const WFVector& wfs) const;
    template<typename IntType>
    inline void print_wfs_low_mem(const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                  const std::deque<Wavefront<IntType>>& wf_buffer, int s) const;
    template<typename IntType>
    inline void print_wfs_tb(const std::deque<Wavefront<IntType>>& traceback_block,
                             const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                             int i) const;
    template <typename StringType, typename WFVector>
    void wavefront_viz(const StringType& seq1, const StringType& seq2,
                       const WFVector& wfs) const;
    
    static_assert(NumPW >= 0, "Must choose a non-negative number of piecewise linear penalities");
};

/*
 * Adapter to avoid string copying for when aligning subintervals
 */
class StringView {
public:
    StringView(const char* seq, size_t i, size_t len) : seq(seq + i), len(len)
    {
        
    }
    
    inline const char& operator[](size_t idx) const {
        return seq[idx];
    }
    
    inline size_t size() const {
        return len;
    }
    
private:
    
    const char* seq;
    size_t len;
};

/*
 * Adapter to avoid string copying for when aligning subintervals
 * that also reverses the string
 */
class RevStringView {
public:
    // note: interval is provided in forward direction
    RevStringView(const char* seq, size_t i, size_t len) : seql(seq + i + len - 1), len(len)
    {
        
    }
    
    inline const char& operator[](size_t idx) const {
        return *(seql - idx);
    }
    
    inline size_t size() const {
        return len;
    }
    
private:
    
    // last char in the seq
    const char* seql;
    size_t len;
};

/*
 * templated SIMD operations by integer width
 */
template<typename IntType>
struct SIMD {
    static inline int vec_size() {
        return sizeof(__m128i) / sizeof(IntType);
    }
    static inline int round_up(int l) {
        return (l + vec_size() - 1) / vec_size();
    }
    static inline int round_down(int l) {
        return l / vec_size();
    }
    static inline bool is_aligned(__m128i* ptr) {
        return !(((uintptr_t) ptr) & (sizeof(__m128i) - 1));
    }
    static inline __m128i load_unaligned(const __m128i* ptr) {
        return _mm_loadu_si128(ptr);
    }
    static inline __m128i broadcast(IntType i);
    static inline __m128i min_inf() {
        return broadcast(std::numeric_limits<IntType>::min());
    }
    // O, 1, ..., L-1
    static inline __m128i range();
    static inline __m128i mask_and(__m128i a, __m128i b) {
        return _mm_and_si128(a, b);
    }
    static inline __m128i mask_not(__m128i a) {
        // there's not bitwise not, so i have to simulate it with xor
        return _mm_xor_si128(a, broadcast(-1));
    }
    static inline __m128i eq(__m128i a, __m128i b);
    static inline __m128i less(__m128i a, __m128i b);
    static inline __m128i leq(__m128i a, __m128i b) {
        return mask_not(less(b, a));
    }
    static inline __m128i max(__m128i a, __m128i b);
    static inline __m128i add(__m128i a, __m128i b);
    static inline __m128i subtract(__m128i a, __m128i b);
    //static inline __m128i div2(__m128i a);
    // t ? a : b
    // note: assumes all 1 bits in condition so that 8-bit wors for all widths
    static inline __m128i ifelse(__m128i t, __m128i a, __m128i b) {
        return _mm_blendv_epi8(b, a, t);
    }
};

// a match-computing function based on direct character comparison
class CompareMatchFunc {
    
protected:
    
    template<bool Reversed, class StringType>
    inline void load(const StringType& seq1, const StringType& seq2,
                     size_t i, size_t j, __m128i& vec1, __m128i& vec2) const {
        if (Reversed) {
            static const __m128i reverser = _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
            vec1 = _mm_shuffle_epi8(SIMD<int8_t>::load_unaligned((__m128i const*)(&seq1[i] - sizeof(__m128i) + 1)), reverser);
            vec2 = _mm_shuffle_epi8(SIMD<int8_t>::load_unaligned((__m128i const*)(&seq2[j] - sizeof(__m128i) + 1)), reverser);
            
        }
        else {
            vec1 = SIMD<int8_t>::load_unaligned((__m128i const*)&seq1[i]);
            vec2 = SIMD<int8_t>::load_unaligned((__m128i const*)&seq2[j]);
        }
    }
    
    // sub-routine of the internal function that computes the length of the prefix match
    // in the first 8 bytes of an equality mask vector
    inline size_t partial_vec_match(const __m128i& eq) const {
        
        // find the mismatch among the first 8 bytes if there is one
        __m128i pos = _mm_minpos_epu16(_mm_cvtepi8_epi16(eq));
        return _mm_extract_epi16(pos, 0) == 0 ? _mm_extract_epi16(pos, 1) : 8;
        
    }
    
    // a template-controlled generic function for the specializations to call
    template<bool Reversed, class StringType>
    inline size_t match_func_internal(const StringType& seq1, const StringType& seq2, size_t i, size_t j) const {
        // remember where we start from
        size_t init = i;
        
        // TODO: is it safe to make this execute within 8 bytes from the end since the first
        // comparison is only of the first 8 bases anyway?
        
        // use vectorized operations to find matches if we can
        bool found_mismatch = false;
        if (i + SIMD<int8_t>::vec_size() <= seq1.size() &&
            j + SIMD<int8_t>::vec_size() <= seq2.size()) {
            
            // try to find a partial match of the first 8 bases
            __m128i vec1, vec2;
            load<Reversed, StringType>(seq1, seq2, i, j, vec1, vec2);
            
            __m128i eq = SIMD<int8_t>::eq(vec1, vec2);
            
            size_t match_len = partial_vec_match(eq);
            i += match_len;
            j += match_len;
            
            found_mismatch = (match_len < 8);
            if (!found_mismatch) {
                // we matched the first 8 bases, so now we transition into a cheaper matching of an
                // entire vector to move quickly in instances of very long matches
                
                // TODO: i could maybe load aligned after the first iteration, but you can't guarantee
                // alignment of both strings...
                
                while (i + SIMD<int8_t>::vec_size() <= seq1.size() &&
                       j + SIMD<int8_t>::vec_size() <= seq2.size()) {
                    
                    load<Reversed, StringType>(seq1, seq2, i, j, vec1, vec2);
                    eq = SIMD<int8_t>::eq(vec1, vec2);
                    
                    if (_mm_test_all_ones(eq)) {
                        // we matched the entire vector
                        i += 16;
                        j += 16;
                    }
                    else {
                        // there's a mismatch in the vector somewhere
                        match_len = partial_vec_match(eq);
                        i += match_len;
                        j += match_len;
                        
                        if (match_len == 8) {
                            // the mismatch is after the first 8 bases, shift over and find it
                            eq = _mm_alignr_epi8(SIMD<int8_t>::broadcast(0), eq, 8);
                            match_len = partial_vec_match(eq);
                            i += match_len;
                            j += match_len;
                        }
                        found_mismatch = true;
                        break;
                    }
                }
            }
        }
        
        if (!found_mismatch) {
            // we got too close to the end of the sequences to do vectorized operations
            // without finding a mismatch, finish out with
            while (i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]) {
                ++i;
                ++j;
            }
        }
        
        return i - init;
    }
};

/**
 * Specialization for forward matching
 */
class FwdCompareMatchFunc : public CompareMatchFunc {
public:
    
    FwdCompareMatchFunc(const StringView& seq1, const StringView& seq2)
        : seq1(seq1), seq2(seq2)
    {
        
    }
    
    // leave the actual implementation to the template specializations
    inline size_t operator()(size_t i, size_t j) const {
        return match_func_internal<false, StringView>(seq1, seq2, i, j);
    }
    
protected:
    
    const StringView& seq1;
    const StringView& seq2;
};

/**
 * Specialization for reverse matching
 */
class RevCompareMatchFunc : public CompareMatchFunc {
public:
    
    RevCompareMatchFunc(const RevStringView& seq1, const RevStringView& seq2)
        : seq1(seq1), seq2(seq2)
    {
        
    }
    
    // leave the actual implementation to the template specializations
    inline size_t operator()(size_t i, size_t j) const {
        return match_func_internal<true, RevStringView>(seq1, seq2, i, j);
    }
    
protected:
    
    const RevStringView& seq1;
    const RevStringView& seq2;
};

template<int NumPW>
template<typename IntType>
struct WFAligner<NumPW>::Wavefront {
    
private:
    
    inline void init_arrays() {
        // we pad the array with a full vector on each end to handle boundary conditions
        interval = SIMD<IntType>::round_up(len) + 2;
        int alloc_size = (2 * NumPW + 1) * interval * sizeof(__m128i);
        if (posix_memalign((void**)&alloced, sizeof(__m128i), alloc_size)) {
            throw std::runtime_error(std::string("error:[WFAligner] could not perform aligned allocation of size "
                                                 + std::to_string(alloc_size)).c_str());
        }
        // M is at the middle of the array
        M = alloced + interval * NumPW + 1;
    }
public:
    
    
    Wavefront() {}
    Wavefront(int32_t diag_begin, int32_t diag_end)
        : diag_begin(diag_begin), len(diag_end - diag_begin)
    {

        init_arrays();
        for (int32_t i = 0, n = interval * (2 * NumPW + 1); i < n; ++i) {
            // TODO: i might be able to get away with only setting the
            // shoulders now that the ifelse condition in DP assigns -infinity
            alloced[i] = SIMD<IntType>::min_inf();
        }
    }
    
    inline Wavefront& operator=(const Wavefront& other) noexcept {
        if (this != &other) {
            diag_begin = other.diag_begin;
            len = other.len;
            init_arrays();
            
            // set the boundaries
            for (int i = -NumPW; i <= NumPW; ++i) {
                auto A = M + interval * i;
                *(A - 1) = SIMD<IntType>::min_inf();
                *(A + interval - 2) = SIMD<IntType>::min_inf();
            }
            
            // copy the interior of the array
            auto vec_n = SIMD<IntType>::round_up(len);
            if (SIMD<IntType>::is_aligned(other.M)) {
                // the arrays are aligned
                for (int i = -NumPW; i <= NumPW; ++i) {
                    auto A = M + interval * i;
                    auto B = other.M + other.interval * i;
                    for (int32_t j = 0; j < vec_n; ++j) {
                        A[i] = B[i];
                    }
                }
            }
            else {
                // the other arrays are not aligned
                for (int i = -NumPW; i <= NumPW; ++i) {
                    auto A = M + interval * i;
                    auto B = other.M + other.interval * i;
                    for (int32_t j = 0; j < vec_n; ++j) {
                        A[i] = SIMD<IntType>::load_unaligned(B + i);
                    }
                }
            }
        }
        return *this;
    }
    
    inline Wavefront& operator=(Wavefront&& other) noexcept {
        if (this != &other) {
            diag_begin = other.diag_begin;
            len = other.len;
            interval = other.interval;
            if (alloced) {
                free(alloced);
            }
            alloced = other.alloced;
            M = other.M;
            
            other.alloced = other.M = nullptr;
            other.len = 0;
            other.interval = 0;
        }
        return *this;
    }
    
    Wavefront(const Wavefront& other) noexcept {
        *this = other;
    }
    
    Wavefront(Wavefront&& other) noexcept {
        *this = std::move(other);
    }
    
    ~Wavefront() {
        free(alloced);
    }
    
    // mat == 0 for M, < 0 for Is, > 0 for Ds
    inline __m128i* vector_at(int32_t diag, int mat) {
        return (__m128i*)(((IntType*) (M + mat * interval)) + (diag - diag_begin));
    }
    inline const __m128i* vector_at(int32_t diag, int mat) const {
        return (__m128i*)(((IntType*) (M + mat * interval)) + (diag - diag_begin));
    }
    inline IntType& int_at(int32_t diag, int mat) {
        return *(((IntType*) (M + mat * interval) ) + (diag - diag_begin));
    }
    inline const IntType& int_at(int32_t diag, int mat) const {
        return *(((IntType*) (M + mat * interval) ) + (diag - diag_begin));
    }
    
    // warning: assumes that this interval is contained within the current one
    // and that the vectors are currently aligned (i.e. should only be called
    // once during wf_prune)
    // it also might not work if pruning to size 0, but we don't need that
    inline void resize(int32_t _diag_begin, int32_t _size) {
        if (_diag_begin + _size != diag_begin + len) {
            // we're pruning the right side
            
            // blend the right transition
            int32_t k = SIMD<IntType>::round_down(_diag_begin + _size - diag_begin - 1);
            __m128i mask = SIMD<IntType>::less(SIMD<IntType>::add(SIMD<IntType>::broadcast(diag_begin +
                                                                                           k * SIMD<IntType>::vec_size()),
                                                                  SIMD<IntType>::range()),
                                               SIMD<IntType>::broadcast(_diag_begin + _size));
            
            for (int i = -NumPW; i <= NumPW; ++i) {
                auto A = M + interval * i;
                A[k] = SIMD<IntType>::ifelse(mask, A[k], SIMD<IntType>::min_inf());
            }
            
            // copy full -infinity vectors
            ++k;
            for (int32_t n = SIMD<IntType>::round_up(diag_begin + len); k < n; ++k) {
                for (int i = -NumPW; i <= NumPW; ++i) {
                    auto A = M + interval * i;
                    A[k] = SIMD<IntType>::min_inf();
                }
            }
        }
        
        if (_diag_begin != diag_begin) {
            // we're pruning the left side
            
            // copy full -infinity vectors
            int32_t k = 0;
            for (int32_t n = SIMD<IntType>::round_down(_diag_begin - diag_begin); k < n; ++k) {
                for (int i = -NumPW; i <= NumPW; ++i) {
                    auto A = M + interval * i;
                    A[k] = SIMD<IntType>::min_inf();
                }
            }
            // blend the left transition
            __m128i mask = SIMD<IntType>::less(SIMD<IntType>::add(SIMD<IntType>::broadcast(diag_begin +
                                                                                           k * SIMD<IntType>::vec_size()),
                                                                  SIMD<IntType>::range()),
                                               SIMD<IntType>::broadcast(_diag_begin));
            for (int i = -NumPW; i <= NumPW; ++i) {
                auto A = M + interval * i;
                A[k] = SIMD<IntType>::ifelse(mask, SIMD<IntType>::min_inf(), A[k]);
            }
            
            // update the pointers
            int32_t diff = _diag_begin - diag_begin;
            M = (__m128i*)(((IntType*) M) + diff);
            diag_begin = _diag_begin;
        }
        
        len = _size;
    }
    
    inline int64_t size() const {
        return len;
    }
    
    inline bool empty() const {
        return len == 0;
    }
    
    inline uint64_t bytes() {
        return sizeof(__m128i) * (2 * NumPW + 1) * interval;
    }
    
    inline void print(std::ostream& strm) {
        for (int i = 0, j = 0; i < 2 * NumPW + 1; ++i) {
            for (int k = 0; k < interval; ++k, ++j) {
                if (k) {
                    strm << "\t";
                }
                IntType* v = (IntType*) &alloced[j];
                strm << "(";
                for (int l = 0; l < SIMD<IntType>::vec_size(); ++l) {
                    if (l) {
                        strm << "\t";
                    }
                    strm << (int) v[l];
                }
                strm << ")";
            }
            strm << std::endl;
        }
    }
    
    __m128i* alloced = nullptr;
    __m128i* M = nullptr;
    
    int32_t diag_begin = 0;
    // in number of integers
    int32_t len = 0;
    // the interval between matrix arrays in number of vectors
    int32_t interval = 0;
};

inline LinearWFAligner make_linear_wfaligner(uint32_t mismatch, uint32_t gap) {
    return WFAligner<0>(mismatch, std::array<uint32_t, 1>{gap}, std::array<uint32_t, 0>());
}

inline LinearWFAligner make_linear_wfaligner(uint32_t match, uint32_t mismatch, uint32_t gap) {
    return WFAligner<0>(match, mismatch, std::array<uint32_t, 1>{gap}, std::array<uint32_t, 0>());
}


inline AffineWFAligner make_affine_wfaligner(uint32_t mismatch, uint32_t gap_extend, uint32_t gap_open) {
    return WFAligner<1>(mismatch, std::array<uint32_t, 1>{gap_extend}, std::array<uint32_t, 1>{gap_open});
}

inline AffineWFAligner make_affine_wfaligner(uint32_t match, uint32_t mismatch, uint32_t gap_extend, uint32_t gap_open) {
    return WFAligner<1>(match, mismatch, std::array<uint32_t, 1>{gap_extend}, std::array<uint32_t, 1>{gap_open});
}

inline ConvexWFAligner make_convex_wfaligner(uint32_t mismatch, uint32_t gap_extend_1, uint32_t gap_open_1,
                                             uint32_t gap_extend_2, uint32_t gap_open_2) {
    return WFAligner<2>(mismatch, std::array<uint32_t, 2>{gap_extend_1, gap_extend_2},
                        std::array<uint32_t, 2>{gap_open_1, gap_open_2});
}

inline ConvexWFAligner make_convex_wfaligner(uint32_t match, uint32_t mismatch,
                                             uint32_t gap_extend_1, uint32_t gap_open_1,
                                             uint32_t gap_extend_2, uint32_t gap_open_2) {
    return WFAligner<2>(match, mismatch, std::array<uint32_t, 2>{gap_extend_1, gap_extend_2},
                        std::array<uint32_t, 2>{gap_open_1, gap_open_2});
}

template<int NumPW>
inline void WFAligner<NumPW>::init(uint32_t _mismatch,
                                   std::array<uint32_t, std::max(1, NumPW)> _gap_extend,
                                   std::array<uint32_t, NumPW> _gap_open) {
    // reduce by common factor to accelerate the alignment
    
    if (_mismatch == 0) {
        std::stringstream strm;
        strm << "error:[WFAligner] mismatch penalty must be non-zero: " << _mismatch << std::endl;
        throw std::runtime_error(strm.str());
    }
    factor = _mismatch;
    for (size_t i = 0; i < _gap_extend.size(); ++i) {
        if (_gap_extend[i] == 0) {
            std::stringstream strm;
            strm << "error:[WFAligner] gap extend penalties must be non-zero: " << _gap_extend[i] << std::endl;
            throw std::runtime_error(strm.str());
        }
        factor = gcd(factor, _gap_extend[i]);
    }
    for (size_t i = 0; i < _gap_open.size(); ++i) {
        if (_gap_open[i] != 0) {
            factor = gcd(factor, _gap_open[i]);
        }
    }
    
    mismatch = _mismatch / factor;
    stripe_width = mismatch;
    
    for (size_t i = 0; i < _gap_open.size(); ++i) {
        gap_open[i] = _gap_open[i] / factor;
    }
    for (size_t i = 0; i < _gap_extend.size(); ++i) {
        gap_extend[i] = _gap_extend[i] / factor;
        
        if (NumPW != 0) {
            stripe_width = std::max<int64_t>(stripe_width, gap_open[i] + gap_extend[i]);
        }
        else {
            stripe_width = std::max<int64_t>(stripe_width, gap_extend[i]);
        }
    }
}

// default
template<int NumPW>
inline WFAligner<NumPW>::WFAligner() {}

template<int NumPW>
inline WFAligner<NumPW>::WFAligner(uint32_t _mismatch,
                                   std::array<uint32_t, std::max(1, NumPW)> _gap_extend,
                                   std::array<uint32_t, NumPW> _gap_open) {
    
    init(_mismatch, _gap_extend, _gap_open);
};

template<int NumPW>
inline WFAligner<NumPW>::WFAligner(uint32_t _match, uint32_t _mismatch,
                                   std::array<uint32_t, std::max(1, NumPW)> _gap_extend,
                                   std::array<uint32_t, NumPW> _gap_open)
{
    if (_match == 0) {
        std::stringstream strm;
        strm << "error:[WFAligner] match score must be an integer > 0 for Smith-Waterman-Gotoh-style score parameters: " << match << "\n";
        throw std::runtime_error(strm.str());
    }
    
    // convert to equivalent WFA-style params
    uint32_t adj_mismatch = 2 * (_match + _mismatch);
    std::array<uint32_t, std::max(1, NumPW)> adj_gap_extend;
    std::array<uint32_t, NumPW> adj_gap_open;
    for (size_t i = 0; i < _gap_extend.size(); ++i) {
        adj_gap_extend[i] = 2 * _gap_extend[i] + _match;
    }
    for (size_t i = 0; i < _gap_open.size(); ++i) {
        adj_gap_open[i] = 2 * _gap_open[i];
    }
    
    // init the WFA-style params
    init(adj_mismatch, adj_gap_extend, adj_gap_open);
        
    // and remember the match score
    match = _match;
}

// 8 bit integers
template<> inline __m128i SIMD<int8_t>::broadcast(int8_t i) {
    return _mm_set1_epi8(i);
}
template<> inline __m128i SIMD<int8_t>::range() {
    static const __m128i r = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    return r;
}
template<> inline __m128i SIMD<int8_t>::eq(__m128i a, __m128i b) {
    return _mm_cmpeq_epi8(a, b);
}
template<> inline __m128i SIMD<int8_t>::less(__m128i a, __m128i b) {
    return _mm_cmplt_epi8(a, b);
}
template<> inline __m128i SIMD<int8_t>::max(__m128i a, __m128i b) {
    return _mm_max_epi8(a, b);
}
template<> inline __m128i SIMD<int8_t>::add(__m128i a, __m128i b) {
    return _mm_add_epi8(a, b);
}
template<> inline __m128i SIMD<int8_t>::subtract(__m128i a, __m128i b) {
    return _mm_sub_epi8(a, b);
}
//template<> inline __m128i SIMD<int8_t>::div2(__m128i a) {
//    // no 8-bit version, use 16 and then mask out the overflow bits
//    return _mm_and_si128(_mm_srli_epi16(a, 1), _mm_set1_epi8(0x3f));
//}

// 16 bit integers
template<> inline __m128i SIMD<int16_t>::broadcast(int16_t i) {
    return _mm_set1_epi16(i);
}
template<> inline __m128i SIMD<int16_t>::range() {
    static const __m128i r = _mm_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7);
    return r;
}
template<> inline __m128i SIMD<int16_t>::eq(__m128i a, __m128i b) {
    return _mm_cmpeq_epi16(a, b);
}
template<> inline __m128i SIMD<int16_t>::less(__m128i a, __m128i b) {
    return _mm_cmplt_epi16(a, b);
}
template<> inline __m128i SIMD<int16_t>::max(__m128i a, __m128i b) {
    return _mm_max_epi16(a, b);
}
template<> inline __m128i SIMD<int16_t>::add(__m128i a, __m128i b) {
    return _mm_add_epi16(a, b);
}
template<> inline __m128i SIMD<int16_t>::subtract(__m128i a, __m128i b) {
    return _mm_sub_epi16(a, b);
}
//template<> inline __m128i SIMD<int16_t>::div2(__m128i a) {
//    return _mm_srli_epi16(a, 1);
//}

// 32 bit integers
template<> inline __m128i SIMD<int32_t>::broadcast(int32_t i) {
    return _mm_set1_epi32(i);
}
template<> inline __m128i SIMD<int32_t>::range() {
    static const __m128i r = _mm_setr_epi32(0, 1, 2, 3);
    return r;
}
template<> inline __m128i SIMD<int32_t>::eq(__m128i a, __m128i b) {
    return _mm_cmpeq_epi32(a, b);
}
template<> inline __m128i SIMD<int32_t>::less(__m128i a, __m128i b) {
    return _mm_cmplt_epi32(a, b);
}
template<> inline __m128i SIMD<int32_t>::max(__m128i a, __m128i b) {
    return _mm_max_epi32(a, b);
}
template<> inline __m128i SIMD<int32_t>::add(__m128i a, __m128i b) {
    return _mm_add_epi32(a, b);
}
template<> inline __m128i SIMD<int32_t>::subtract(__m128i a, __m128i b) {
    return _mm_sub_epi32(a, b);
}
//template<> inline __m128i SIMD<int32_t>::div2(__m128i a) {
//    return _mm_srli_epi32(a, 1);
//}

template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align(const char* seq1, size_t len1,
                                  const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, StdMem, FwdCompareMatchFunc>(str1, str2,
                                                                  std::numeric_limits<uint64_t>::max());
}


template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_low_mem(const char* seq1, size_t len1,
                                          const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, LowMem, FwdCompareMatchFunc>(str1, str2,
                                                                  std::numeric_limits<uint64_t>::max());
}


template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_recursive(const char* seq1, size_t len1,
                                            const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, Recursive, FwdCompareMatchFunc>(str1, str2,
                                                                     std::numeric_limits<uint64_t>::max());
}

template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_adaptive(const char* seq1, size_t len1,
                                           const char* seq2, size_t len2,
                                           uint64_t max_mem) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, AdaptiveMem, FwdCompareMatchFunc>(str1, str2, max_mem);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner<NumPW>::wavefront_align_local(const char* seq1, size_t len1,
                                        const char* seq2, size_t len2,
                                        size_t anchor_begin_1, size_t anchor_end_1,
                                        size_t anchor_begin_2, size_t anchor_end_2,
                                        bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef StandardSemilocalWFA<RevCompareMatchFunc, RevStringView> PrefWFA;
    typedef StandardGlobalWFA<FwdCompareMatchFunc, StringView> AnchorWFA;
    typedef StandardSemilocalWFA<FwdCompareMatchFunc, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                   std::numeric_limits<uint64_t>::max(),
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner<NumPW>::wavefront_align_local_low_mem(const char* seq1, size_t len1,
                                                const char* seq2, size_t len2,
                                                size_t anchor_begin_1, size_t anchor_end_1,
                                                size_t anchor_begin_2, size_t anchor_end_2,
                                                bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef LowMemSemilocalWFA<RevCompareMatchFunc, RevStringView> PrefWFA;
    typedef LowMemGlobalWFA<FwdCompareMatchFunc, StringView> AnchorWFA;
    typedef LowMemSemilocalWFA<FwdCompareMatchFunc, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                   std::numeric_limits<uint64_t>::max(),
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner<NumPW>::wavefront_align_local_recursive(const char* seq1, size_t len1,
                                                  const char* seq2, size_t len2,
                                                  size_t anchor_begin_1, size_t anchor_end_1,
                                                  size_t anchor_begin_2, size_t anchor_end_2,
                                                  bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef RecursiveSemilocalWFA<RevCompareMatchFunc, RevStringView> PrefWFA;
    typedef RecursiveGlobalWFA<FwdCompareMatchFunc, StringView> AnchorWFA;
    typedef RecursiveSemilocalWFA<FwdCompareMatchFunc, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                   std::numeric_limits<uint64_t>::max(),
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner<NumPW>::wavefront_align_local_adaptive(const char* seq1, size_t len1,
                                                 const char* seq2, size_t len2,
                                                 uint64_t max_mem,
                                                 size_t anchor_begin_1, size_t anchor_end_1,
                                                 size_t anchor_begin_2, size_t anchor_end_2,
                                                 bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef AdaptiveSemilocalWFA<RevCompareMatchFunc, RevStringView> PrefWFA;
    typedef AdaptiveGlobalWFA<FwdCompareMatchFunc, StringView> AnchorWFA;
    typedef AdaptiveSemilocalWFA<FwdCompareMatchFunc, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2, max_mem,
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

template<int NumPW>
template<bool Local, int Mem, typename MatchFunc, typename StringType>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_dispatch(const StringType& seq1, const StringType& seq2, uint64_t max_mem) const {
    
    if (seq1.size() + seq2.size() > std::numeric_limits<int32_t>::max()) {
        std::stringstream strm;
        strm << "error:[WFAligner] WFA implementation can only align sequences of combined length <= "  << std::numeric_limits<int32_t>::max() << '\n';
        throw std::runtime_error(strm.str());
    }
    
    // make a match function as directed by the calling environment
    MatchFunc match_func(seq1, seq2);
    
    // decide how small of integers we can get away using, then do either the
    // low memory or standard algorithm with the desired configuration
    std::pair<std::vector<CIGAROp>, int32_t> result;
    // this is larger than the max DP value, we need it to avoid overflow in the boundary checking
    size_t max_anti_diag = 2 * std::max(seq1.size(), seq2.size());
    if (max_anti_diag <= (size_t) std::numeric_limits<int8_t>::max()) {
        if (Mem == AdaptiveMem) {
            result = wavefront_align_core<Local, true, int8_t>(seq1, seq2, match_func, max_mem);
        }
        else if (Mem == StdMem) {
            result = wavefront_align_core<Local, false, int8_t>(seq1, seq2, match_func,
                                                                std::numeric_limits<uint64_t>::max());
        }
        else if (Mem == LowMem) {
            result = wavefront_align_low_mem_core<Local, int8_t>(seq1, seq2, match_func);
        }
        else {
            result = wavefront_align_recursive_core<Local, int8_t>(seq1, seq2, match_func);
        }
    }
    else if (max_anti_diag <= (size_t) std::numeric_limits<int16_t>::max()) {
        if (Mem == AdaptiveMem) {
            result = wavefront_align_core<Local, true, int16_t>(seq1, seq2, match_func, max_mem);
        }
        else if (Mem == StdMem) {
            result = wavefront_align_core<Local, false, int16_t>(seq1, seq2, match_func,
                                                                 std::numeric_limits<uint64_t>::max());
        }
        else if (Mem == LowMem) {
            result = wavefront_align_low_mem_core<Local, int16_t>(seq1, seq2, match_func);
        }
        else {
            result = wavefront_align_recursive_core<Local, int16_t>(seq1, seq2, match_func);
        }
    }
    else {
        if (Mem == AdaptiveMem) {
            result = wavefront_align_core<Local, true, int32_t>(seq1, seq2, match_func, max_mem);
        }
        else if (Mem == StdMem) {
            result = wavefront_align_core<Local, false, int32_t>(seq1, seq2, match_func,
                                                                 std::numeric_limits<uint64_t>::max());
        }
        else if (Mem == LowMem) {
            result = wavefront_align_low_mem_core<Local, int32_t>(seq1, seq2, match_func);
        }
        else {
            result = wavefront_align_recursive_core<Local, int32_t>(seq1, seq2, match_func);
        }
    }
    // TODO: not doing int64_t for now -- we have bigger problems if we're aligning
    // gigabase sequences
    
    // unreduce the score
    result.second *= factor;
    
    // for local alignment, the conversion back will take place in the local core
    // because we need to measure the length of the aligned portion
    if (match && !Local) {
        result.second = convert_score(seq1.size(), seq2.size(), result.second);
    }
    
    return result;
}

template<int NumPW>
template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_core(const StringType& seq1, const StringType& seq2,
                                       const MatchFunc& match_func, uint64_t max_mem) const {

    // the return value
    std::pair<std::vector<CIGAROp>, int32_t> res;
    
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    
    uint64_t curr_mem = 0;
    bool hit_max_mem = false;
    
    // init wavefront to the upper left of both sequences' starts
    std::vector<Wavefront<IntType>> wfs;
    wfs.emplace_back(0, 1);
    wfs.front().int_at(0, 0) = 0;
    
    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wfs.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, 0, wfs.front(),
                       opt, opt_diag, opt_s, max_s);
    }
    if (Adaptive) {
        curr_mem += wfs.front().bytes();
    }
    
    while (!wavefront_reached(wfs.back(), final_diag, final_anti_diag) &&
           wfs.size() - 1 <= max_s) {
        if (Adaptive && curr_mem > max_mem) {
            hit_max_mem = true;
            break;
        }
        
        wfs.emplace_back(wavefront_next<IntType>(seq1, seq2, wfs));
        
        wavefront_extend(seq1, seq2, wfs.back(), match_func);
        
        if (Local) {
            find_local_opt(seq1, seq2, wfs.size() - 1, wfs.back(),
                           opt, opt_diag, opt_s, max_s);
        }
        // prune lagging diagonals
        wavefront_prune<Local>(wfs.back());
        
        if (Adaptive) {
            curr_mem += wfs.back().bytes();
        }
    }
    
#ifdef debug_viz
    wavefront_viz(seq1, seq2, wfs, scores);
#endif
    
    if (Adaptive && hit_max_mem) {
        // we hit the memory maximum without completing, fall back to the low memory algorithm
        res = fall_back_to_low_mem<Local, IntType>(seq1, seq2, match_func, max_mem,
                                                   opt, opt_diag, opt_s, max_s, wfs);
    }
    else {
        // we completed the algorithm
        int64_t traceback_s, traceback_d;
        if (Local) {
            traceback_d = opt_diag;
            traceback_s = opt_s;
        }
        else {
            traceback_d = seq1.size() - seq2.size();
            traceback_s = wfs.size() - 1;
        }
        res.second = traceback_s;
        res.first = wavefront_traceback(seq1, seq2, wfs, traceback_s, traceback_d);
    }
    
    return res;
}

template<int NumPW>
template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::fall_back_to_low_mem(const StringType& seq1, const StringType& seq2,
                                       const MatchFunc& match_func, uint64_t max_mem,
                                       int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                       std::vector<WFAligner::Wavefront<IntType>>& wfs) const {
    
    std::deque<Wavefront<IntType>> wf_buffer;
    std::vector<std::pair<int32_t, Wavefront<IntType>>> wf_bank;
    
    // subsampling params
    int64_t epoch_len = 1;
    int64_t epoch_end = 1;
    int64_t sample_rate = 1;
    
    uint64_t curr_mem_bank = 0;
    uint64_t curr_mem_block = 0;
    uint64_t max_mem_block = 0;
    uint64_t curr_mem_buffer = 0;
    
    size_t s = wfs.size() - 1;
    if (wfs.size() < 2 * stripe_width) {
        // we wouldn't have started evicting yet, everything should end up in the buffer
        for (auto& wf : wfs) {
            curr_mem_buffer += wf.bytes();
            wf_buffer.emplace_back(std::move(wf));
        }
    }
    else {
        // grab the stripe(s) that would have been found in the buffer
        size_t buffer_stripe_num = wfs.size() / stripe_width - 1;
        for (size_t i = buffer_stripe_num * stripe_width; i < wfs.size(); ++i) {
            wf_buffer.emplace_back(std::move(wfs[i]));
        }
        
        // sub-sample the stripes according the the epoch schedule
        // TODO: kind of a lot of copy-pasta code here...
        for (int64_t stripe_num = 0; stripe_num < buffer_stripe_num; ++stripe_num) {
            
            // is this the first wavefront of a new epoch?
            if (stripe_num == epoch_end) {
                // we're entering a new epoch, adust the sub-sampling parameters
                // accordingly
                epoch_len *= 4;
                epoch_end += epoch_len;
                sample_rate *= 2;
            }
            
            // bit-hacky sub for % that works because sample_rate is a power of 2
            bool keep = ((stripe_num & (sample_rate - 1)) == 0);
            for (size_t i = 0, k = stripe_num * stripe_width; i < stripe_width; ++i) {
                if (keep) {
                    // we're ejecting a stripe that's being retained,
                    wf_bank.emplace_back(k + i, std::move(wfs[k + i]));
                    curr_mem_bank += wf_bank.back().second.bytes();
                    curr_mem_block = 0;
                }
                else {
                    curr_mem_block += wfs[k + i].bytes();
                    max_mem_block = std::max(max_mem_block, curr_mem_block);
                }
            }
        }
    }
    
    // clear the memory that we're no longer using
    wfs = std::vector<Wavefront<IntType>>();
    
    return wavefront_align_low_mem_core_internal<Local, true, IntType>(seq1, seq2, match_func, max_mem,
                                                                       epoch_len, epoch_end, sample_rate,
                                                                       opt, opt_diag, opt_s, max_s,
                                                                       wf_buffer, s, wf_bank,
                                                                       curr_mem_buffer, curr_mem_bank, curr_mem_block, max_mem_block);
}

template<int NumPW>
template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_low_mem_core(const StringType& seq1, const StringType& seq2,
                                               const MatchFunc& match_func) const {

    // subsampling params
    int64_t epoch_len = 1;
    int64_t epoch_end = 1;
    int64_t sample_rate = 1;
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    
    // init wavefront to the upper left of both sequences' starts
    std::deque<Wavefront<IntType>> wf_buffer;
    wf_buffer.emplace_back(0, 1);
    wf_buffer.front().int_at(0, 0) = 0;
    size_t s = 0;
    
    // make sure the first wavefront is extended
    wavefront_extend(seq1, seq2, wf_buffer.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, s, wf_buffer.front(),
                       opt, opt_diag, opt_s, max_s);
    }
    
    // the wavefronts we've decided to keep after they leave the buffer
    std::vector<std::pair<int32_t, Wavefront<IntType>>> wf_bank;
    
    uint64_t curr_mem = wf_buffer.front().bytes();
    
    // TODO: i'm not sure if all of the tracking variables will be able to get optimized out
    // since it will be hard to identify them as constant. compiler might not respect inlining
    // of internal method
    return wavefront_align_low_mem_core_internal<Local, false, IntType>(seq1, seq2, match_func,
                                                                        std::numeric_limits<uint64_t>::max(),
                                                                        epoch_len, epoch_end, sample_rate,
                                                                        opt, opt_diag, opt_s, max_s,
                                                                        wf_buffer, s, wf_bank, curr_mem, 0, 0, 0);
};

template<int NumPW>
template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_low_mem_core_internal(const StringType& seq1, const StringType& seq2,
                                                        const MatchFunc& match_func, uint64_t max_mem,
                                                        int64_t epoch_len, int64_t epoch_end, int64_t sample_rate,
                                                        int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                                        std::deque<WFAligner<NumPW>::Wavefront<IntType>>& wf_buffer, size_t s,
                                                        std::vector<std::pair<int32_t, WFAligner<NumPW>::Wavefront<IntType>>>& wf_bank,
                                                        uint64_t curr_mem_buffer, uint64_t curr_mem_bank,
                                                        uint64_t curr_mem_block, uint64_t max_mem_block) const {

    // the return value
    std::pair<std::vector<CIGAROp>, int32_t> res;
    
    // use these to know when we've finished the alignment
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();
    
    bool hit_max_mem = false;
    
    // do wavefront iterations until hit max score or alignment finishes
    while (!wavefront_reached(wf_buffer.back(), final_diag, final_anti_diag) &&
           s <= max_s) {
        
        // increase score and get the next wavefront
        ++s;
        wf_buffer.emplace_back(wavefront_next<IntType>(seq1, seq2, wf_buffer));
        if (Adaptive) {
            curr_mem_buffer += wf_buffer.back().bytes();
        }
        
        // note: the method of keeping up to 2 stripes in memory at a time is a bit strange
        // here, but it greatly simplifies the fall back to the recursive alignment
        if (wf_buffer.size() == 2 * stripe_width) {
            // we need to evict a stripe from the buffer to continue
            
            // the stripe of the wavefront that is falling out of the buffer
            int32_t stripe_num = s / stripe_width - 1; // formula works because s must be k*w-1
            
            // is this the first wavefront of a new epoch?
            if (stripe_num == epoch_end) {
                // we're entering a new epoch, adust the sub-sampling parameters
                // accordingly
                epoch_len *= 4;
                epoch_end += epoch_len;
                sample_rate *= 2;
            }
            
            // bit-hacky sub for % that works because sample_rate is a power of 2
            bool keep = ((stripe_num & (sample_rate - 1)) == 0);
            for (size_t i = 0; i < stripe_width; ++i) {
                if (keep) {
                    // we're ejecting a stripe that's being retained,
                    wf_bank.emplace_back(stripe_num * stripe_width + i, std::move(wf_buffer.front()));
                    if (Adaptive) {
                        curr_mem_block = 0;
                        curr_mem_bank += wf_buffer.front().bytes();
                    }
                }
                else if (Adaptive) {
                    // record that we're freeing memory
                    auto b = wf_buffer.front().bytes();
                    curr_mem_buffer -= b;
                    curr_mem_block += b;
                    max_mem_block = std::max(curr_mem_block, max_mem_block);
                }
                wf_buffer.pop_front();
            }
        }
        
        // follow matches
        wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
        
        if (Local) {
            find_local_opt(seq1, seq2, s, wf_buffer.back(),
                           opt, opt_diag, opt_s, max_s);
        }
        // prune lagging diagonals
        wavefront_prune<Local>(wf_buffer.back());
        
        if (Adaptive) {
            curr_mem_bank += wf_buffer.back().bytes();
            if (curr_mem_bank + max_mem_block + curr_mem_buffer >= max_mem) {
                hit_max_mem = true;
                break;
            }
        }
    }
    
    if (Adaptive && hit_max_mem) {
        // we hit the maximum memory without completing, fall back to the recursive algorithm
        res = fall_back_to_recursive<Local, IntType>(seq1, seq2, match_func,
                                                     opt, opt_diag, opt_s, max_s,
                                                     wf_buffer, s, wf_bank);
    }
    else {
        // we completed successfully
        
        // move the final wavefronts from the buffer to the bank
        for (size_t i = 0; i < wf_buffer.size(); ++i) {
            wf_bank.emplace_back(s - wf_buffer.size() + i + 1, std::move(wf_buffer[i]));
        }
        
        int64_t traceback_s, traceback_d;
        if (Local) {
            traceback_d = opt_diag;
            traceback_s = opt_s;
            s = opt_s;
        }
        else {
            traceback_d = seq1.size() - seq2.size();
            traceback_s = s;
        }
        
        res.first = wavefront_traceback_low_mem<Local>(seq1, seq2, wf_bank, traceback_s, traceback_d, match_func);
        res.second = s;
    }
    
    return res;
}

template<int NumPW>
template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::fall_back_to_recursive(const StringType& seq1, const StringType& seq2,
                                         const MatchFunc& match_func,
                                         int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                         std::deque<WFAligner<NumPW>::Wavefront<IntType>>& wf_buffer, size_t s,
                                         std::vector<std::pair<int32_t, WFAligner<NumPW>::Wavefront<IntType>>>& wf_bank) const {
    
    // convert the low memory data structures into the recursive data structures
    
    std::vector<Wavefront<IntType>> opt_stripe;
    std::vector<Wavefront<IntType>> first_stripe;
    
    if (!wf_bank.empty()) {
        // we ejected a stripe from the buffer, which we guarantee must include
        // the first stripe, so now we remember it
        first_stripe.reserve(stripe_width);
        for (size_t i = 0; i < stripe_width; ++i) {
            first_stripe.emplace_back(std::move(wf_bank[i].second));
        }
        
        if (Local) {
            
            int64_t opt_stripe_num = opt_s / stripe_width;
            if (opt_stripe_num != 0 && opt_stripe_num < (s - wf_buffer.size() + 1) / stripe_width) {
                // the opt isn't in the first stripe or the buffer, so we need to either find or recompute it
                
                size_t i = 0;
                while (wf_bank[i].first / stripe_width < opt_stripe_num) {
                    i += stripe_width;
                }
                
                if (wf_bank[i].first / stripe_width != opt_stripe_num) {
                    // we didn't sub-sample the stripe with the opt, so we have recompute it
                    
                    // eliminate the portion after the opt, which we don't need
                    wf_bank.resize(i);
                    // and add another stripe to the bank
                    WFBankAdapter<IntType> adapter(wf_bank);
                    for (size_t j = 0; j < stripe_width; ++j) {
                        wf_bank.emplace_back(wf_bank.back().first + 1,
                                             wavefront_next<IntType>(seq1, seq2, adapter));
                        wavefront_extend(seq1, seq2, wf_bank.back().second, match_func);
                        wavefront_prune<Local>(wf_bank.back().second);
                    }
                }
                
                // retrieve the opt's stripe from the bank
                opt_stripe.reserve(stripe_width);
                for (size_t n = i + stripe_width; i < n; ++i) {
                    opt_stripe.emplace_back(std::move(wf_bank[i].second));
                }
            }
        }
    }
    
    // discard the wavefronts we're no longer retaining
    wf_bank = std::vector<std::pair<int32_t, Wavefront<IntType>>>();
    
    return wavefront_align_recursive_core_internal<Local, IntType>(seq1, seq2, match_func, opt, opt_diag, opt_s, max_s,
                                                                   opt_stripe, first_stripe, wf_buffer, s);
}

template<int NumPW>
template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_recursive_core(const StringType& seq1, const StringType& seq2,
                                                 const MatchFunc& match_func) const {
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    std::vector<Wavefront<IntType>> opt_stripe;
    
    // init wavefront to the upper left of both sequences' starts
    std::deque<Wavefront<IntType>> wf_buffer;
    wf_buffer.emplace_back(0, 1);
    wf_buffer.front().int_at(0, 0) = 0;
    size_t s = 0;
    
    wavefront_extend(seq1, seq2, wf_buffer.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, s, wf_buffer.front(),
                       opt, opt_diag, opt_s, max_s);
    }
    
    // we need to remember the first stripe for the recursive step
    std::vector<Wavefront<IntType>> first_stripe;
    
    return wavefront_align_recursive_core_internal<Local, IntType>(seq1, seq2, match_func, opt, opt_diag, opt_s, max_s,
                                                                   opt_stripe, first_stripe, wf_buffer, s);
}

template<int NumPW>
template<bool Local, typename IntType, typename StringType, typename MatchFunc>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner<NumPW>::wavefront_align_recursive_core_internal(const StringType& seq1, const StringType& seq2,
                                                          const MatchFunc& match_func,
                                                          int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                                          std::vector<WFAligner<NumPW>::Wavefront<IntType>>& opt_stripe,
                                                          std::vector<WFAligner<NumPW>::Wavefront<IntType>>& first_stripe,
                                                          std::deque<WFAligner<NumPW>::Wavefront<IntType>>& wf_buffer, size_t s) const {
    
    // use these to know when we've finished the alignment
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();
    
    // do wavefront iterations until hit max score or alignment finishes
    while (!wavefront_reached(wf_buffer.back(), final_diag, final_anti_diag) &&
           s <= max_s) {
        // increase score and get the next wavefront
        ++s;
        wf_buffer.emplace_back(wavefront_next<IntType>(seq1, seq2, wf_buffer));
        
        if (wf_buffer.size() == 2 * stripe_width) {
            // we need to evict a stripe wavefront from the buffer to continue
            for (int64_t i = 0; i < stripe_width; ++i) {
                if (s + 1 == 2 * stripe_width) {
                    // we need to remember the first stripe when we evict it
                    first_stripe.emplace_back(std::move(wf_buffer.front()));
                }
                else if (Local && opt_s > s - 2 * stripe_width && opt_s <= s - stripe_width) {
                    // we're ejecting the stripe that contains the local opt, so
                    // we remember it for later
                    if (i == 0) {
                        // we might have an existing stripe in the opt stripe
                        opt_stripe.clear();
                    }
                    opt_stripe.emplace_back(std::move(wf_buffer.front()));
                }
                wf_buffer.pop_front();
            }
            
        }
        
        // follow matches
        wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
        
        if (Local) {
            find_local_opt(seq1, seq2, s, wf_buffer.back(),
                           opt, opt_diag, opt_s, max_s);
        }
        // prune lagging diagonals
        wavefront_prune<Local>(wf_buffer.back());
    }
    
    int64_t traceback_d;
    int64_t traceback_s;
    int32_t score;
    if (Local) {
        traceback_d = opt_diag;
        traceback_s = opt_s;
        score = opt_s;
        
        int64_t opt_stripe_num = opt_s / stripe_width;
        if (opt_stripe_num != 0 && opt_stripe_num < (s - wf_buffer.size() + 1) / stripe_width) {
            // the local optimum did not occur in the first or last stripe, so we replace the contents
            // of the last stripe with the opt's stripe, which we should have retained
            // we maintain the final position in the buffer is where traceback should begin
            wf_buffer.clear();
            for (size_t i = 0, n = opt_s - opt_stripe_num * stripe_width; i <= n; ++i) {
                wf_buffer.emplace_back(std::move(opt_stripe[i]));
            }
        }
    }
    else {
        traceback_d = final_diag;
        traceback_s = s;
        score = s;
    }
    
    std::vector<CIGAROp> cigar;
    if (s < 2 * stripe_width) {
        // we finished the entire alignment without evicting, so we just do normal traceback
        // TODO: i should be able to handle this with a template, but it's having trouble with
        // type inference...
        std::vector<Wavefront<IntType>> traceback_wfs;
        traceback_wfs.reserve(wf_buffer.size());
        for (auto& wf : wf_buffer) {
            traceback_wfs.emplace_back(std::move(wf));
        }
        cigar = wavefront_traceback(seq1, seq2, traceback_wfs,
                                    traceback_s, traceback_d);
    }
    else if (Local && traceback_s < stripe_width) {
        // we encountered the local optimum within the first stripe, so we can skip complicated
        // traceback logic
        cigar = wavefront_traceback(seq1, seq2, first_stripe,
                                    traceback_s, traceback_d);
    }
    else {
        // we enter the recursive portion of the algorithm
        int64_t stripe_num = (traceback_s - wf_buffer.size() + 1) / stripe_width;
        int traceback_mat = 0;
        int64_t lead_matches = 0;
        wavefront_align_recursive_internal<Local, IntType>(seq1, seq2,
                                                           0, first_stripe,
                                                           stripe_num, wf_buffer,
                                                           traceback_s, traceback_d,
                                                           traceback_mat, lead_matches,
                                                           cigar, match_func);
        assert(traceback_s == 0 && traceback_d == 0 && traceback_mat == 0);
        // the CIGAR is constructed in reverse, unreverse it
        std::reverse(cigar.begin(), cigar.end());
        // also merge operations from adjacent blocks
        coalesce_cigar(cigar);
    }
    
    
    return std::make_pair(cigar, int32_t(score));
}


template<int NumPW>
template<bool Local, typename IntType, typename StringType, typename MatchFunc, typename WFVector1, typename WFVector2>
void WFAligner<NumPW>::wavefront_align_recursive_internal(const StringType& seq1, const StringType& seq2,
                                                          int64_t lower_stripe_num, WFVector1& lower_stripe,
                                                          int64_t upper_stripe_num, WFVector2& upper_stripe,
                                                          int64_t& traceback_s, int64_t& traceback_d,
                                                          int& traceback_mat, int64_t& lead_matches,
                                                          std::vector<CIGAROp>& cigar,
                                                          const MatchFunc& match_func) const {
    
    if (lower_stripe_num + 1 == upper_stripe_num) {
        // do a step of traceback
        
        // collect all of the wavefronts in one vector
        for (auto& wf : upper_stripe) {
            lower_stripe.emplace_back(std::move(wf));
        }
        
        
        // convert the current s into an index within this block
        int64_t relative_s = traceback_s - lower_stripe_num * stripe_width;
        
        // traceback the stripe
        wavefront_traceback_internal(seq1, seq2, lower_stripe,
                                     traceback_d, lead_matches, traceback_mat,
                                     relative_s, cigar);
        
        // convert the final within-block s back into the real s
        traceback_s = lower_stripe_num * stripe_width + relative_s;
        
        // we no longer need the upper stripe for anything, we can discard it
        while (lower_stripe.size() > stripe_width) {
            lower_stripe.pop_back();
        }
    }
    else {
        // borrow the parent call's wavefronts to do our own DP with
        std::deque<Wavefront<IntType>> wf_buffer;
        //wf_buffer.reserve(2 * stripe_width);
        for (auto& wf : lower_stripe) {
            wf_buffer.emplace_back(std::move(wf));
        }
        
        // do DP out to the middle of the row interval
        int64_t target_stripe_num = (lower_stripe_num + upper_stripe_num) / 2;
        for (int64_t i = 0, n = stripe_width * (target_stripe_num - lower_stripe_num); i < n; ++i) {
            // one step of DP
            wf_buffer.emplace_back(wavefront_next<IntType>(seq1, seq2, wf_buffer));
            wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
            // prune lagging diagonals
            wavefront_prune<Local>(wf_buffer.back());
            
            if (wf_buffer.size() == 2 * stripe_width) {
                // we've completed a stripe, we can now evict the previous one
                for (int64_t j = 0; j < stripe_width; ++j) {
                    if (i < stripe_width) {
                        // we're evicting the first stripe, so we have to give it back
                        // to the parent call
                        lower_stripe[j] = std::move(wf_buffer.front());
                        
                    }
                    wf_buffer.pop_front();
                }
            }
        }
        
        // recurse into upper half
        wavefront_align_recursive_internal<Local, IntType>(seq1, seq2,
                                                           target_stripe_num, wf_buffer,
                                                           upper_stripe_num, upper_stripe,
                                                           traceback_s, traceback_d,
                                                           traceback_mat, lead_matches,
                                                           cigar, match_func);
        // recurse into lower half
        wavefront_align_recursive_internal<Local, IntType>(seq1, seq2,
                                                           lower_stripe_num, lower_stripe,
                                                           target_stripe_num, wf_buffer,
                                                           traceback_s, traceback_d,
                                                           traceback_mat, lead_matches,
                                                           cigar, match_func);
    }
}

template<int NumPW>
template<typename PrefSemilocalWFA, typename AnchorGlobalWFA, typename SuffSemilocalWFA>
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner<NumPW>::wavefront_align_local_core(const char* seq1, size_t len1,
                                             const char* seq2, size_t len2,
                                             uint64_t max_mem,
                                             size_t anchor_begin_1, size_t anchor_end_1,
                                             size_t anchor_begin_2, size_t anchor_end_2,
                                             bool anchor_is_match) const {
    
    if (match == 0) {
        std::stringstream strm;
        strm << "error:[WFAligner] aligner must be initialized with Smith-Waterman-Gotoh-style scoring parameters to perform local alignment" << std::endl;
        throw std::runtime_error(strm.str());
    }
    
    // tail align the prefix
    RevStringView pref1(seq1, 0, anchor_begin_1);
    RevStringView pref2(seq2, 0, anchor_begin_2);
    
    auto pref_result = PrefSemilocalWFA(this, max_mem)(pref1, pref2);
    
    std::vector<CIGAROp> cigar = std::move(pref_result.first);
    std::reverse(cigar.begin(), cigar.end());
    auto aligned_pref_len = cigar_base_length(cigar);
    int32_t pref_score = convert_score(aligned_pref_len.first, aligned_pref_len.second,
                                       pref_result.second);
    
    int32_t anchor_score = 0;
    if (!anchor_is_match) {
        
        // align the anchor
        StringView anchor1(seq1, anchor_begin_1, anchor_end_1 - anchor_begin_1);
        StringView anchor2(seq2, anchor_begin_2, anchor_end_2 - anchor_begin_2);
        
        auto anchor_res = AnchorGlobalWFA(this, max_mem)(anchor1, anchor2);
        
        // get the score and add to the CIGAR
        anchor_score = convert_score(anchor_end_1 - anchor_begin_1, anchor_end_2 - anchor_begin_2,
                                     anchor_res.second);
        for (auto& op : anchor_res.first) {
            cigar.emplace_back(std::move(op));
        }
    }
    else {
        // infer the anchor match
        auto len = anchor_end_1 - anchor_begin_1;
        if (len != anchor_end_2 - anchor_begin_2) {
            throw std::runtime_error("error: match anchor for local alignment has mismatched lengths");
        }
        anchor_score = match * len;
        cigar.emplace_back('M', len);
    }
    
    // tail align the suffix
    StringView suff1(seq1, anchor_end_1, len1 - anchor_end_1);
    StringView suff2(seq2, anchor_end_2, len2 - anchor_end_2);
    
    auto suff_result = SuffSemilocalWFA(this, max_mem)(suff1, suff2);
    
    auto aligned_suff_len = cigar_base_length(suff_result.first);
    int32_t suff_score = convert_score(aligned_suff_len.first, aligned_suff_len.second,
                                       suff_result.second);
    
    for (auto& op : suff_result.first) {
        cigar.emplace_back(std::move(op));
    }
    
    // CIGAR operations might not be merged over the boundaries
    coalesce_cigar(cigar);
    return std::make_tuple(std::move(cigar),
                           pref_score + anchor_score + suff_score,
                           std::make_pair(anchor_begin_1 - aligned_pref_len.first,
                                          anchor_end_1 + aligned_suff_len.first),
                           std::make_pair(anchor_begin_2 - aligned_pref_len.second,
                                          anchor_end_2 + aligned_suff_len.second));
}

template<int NumPW>
template<typename StringType, typename MatchFunc, typename IntType>
inline void WFAligner<NumPW>::wavefront_extend(const StringType& seq1, const StringType& seq2,
                                               Wavefront<IntType>& wf, const MatchFunc& match_func) const {
    // TODO: vectorize i and j computation
    // this currently doesn't work because i have signed values, so the bitshift doesn't accomplish
    // a divide. i'm also not sure it will actually be faster because of the serial final step that
    // can overshoot in narrow bands
//    for (int k = 0, n = SIMD<IntType>::round_up(wf.size()), diag = wf.diag_begin;
//         k < n; ++k, diag += SIMD<IntType>::vec_size()) {
//
//        // diagonal values
//        __m128i dv = SIMD<IntType>::add(SIMD<IntType>::broadcast(diag), SIMD<IntType>::range());
//        // seq1 indexes
//        __m128i iv = SIMD<IntType>::div2(SIMD<IntType>::add(wf.M[k], dv));
//        // seq2 indexes
//        __m128i jv = SIMD<IntType>::div2(SIMD<IntType>::subtract(wf.M[k], dv));
//        // TODO: this can overshoot and slow it down in a narrow band...
//        for (int l = 0; l < SIMD<IntType>::vec_size(); ++l) {
//            int64_t i = *(((IntType*) &iv) + l);
//            int64_t j = *(((IntType*) &jv) + l);
//            if (i >= 0 && j >= 0 && i < seq1.size() && j < seq2.size()) {
//                wf.int_at(*(((IntType*) &dv) + l), 0) += 2 * match_func(i, j);
//            }
//        }
//    }
    for (int64_t k = 0; k < wf.size(); ++k) {
        int64_t diag = wf.diag_begin + k;
        int64_t anti_diag = wf.int_at(diag, 0);
        int64_t i = (diag + anti_diag) / 2;
        int64_t j = (anti_diag - diag) / 2;
        if (i >= 0 && j >= 0) {
            wf.int_at(diag, 0) += 2 * match_func(i, j);
        }
    }
}

template<int NumPW>
template<typename IntType, typename WFVector, typename StringType>
inline typename WFAligner<NumPW>::template Wavefront<IntType>
WFAligner<NumPW>::wavefront_next(const StringType& seq1, const StringType& seq2,
                                 const WFVector& wfs) const {
    
    // find the range of incoming wavefronts
    int32_t hi = std::numeric_limits<int32_t>::min();
    int32_t lo = std::numeric_limits<int32_t>::max();
    // mismatch stays in the same diagonals
    if (wfs.size() >= mismatch) {
        const auto& wf_prev = wfs[wfs.size() - mismatch];
        if (!wf_prev.empty()) {
            lo = std::min<int32_t>(lo, wf_prev.diag_begin);
            hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.size());
        }
    }
    // TODO: but sometimes these bounds are slightly too large because of the gap extend
    // but it should be only ~s * o extra work
    // insertion/deletion expand the diagonals by 1
    if (NumPW == 0) {
        // linear scoring
        if (wfs.size() >= gap_extend[0]) {
            const auto& wf_prev = wfs[wfs.size() - gap_extend[0]];
            if (!wf_prev.empty()) {
                lo = std::min<int32_t>(lo, wf_prev.diag_begin - 1);
                hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.size() + 1);
            }
        }
    }
    else {
        // piecewise affine scoring
        for (int i = 0; i < NumPW; ++i) {
            auto s_back = gap_extend[i] + gap_open[i];
            if (wfs.size() >= s_back) {
                const auto& wf_prev = wfs[wfs.size() - s_back];
                if (!wf_prev.empty()) {
                    lo = std::min<int32_t>(lo, wf_prev.diag_begin - 1);
                    hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.size() + 1);
                }
            }
            s_back = gap_extend[i];
            if (wfs.size() >= s_back) {
                const auto& wf_prev = wfs[wfs.size() - s_back];
                if (!wf_prev.empty()) {
                    lo = std::min<int32_t>(lo, wf_prev.diag_begin - 1);
                    hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.size() + 1);
                }
            }
        }
    }
    
    // some twiddly code so that we can call constructor only once and preserve RVO
    int32_t ctor_lo, ctor_hi;
    if (hi < lo) {
        ctor_lo = 0;
        ctor_hi = 0;
    }
    else {
        ctor_lo = lo;
        ctor_hi = hi;
    }
    Wavefront<IntType> wf(ctor_lo, ctor_hi);
    
    if (lo < hi) {
        
        __m128i seq1_lim = SIMD<IntType>::broadcast(2 * seq1.size());
        __m128i seq2_lim = SIMD<IntType>::broadcast(2 * seq2.size());
        
        // a general macro for a vectorized DP step in wf_next
        #define __wfalm__wavefront_next_DP(MAT_FROM_ID, MAT_INTO_NAME, SHIFT, INCR, MERGE) \
        for (int begin = wf_prev.diag_begin - wf.diag_begin + (SHIFT), \
             k = SIMD<IntType>::round_down(begin), \
             n = SIMD<IntType>::round_up(begin + wf_prev.size()), \
             first_diag = wf.diag_begin + k * SIMD<IntType>::vec_size(); \
             k < n; ++k, first_diag += SIMD<IntType>::vec_size()) { \
            /* the DP values (we assume that the previous vector won't be aligned to the same index) */ \
            __m128i a = SIMD<IntType>::add(SIMD<IntType>::load_unaligned(wf_prev.vector_at(first_diag - (SHIFT), (MAT_FROM_ID))), \
                                           SIMD<IntType>::broadcast(INCR)); \
            /* the diagonal values */ \
            __m128i d = SIMD<IntType>::add(SIMD<IntType>::broadcast(first_diag), SIMD<IntType>::range()); \
            /* check that it doesn't go outside the matrix */ \
            __m128i bnd = SIMD<IntType>::mask_and(SIMD<IntType>::leq(SIMD<IntType>::add(a, d), seq1_lim), \
                                                  SIMD<IntType>::leq(SIMD<IntType>::subtract(a, d), seq2_lim)); \
            /* store all of the values that were inside the matrix */ \
            if (MERGE) { \
                (MAT_INTO_NAME)[k] = SIMD<IntType>::max(SIMD<IntType>::ifelse(bnd, a, SIMD<IntType>::min_inf()), (MAT_INTO_NAME)[k]); \
            } \
            else { \
                (MAT_INTO_NAME)[k] = SIMD<IntType>::ifelse(bnd, a, SIMD<IntType>::min_inf()); \
            } \
        }
        
        // mismatches
        if (wfs.size() >= mismatch) {
            const auto& wf_prev = wfs[wfs.size() - mismatch];
            
            // from 0 (i.e. M)
            // into wf.M
            // shift diag by 0
            // increment by 2
            // do not merge with existing values
            __wfalm__wavefront_next_DP(0, wf.M, 0, 2, false);
        }
        
        if (NumPW == 0) {
            // single matrix linear gaps
            
            // TODO: pretty repetitive with piecewise affine
            // take gaps
            if (wfs.size() >= gap_extend[0]) {
                
                const auto& wf_prev = wfs[wfs.size() - gap_extend[0]];
                
                // from 0 (i.e. M)
                // into wf.M
                // shift diag by 1
                // increment by 1
                // merge with existing values
                __wfalm__wavefront_next_DP(0, wf.M, 1, 1, true);
                
                // from 0 (i.e. M)
                // into wf.M
                // shift diag by -1
                // increment by 1
                // merge with existing values
                __wfalm__wavefront_next_DP(0, wf.M, -1, 1, true);
            }
        }
        else {
            // piecewise affine gaps in separate matrices
            
            for (int i = 0; i < NumPW; ++i) {
                
                int mat_num = i + 1;
                auto I = wf.M - mat_num * wf.interval;
                auto D = wf.M + mat_num * wf.interval;
                
                // open gaps
                if (wfs.size() >= gap_extend[i] + gap_open[i]) {
                    
                    const auto& wf_prev = wfs[wfs.size() - gap_extend[i] - gap_open[i]];
                    
                    // from 0 (i.e. M)
                    // into I
                    // shift diag by 1
                    // increment by 1
                    // do not merge with existing values
                    __wfalm__wavefront_next_DP(0, I, 1, 1, false);
                    
                    // from 0 (i.e. M)
                    // into D
                    // shift diag by -1
                    // increment by 1
                    // do not merge with existing values
                    __wfalm__wavefront_next_DP(0, D, -1, 1, false);
                }
                
                // extend gaps
                if (wfs.size() >= gap_extend[i]) {
                    
                    const auto& wf_prev = wfs[wfs.size() - gap_extend[i]];
                    
                    // from -mat_num (i.e. I)
                    // into I
                    // shift diag by 1
                    // increment by 1
                    // merge with existing values
                    __wfalm__wavefront_next_DP(-mat_num, I, 1, 1, true);
                    
                    // from mat_num (i.e. D)
                    // into D
                    // shift diag by -1
                    // increment by 1
                    // merge with existing values
                    __wfalm__wavefront_next_DP(mat_num, D, -1, 1, true);
                }
            }
            
            // calculate the opt
            for (int k = 0, n = SIMD<IntType>::round_up(wf.size()); k < n; ++k) {
                for (int i = 1; i <= NumPW; i++) {
                    auto I = wf.M - i * wf.interval;
                    auto D = wf.M + i * wf.interval;
                    wf.M[k] = SIMD<IntType>::max(wf.M[k], SIMD<IntType>::max(I[k], D[k]));
                }
            }
        }
        #undef __wfalm__wavefront_next_DP
    }
    
    return wf;
}

template<int NumPW>
template<bool Local, typename IntType>
inline void WFAligner<NumPW>::wavefront_prune(WFAligner<NumPW>::Wavefront<IntType>& wf) const {
    // TODO: for now, just skipping pruning for local alignment, but actually should come up
    // with a better pruning mechanism for it...
    if (!wf.empty() && lagging_diagonal_prune >= 0 && !Local) {

        // TODO: do this with horizontal max
        IntType max_anti_diag = std::numeric_limits<IntType>::min();
        for (int32_t i = 0; i < wf.size(); ++i) {
            max_anti_diag = std::max<IntType>(max_anti_diag, wf.int_at(i, 0));
        }

        // squeeze the diagonals
        // TODO: figure out a vectorized way to do this
        int64_t d_min = wf.diag_begin;
        int64_t d_max = wf.diag_begin + wf.size() - 1;
        while (wf.int_at(d_min, 0) < max_anti_diag - lagging_diagonal_prune) {
            ++d_min;
        }
        while (wf.int_at(d_max, 0) < max_anti_diag - lagging_diagonal_prune) {
            --d_max;
        }
        
        if (d_min != wf.diag_begin || d_max + 1 != wf.diag_begin + wf.size()) {
            wf.resize(d_min, d_max - d_min + 1);
        }
    }
}

template<int NumPW>
template<typename IntType>
inline bool WFAligner<NumPW>::wavefront_reached(const WFAligner<NumPW>::Wavefront<IntType>& wf,
                                                int32_t diag, int32_t anti_diag) const {
    if (diag >= wf.diag_begin && diag < wf.diag_begin + wf.size()) {
        return (wf.int_at(diag, 0) == anti_diag);
    }
    return false;
}

template<int NumPW>
template<typename StringType, typename IntType>
inline std::vector<CIGAROp> WFAligner<NumPW>::wavefront_traceback(const StringType& seq1, const StringType& seq2,
                                                                  const std::vector<WFAligner<NumPW>::Wavefront<IntType>>& wfs,
                                                                  int64_t s, int64_t d) const {
    
    // always begin traceback in the match matrix
    int mat = 0;
    int64_t lead_matches = 0;
    
    std::vector<CIGAROp> cigar;
    wavefront_traceback_internal(seq1, seq2, wfs, d, lead_matches, mat, s, cigar);
    assert(s == 0 && d == 0 && mat == 0);
    std::reverse(cigar.begin(), cigar.end());
    
    return cigar;
}

template<int NumPW>
template <bool Local, typename StringType, typename MatchFunc, typename IntType>
std::vector<CIGAROp> WFAligner<NumPW>::wavefront_traceback_low_mem(const StringType& seq1, const StringType& seq2,
                                                                   std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                                                   int64_t s, int64_t d, const MatchFunc& match_func) const {
    
    // the return value
    std::vector<CIGAROp> cigar;
    
    if (Local) {
        // if we overshot the opt traceback location, clear the way for it to be
        // in the last position of the bank
        while (wf_bank.back().first > s) {
            wf_bank.pop_back();
        }
        // and recompute out to the opt if necessary
        WFBankAdapter<IntType> adapter(wf_bank);
        while (wf_bank.back().first < s) {
            wf_bank.emplace_back(wf_bank.back().first + 1,
                                 wavefront_next<IntType>(seq1, seq2, adapter));
            wavefront_extend(seq1, seq2, wf_bank.back().second, match_func);
            wavefront_prune<Local>(wf_bank.back().second);
        }
        // the traceback starting location should now be in the final wavefront
        // in the bank, regardless of global/local
    }
    
    // a moving window over the DP that we will incrementally do traceback through
    std::deque<Wavefront<IntType>> traceback_block;
    // a pointer to the first position in wf_bank that has been moved into
    // the traceback block
    int64_t i = wf_bank.size() - 1;
    traceback_block.emplace_front(std::move(wf_bank[i].second));
    
    while (i - 1 >= 0 && wf_bank[i - 1].first == wf_bank[i].first - 1) {
        // the scores are contiguous
        --i;
        traceback_block.emplace_front(std::move(wf_bank[i].second));
    }
    // begin traceback in the match matrix
    int mat = 0;
    int64_t lead_matches = 0;
    // the index within the traceback block that traceback ends, initial value is arbitrary
    int64_t relative_s = traceback_block.size() - 1;

    // traceback as far as possible in this block
    wavefront_traceback_internal(seq1, seq2, traceback_block,
                                 d, lead_matches, mat, relative_s, cigar);
    // remove anything we don't need anymore
    traceback_block.resize(relative_s + 1);
    
    while (i > 0) {
        // there's an earlier discontiguous stripe we need to connect to
        
        // figure out where it starts
        int64_t j = i - 1;
        while (j - 1 >= 0 && wf_bank[j - 1].first == wf_bank[j].first - 1) {
            --j;
        }
        
        // initialize the full block with stripe
        std::vector<Wavefront<IntType>> fillin_wfs;
        fillin_wfs.reserve(wf_bank[i].first - wf_bank[j].first);
        for (int64_t k = j; k < i; ++k) {
            fillin_wfs.emplace_back(std::move(wf_bank[k].second));
        }
        
        // redo the DP between the stripes
        for (int32_t t = wf_bank[i - 1].first + 1, end = wf_bank[i].first; t < end; ++t) {
            fillin_wfs.emplace_back(wavefront_next<IntType>(seq1, seq2, fillin_wfs));
            wavefront_extend(seq1, seq2, fillin_wfs.back(), match_func);
            wavefront_prune<Local>(fillin_wfs.back());
        }
        
        // and move the redone wavefronts into the block buffer for traceback
        for (int64_t k = fillin_wfs.size() - 1; k >= 0; --k) {
            traceback_block.emplace_front(std::move(fillin_wfs[k]));
        }
        i = j;
        
        relative_s = traceback_block.size() - 1;
        wavefront_traceback_internal(seq1, seq2, traceback_block,
                                     d, lead_matches, mat, relative_s, cigar);
        // remove anything we don't need anymore
        traceback_block.resize(relative_s + 1);
    }
    
    // we sometimes need to coalesce CIGAR operations that occurred across
    // the boundaries of blocks
    coalesce_cigar(cigar);
    // the CIGAR is built backwards
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

template<int NumPW>
template<typename WFVector, typename StringType>
void WFAligner<NumPW>::wavefront_traceback_internal(const StringType& seq1, const StringType& seq2,
                                                    const WFVector& wfs, int64_t& d, int64_t& lead_matches,
                                                    int& mat, int64_t& s, std::vector<CIGAROp>& cigar) const {
    
    // TODO: i wish there were a way to unify the linear and affine in one traceback, but it seems pretty
    // hard actually, since the matrix can't be treated as synonymous with the CIGAR operation
    
    if (NumPW != 0) {
        uint32_t op_len = 0;
        while (true) {
            // we're in the match/mismatch matrix
            const auto& wf = wfs[s];
            if (mat == 0) {
                int64_t a = wf.int_at(d, 0);
                // TODO: i don't love this solution for the partial traceback problem...
                a -= 2 * lead_matches;
                lead_matches = 0;
                // extend through any matches
                int64_t i = (d + a) / 2 - 1;
                int64_t j = (a - d) / 2 - 1;
                int gap_close_from = 0;
                while (true) {
                    // check if this is where a gap closed
                    for (int k = 1; k <= NumPW; ++k) {
                        if (a == wf.int_at(d, -k)) {
                            gap_close_from = -k;
                            break;
                        }
                        if (a == wf.int_at(d, k)) {
                            gap_close_from = k;
                            break;
                        }
                    }
                    if (gap_close_from != 0 || i < 0 || j < 0 || seq1[i] != seq2[j]) {
                        break;
                    }
                    op_len += 1;
                    --i;
                    --j;
                    a -= 2;
                    ++lead_matches;
                }
                if (s == 0) {
                    // this condition handles the initial wf_extend, which was not preceded
                    // by a wf_next, so we skip that part of the traceback
                    break;
                }
                if (gap_close_from != 0) {
                    // this is where an insertion closed
                    if (op_len != 0) {
                        cigar.emplace_back('M', op_len);
                    }
                    op_len = 0;
                    lead_matches = 0;
                    mat = gap_close_from;
                    continue;
                }
                else if (s >= mismatch) {
                    // this must a mismatch
                    op_len += 1;
                    s -= mismatch;
                    lead_matches = 0;
                    continue;
                }
                // the traceback goes outside of the DP structure (should only happen
                // in the low-memory code paths)
                break;
            }
            else if (mat < 0) {
                // we're in a insertion matrix
                auto extend = gap_extend[-mat - 1];
                auto open = gap_open[-mat - 1];
                if (s >= extend) {
                    const auto& wf_prev = wfs[s - extend];
                    if (d - 1 >= wf_prev.diag_begin && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()) {
                        if (wf.int_at(d, mat) == wf_prev.int_at(d - 1, mat) + 1) {
                            // an insert extended here
                            s -= extend;
                            op_len += 1;
                            d -= 1;
                            continue;
                        }
                    }
                }
                if (s >= extend + open) {
                    const auto& wf_prev = wfs[s - extend - open];
                    if (d - 1 >= wf_prev.diag_begin  && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()) {
                        if (wf.int_at(d, mat) == wf_prev.int_at(d - 1, 0) + 1) {
                            // an insert opened here
                            s -= extend + open;
                            cigar.emplace_back('I', op_len + 1);
                            d -= 1;
                            mat = 0;
                            op_len = 0;
                            continue;
                        }
                    }
                }
                // the traceback goes outside of the DP structure (should only happen
                // in the low-memory code paths)
                break;
            }
            else {
                // we're in a deletion matrix
                auto extend = gap_extend[mat - 1];
                auto open = gap_open[mat - 1];
                if (s >= extend) {
                    const auto& wf_prev = wfs[s - extend];
                    if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()) {
                        if (wf.int_at(d, mat) == wf_prev.int_at(d + 1, mat) + 1) {
                            // a deletion extended here
                            s -= extend;
                            op_len += 1;
                            d += 1;
                            continue;
                        }
                    }
                }
                if (s >= extend + open) {
                    const auto& wf_prev = wfs[s - extend - open];
                    if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()){
                        if (wf.int_at(d, mat) == wf_prev.int_at(d + 1, 0) + 1) {
                            // a deletion opened her
                            s -= extend + open;
                            d += 1;
                            cigar.emplace_back('D', op_len + 1);
                            mat = 0;
                            op_len = 0;
                            continue;
                        }
                    }
                }
                // the traceback goes outside of the DP structure (should only happen
                // in the low-memory code path)
                break;
            }
        }
        
        // handle the final operation
        if (op_len != 0) {
            if (mat == 0) {
                cigar.emplace_back('M', op_len);
            }
            else if (mat < 0) {
                cigar.emplace_back('I', op_len);
            }
            else {
                cigar.emplace_back('D', op_len);
            }
        }
    }
    else {
        
        int64_t a = wfs[s].int_at(d, 0);
        lead_matches = 0;
        
        while (a != 0) {
            // we're in the match/mismatch matrix
            const auto& wf = wfs[s];
            if (s >= gap_extend[0]) {
                const auto& wf_prev = wfs[s - gap_extend[0]];
                if (d - 1 >= wf_prev.diag_begin && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()){
                    if (a == wf_prev.int_at(d - 1, 0) + 1) {
                        s -= gap_extend[0];
                        a -= 1;
                        d -= 1;
                        lead_matches = 0;
                        if (cigar.empty() || cigar.back().op != 'I') {
                            cigar.emplace_back('I', 1);
                        }
                        else {
                            ++cigar.back().len;
                        }
                        continue;
                    }
                }
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()){
                    if (a == wf_prev.int_at(d + 1, 0) + 1) {
                        s -= gap_extend[0];
                        a -= 1;
                        d += 1;
                        lead_matches = 0;
                        if (cigar.empty() || cigar.back().op != 'D') {
                            cigar.emplace_back('D', 1);
                        }
                        else {
                            ++cigar.back().len;
                        }
                        continue;
                    }
                }
            }
            
            // extend through any matches
            int64_t i = (d + a) / 2 - 1;
            int64_t j = (a - d) / 2 - 1;
    
            if (i >= 0 && j >= 0) {
                if (seq1[i] == seq2[j]) {
                    a -= 2;
                    ++lead_matches;
                    if (cigar.empty() || cigar.back().op != 'M') {
                        cigar.emplace_back('M', 1);
                    }
                    else {
                        ++cigar.back().len;
                    }
                    continue;
                }
                else if (s >= mismatch) {
                    s -= mismatch;
                    a -= 2;
                    lead_matches = 0;
                    if (cigar.empty() || cigar.back().op != 'M') {
                        cigar.emplace_back('M', 1);
                    }
                    else {
                        ++cigar.back().len;
                    }
                    continue;
                }
            }
            
            // we may need to walk the most recent match again in case we overshot
            // at the boundary of the block
            if (lead_matches) {
                cigar.back().len -= lead_matches;
                if (cigar.back().len == 0) {
                    cigar.pop_back();
                }
                lead_matches = 0;
            }
            
            // we didn't find a backtrace, so it must be outside the current
            // computed block (for low memory code paths)
            break;
        }
    }
}


template<int NumPW>
uint32_t WFAligner<NumPW>::gcd(uint32_t a, uint32_t b) const {
    // euclid's algorithm
    if (a < b) {
        std::swap(a, b);
    }
    auto r = a % b;
    if (r == 0) {
        return b;
    }
    else {
        return gcd(b, r);
    }
}

template<int NumPW>
template<typename StringType, typename IntType>
inline void WFAligner<NumPW>::find_local_opt(const StringType& seq1, const StringType& seq2,
                                             int64_t s, const Wavefront<IntType>& wf,
                                             int64_t& opt, int64_t& opt_diag, int64_t& opt_s, size_t& max_s) const {
    // TODO: this should be vectorized as well
    for (int64_t d = wf.diag_begin, n = d + wf.size(); d < n; ++d) {
        // the score of the corresponding local alignment
        auto a = wf.int_at(d, 0);
        int32_t local_s = match * a - s;
        if (local_s > opt && a > std::abs(d)) {
            // we found a new best that is inside the matrix
            opt = local_s;
            opt_diag = d;
            opt_s = s;
            // update our bound on how far out we need to look
            max_s = s + match * (2 * std::min(seq1.size(), seq2.size()) - a);
        }
    }
}

template<int NumPW>
inline int32_t WFAligner<NumPW>::convert_score(size_t len1, size_t len2, int32_t score) const {
    return (match * (len1 + len2) - score) / 2;
}

template<int NumPW>
inline std::pair<size_t, size_t> WFAligner<NumPW>::cigar_base_length(const std::vector<CIGAROp>& cigar) const {
    size_t len1 = 0, len2 = 0;
    for (const auto& c : cigar) {
        switch (c.op) {
            case 'M':
            case '=':
            case 'X':
                len1 += c.len;
                len2 += c.len;
                break;
            case 'I':
                len1 += c.len;
                break;
            case 'D':
                len2 += c.len;
                break;
                
            default:
                break;
        }
    }
    return std::make_pair(len1, len2);
}

// merge all adjacent, equivalent operations
template<int NumPW>
inline void WFAligner<NumPW>::coalesce_cigar(std::vector<CIGAROp>& cigar) const {
    
    size_t into = 0;
    for (size_t j = 1; j < cigar.size(); ++j) {
        if (cigar[j].op == cigar[into].op) {
            // merge
            cigar[into].len += cigar[j].len;
        }
        else {
            // don't merge and move the cursor
            ++into;
            if (j != into) {
                // we merged earlier, need to compact
                cigar[into] = cigar[j];
            }
        }
    }
    // remove the excess CIGAR operations
    cigar.resize(into + 1);
}

template<int NumPW>
template<typename IntType>
inline void WFAligner<NumPW>::print_wf(const WFAligner<NumPW>::Wavefront<IntType>& wf, int diag_begin, int diag_end) const {
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.size()) {
            int a = wf.int_at(d, 0);
            if (a > -100) {
                std::cerr << (int) a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.size()) {
            int a = wf.int_at(d, -1);
            if (a > -100) {
                std::cerr << (int) a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.size()) {
            int a = wf.int_at(d, 1);
            if (a > -100) {
                std::cerr << (int) a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
}

template<int NumPW>
template<typename WFVector>
void WFAligner<NumPW>::print_wfs(const WFVector& wfs) const {
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int s = 0; s < wfs.size(); ++s) {
        min_begin = std::min<int>(wfs[s].diag_begin, min_begin);
        max_end = std::max<int>(wfs[s].diag_begin + wfs[s].size(), max_end);
    }
    std::cerr << "# print wfs #" << std::endl;
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int s = 0; s < wfs.size(); ++s) {
        std::cerr << "score " << s << std::endl;
        print_wf(wfs[s], min_begin, max_end);
    }
}

template<int NumPW>
template<typename IntType>
inline void WFAligner<NumPW>::print_wfs_low_mem(const std::vector<std::pair<int32_t, WFAligner::Wavefront<IntType>>>& wf_bank,
                                                const std::deque<WFAligner<NumPW>::Wavefront<IntType>>& wf_buffer, int s) const {
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int i = 0; i < wf_bank.size(); ++i) {
        min_begin = std::min<int>(wf_bank[i].second.diag_begin, min_begin);
        max_end = std::max<int>(wf_bank[i].second.diag_begin + wf_bank[i].second.size(), max_end);
    }
    for (int i = 0; i < wf_buffer.size(); ++i) {
        min_begin = std::min<int>(wf_buffer[i].diag_begin, min_begin);
        max_end = std::max<int>(wf_buffer[i].diag_begin + wf_buffer[i].size(), max_end);
    }
    std::cerr << "### print wf bank ###" << std::endl;
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int t = 0, i = 0; !wf_bank.empty() && t <= wf_bank.back().first; ++t) {
        std::cerr << "score " << t << std::endl;
        while (t > wf_bank[i].first) {
            ++i;
        }
        if (t == wf_bank[i].first) {
            print_wf(wf_bank[i].second, min_begin, max_end);
        }
        else {
            for (int m = 0; m < 3; ++m) {
                for (int d = min_begin; d < max_end; ++d) {
                    std::cerr << "-\t";
                }
                std::cerr << std::endl;
            }
        }
    }
    std::cerr << "### print wf buffer ###" << std::endl;
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int i = 0; i < wf_buffer.size(); ++i) {
        std::cerr << "score " << s + i - wf_buffer.size() + 1 << std::endl;
        print_wf(wf_buffer[i], min_begin, max_end);
    }
}

template<int NumPW>
template<typename IntType>
inline void WFAligner<NumPW>::print_wfs_tb(const std::deque<WFAligner<NumPW>::Wavefront<IntType>>& traceback_block,
                                           const std::vector<std::pair<int32_t, WFAligner<NumPW>::Wavefront<IntType>>>& wf_bank,
                                           int i) const {

    std::cerr << "### print traceback block " << wf_bank[i].first << ":" << wf_bank[i].first + traceback_block.size() << " ###" << std::endl;
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int k = 0; k < traceback_block.size(); ++k) {
        min_begin = std::min<int>(min_begin, traceback_block[k].diag_begin);
        max_end = std::max<int>(max_end, traceback_block[k].diag_begin + traceback_block[k].size());
    }
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int k = 0; k < traceback_block.size(); ++k) {
        std::cerr << "score " << wf_bank[i].first + k << std::endl;
        print_wf(traceback_block[k], min_begin, max_end);
    }
}

template<int NumPW>
template <typename StringType, typename WFVector>
void WFAligner<NumPW>::wavefront_viz(const StringType& seq1, const StringType& seq2,
                                     const WFVector& wfs) const {

    std::vector<std::vector<int>> matrix(seq1.size() + 1, std::vector<int>(seq2.size() + 1, -1));

    std::vector<std::vector<bool>> opt(seq1.size() + 1, std::vector<bool>(seq2.size() + 1, false));


    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();

    if (wavefront_reached(wfs[wfs.size() - 1], final_diag, final_anti_diag)) {
        auto cigar = wavefront_traceback(seq1, seq2, wfs, wfs.size() - 1, final_diag);
        int i = 0, j = 0;
        opt[i][j] = true;
        for (auto& op : cigar) {
            for (int k = 0; k < op.len; ++k) {
                if (op.op == 'M' || op.op == 'I') {
                    ++i;
                }
                if (op.op == 'M' || op.op == 'D') {
                    ++j;
                }
                opt[i][j] = true;
            }
        }
    }

    for (int s = 0; s < wfs.size(); ++s) {
        const auto& wf = wfs[s];
        for (int k = 0; k < wf.size(); ++k) {
            int d = wf.diag_begin + k;
            int a = wf.int_at(d, 0);
            if (a < -2) {
                continue;
            }
            int i = (a + d) / 2 - 1;
            int j = (a - d) / 2 - 1;
            while (i >= 0 && j >= 0 && i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]
                   && a != wf.I[k] && a != wf.D[k]) {
                matrix[i + 1][j + 1] = s;
                --i;
                --j;
                a -= 2;
            }
            matrix[i + 1][j + 1] = s;
        }
    }



    std::cerr << "\t";
    for (auto c : seq2) {
        std::cerr << "\t" << c;
    }
    std::cerr << std::endl;
    for (int i = 0; i < seq1.size() + 1; ++i) {
        if (i != 0) {
            std::cerr << seq1[i - 1];
        }
        for (int j = 0; j < seq2.size() + 1; ++j) {
            std::cerr << '\t';
            if (matrix[i][j] >= 0) {
                std::cerr << matrix[i][j];
            }
            else {
                std::cerr << '.';
            }
            if (opt[i][j]) {
                std::cerr << '*';
            }
        }
        std::cerr << std::endl;
    }
}

} // end namespace wfalm


#endif /* wfa_lm.hpp */
