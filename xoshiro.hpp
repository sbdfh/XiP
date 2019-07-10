/**
 * \file xoshiro.hpp
 *
 * \brief Definition of the class xoshiro256
 */
#pragma once

#include <cstdint>
#include <vector>
#include <limits>
#include <iostream>

/**
 * \class xoshiro256
 *
 * \brief Base class for implementations of xoshiro256** and xoshiro256+.
 *
 * Derived classes should fulfil the standard lib named requirements for a RandomNumberEngine.
 *
 * \see http://xoshiro.di.unimi.it/
 *
 * \version 1.0
 *
 * \date 2019/07/10
 */
class xoshiro256 {

public:

  /// Type of values returned by the generator.
  /**
   * \typedef result_type
   */
  typedef uint64_t result_type;

  /// Gets the smallest possibly generated value.
  /**
   * \return The smallest value in the range of the generator (0L).
   */
  static uint64_t min() {return std::numeric_limits<uint64_t>::min();};

  /// Gets the largest possibly generated value.
  /**
   * \return The largest value in the range of the generator (18446744073709551615L).
   */
  static uint64_t max() {return std::numeric_limits<uint64_t>::max();};

  /// Default constructor.
  /**
   * Creates a xoshiro256 seeded by a default seed.
   */
  xoshiro256() : xoshiro256(default_seed) {}

  /// Copy constructor.
  /**
   * Creates a xoshiro256 with the same internal state as the provided engine.
   * \param e Engine to copy the internal state from.
   */
  xoshiro256(const xoshiro256& e) : state(e.state) {}

  /// Construct from a given seed.
  /**
   * Creates a xoshiro256 with the initial seed provided.
   * \param s The value to seed the engine with.
   */
  explicit xoshiro256(uint64_t s) {seed(s);}

  /// Construct using a given SeedSequence.
  /**
   * Creates a xoshiro256 and feeds its internal state width values from a SeedSequence.
   * \tparam SeedSequence must meet the standard lib named requirements of an SeedSequence.
   * \param sseq The SeedSequence to generate values from.
   * \see https://en.cppreference.com/w/cpp/named_req/SeedSequence
   */
  template<typename SeedSequence>
  explicit xoshiro256(SeedSequence& sseq) {seed(sseq);}

  /// Reset the engines internal state based on a seed.
  /**
   * The engine is in the same state as it would be if it was constructed using the provided seed.
   * \param seed The value to seed the engine with.
   *             If no seed is provided, then the default seed is used.
   */
  void seed(uint64_t seed = default_seed){
    state = std::vector<uint64_t>(4, seed);
  }

  /// Reset the engines internal state based on a SeedSequence.
  /**
   * The engine is in the same state as it would be if it was constructed using the provided SeedSequence.
   * \tparam SeedSequence must meet the standard lib named requirements of an SeedSequence.
   * \param sseq The SeedSequence to generate values from.
   * \see https://en.cppreference.com/w/cpp/named_req/SeedSequence
   */
  template<typename SeedSequence>
  void seed(SeedSequence& sseq){
    state = std::vector<uint64_t>(4);
    sseq.generate(state.begin(), state.end());
  }

  // Copy assignment operator.
  /**
  * Creates a xoshiro256 with the same internal state as the right hand side of the assignment.
  * \param e Engine to copy the internal state from.
   */
  xoshiro256& operator=(const xoshiro256& e){
    if (this == &e)
      return *this;
    state = e.state;
    return *this;
  }

  /**
   * We want derived classes to fulfil the standard lib named requirements for a RandomNumberEngine, hence they have to provide the () operator.
   */
  virtual uint64_t operator()() = 0;

  /// Advance the engine by a given number of steps.
  /**
   * Alias for calling operator() the same number of times and discarding the result.
   * \param z The number of steps to advance the internal state by.
   */
  void discard(unsigned long long z){
    for (unsigned long long i = 0; i < z; ++i)
      step();
  }

  /// Advance the engine by a large number of steps.
  /**
   * Alias for discard(2^128), but with constant complexity.
   */
  void jump(){
    shift_state(std::vector<uint64_t>({0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c}));
  }

  /// Advance the engine by a huge number of steps.
  /**
   * Alias for discard(2^192), but with constant complexity.
   */
  void longjump(){
    shift_state(std::vector<uint64_t>({0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635}));
  }

  /// Compares two engines
  /**
   * Two engines are equal if the they have the same internal state.
   * Engines with the same internal state produce the same sequences of numbers.
   * \param lhs,rhs The engines to compare.
   * \return true if engines are equal, false if not.
   */
  friend bool operator==(const xoshiro256& lhs, const xoshiro256& rhs){
    return (lhs.state == rhs.state);
  }

  /// Compares two engines
  /**
   * Two engines are equal if the they have the same internal state.
   * Engines with the same internal state produce the same sequences of numbers.
   * \param lhs,rhs The engines to compare.
   * \return false if engines are equal, true if not.
   */
  friend bool operator!=(const xoshiro256& lhs, const xoshiro256& rhs){
    return !(lhs == rhs);
  }

  /// Generates a textual representation of the internal state of an engine.
  /**
   * Serializes the internal state of an engine as decimal numbers separated by spaces.
   * \param os Output stream to write to.
   * \param e engine to serialize.
   * \return os
   */
  template<typename CharT, typename Traits>
  friend std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, const xoshiro256& e){
    for (auto it = e.state.begin(); it != e.state.end(); ++it)
      os << *it << " ";
    return os;
  }

  /// Tries to read an engine's internal state from a stream.
  /**
   * Tries to read decimal numbers separated by spaces from a given stream and take them as the internal state of an engine.
   * If the result of operator<< called on an engine f is used, then e is equal to f.
   * \param is Input stream to try to read from.
   * \param e engine to insert the read state into.
   * \return is
   */
  template<typename CharT, typename Traits>
  friend std::basic_istream<CharT, Traits>& operator>> (std::basic_istream<CharT, Traits>& is, xoshiro256& e){
    for (auto it = e.state.begin(); it != e.state.end(); ++it)
      is >> *it;
    return is;
  }

protected:

  std::vector<uint64_t> state;
  static const uint64_t default_seed = 0x4f9db546677ecc8d;

  uint64_t rotl (uint64_t x, uint_fast8_t k){
    return (x << k) | (x >> (64 - k));
  }

  void step(){
    uint64_t t = state[1] << 17;
    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];
    state[2] ^= t;
    state[3] = rotl(state[3], 45);
  }

  void shift_state(const std::vector<uint64_t>& jump){
    std::vector<uint64_t> s (4,0);
    for (auto it = jump.begin(); it != jump.end(); ++it){
      for (uint_fast8_t b = 0; b < 64; ++b){
        if (*it & (1 << b))
          for (uint_fast8_t i = 0; i < 4; ++i)
            s[i] ^= state[i];
        step();
      }
    }
    for (uint_fast8_t i = 0; i < 4; ++i)
      state[i] = s[i];
  }

};
