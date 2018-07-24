/**
 * \file xoshiro.hpp
 *
 * \brief Definition of the class xoshiro
 */
#pragma once

#include <cstdint>
#include <vector>
#include <limits>
#include <iostream>

/**
 * \class xoshiro
 *
 * \brief Implementations of xoshiro** and xoshiro+ with different state sizes.
 *
 * A RandomNumberEngine implementation of xoshiro style generators.
 *
 * Generates equidistributed bit vectors of length 64.
 * NOT CRYPTO SECURE!
 *
 * Fulfils the standard lib named requirements for a RandomNumberEngine.
 * Additionally offers a jump() and a longjump() function.
 *
 * Might fail linearity tests on the lowest bits in plus mode.
 * Should only be used by extracting the upper bits (for example for floating point generation).
 *
 * \see https://en.cppreference.com/w/cpp/named_req/RandomNumberEngine
 * \see http://xoshiro.di.unimi.it/
 *
 * \version 1.0
 *
 * \date 2018/07/23
 *
 * \see http://xoshiro.di.unimi.it/
 *
 * \version 1.0
 *
 * \date 2018/07/24
 */
class xoshiro {

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
  static uint64_t min() {return std::numeric_limits<uint64_t>::min();}

  /// Gets the largest possibly generated value.
  /**
   * \return The largest value in the range of the generator (18446744073709551615L).
   */
  static uint64_t max() {return std::numeric_limits<uint64_t>::max();}

  /// Default constructor.
  /**
   * Creates a xoshiro256** seeded by a default seed.
   */
  xoshiro() : xoshiro(default_seed, default_mode, default_size) {}

  /// Copy constructor.
 /**
  * Creates a xoshiro with the same internal state as the provided engine.
  * \param e Engine to copy the internal state from.
  */
 xoshiro(const xoshiro& e) : mode(e.mode), state_size(e.state_size), state(e.state) {}

  /// Construct using a given seed.
  /**
   * Creates a xoshiro256** and feeds its internal state width values from the provided seed.
   * \tparam SeedSequence must meet the standard lib named requirements of an SeedSequence or be an int of at most 64-bits width.
   * \param seed The SeedSequence to generate values from or the int to feed to the state.
   * \see https://en.cppreference.com/w/cpp/named_req/SeedSequence
   */
  template<typename SeedSequence>
  explicit xoshiro(SeedSequence& seed) : xoshiro(seed, default_mode, default_size) {}

  /// Construct with custom seed, genmode, and size.
  /**
   * Creates a xoshiro with the provided generation mode (** or +), state size (128, 256, or 512), and feeds its internal state width values from the provided seed.
   * \tparam SeedSequence must meet the standard lib named requirements of an SeedSequence or be an int of at most 64-bits width.
   * \param seed The SeedSequence to generate values from or the int to feed to the state.
   * \param genmode Either 0 for ** mode or 1 for + mode.
   * \param size Either 128, 256, or 512;
   * \see https://en.cppreference.com/w/cpp/named_req/SeedSequence
   */
  template<typename SeedSequence>
  explicit xoshiro(SeedSequence& s, size_t genmode, size_t size) {
    if (genmode != 0 and genmode != 1)
      throw std::logic_error("invalid genmode");
    if (size != 128 && size != 256 && size != 512)
      throw std::logic_error("invalid state size");
    mode = genmode;
    state_size = size / 64;
    state = std::vector<uint64_t> (state_size);
    seed(s);
  }

  /// Reset the engines internal state based on a seed.
  /**
   * The engine is in the same state as it would be if it was constructed using the provided seed.
   * \param seed The value to seed the engine with.
   *             If no seed is provided, then the default seed is used.
   */
  void seed(uint64_t seed = default_seed){
    for (auto it = state.begin(); it != state.end(); ++it)
      *it = seed;
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
    sseq.generate(state.begin(), state.end());
  }

  /// Generates a pseudorandom number.
  /**
   * The engines internal state is advanced by one step.
   *
   * Has constant complexity.
   * \return A pseudorandom number of the result_type from the range [min(), max()].
   */
  uint64_t operator()(){
    uint64_t ret;
    uint64_t t;
    switch (state_size){
      case 2:
        switch (mode){
          case 0:
            ret = rotl(state[0] * 5 , 7) * 9;
            break;
          case 1:
            ret = state[0] + state[1];
            break;
        }
        t = state[0] ^ state[1];
        state[0] = rotl(state[0],24) ^ t ^ (t << 16);
        state[1] = rotl(t,37);
        break;
      case 4:
        switch (mode){
          case 0:
            ret = rotl(state[1] * 5 , 7) * 9;
            break;
          case 1:
            ret = state[0] + state[3];
            break;
        }
        t = state[1] << 17;
        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];
        state[2] ^= t;
        state[3] = rotl(state[3], 45);
        break;
      case 8:
        switch (mode){
          case 0:
            ret = rotl(state[1] * 5 , 7) * 9;
            break;
          case 1:
            ret = state[0] + state[2];
            break;
        }
        t = state[1] << 11;
        state[2] ^= state[0];
        state[5] ^= state[1];
        state[1] ^= state[2];
        state[7] ^= state[3];
        state[3] ^= state[4];
        state[4] ^= state[5];
        state[0] ^= state[6];
        state[6] ^= state[7];
        state[6] ^= t;
        state[7] = rotl(state[7], 21);
        break;
    }
    return ret;
  }

  /// Advance the engine by a given number of steps.
  /**
   * Alias for calling operator() \p z number of times and discarding the result.
   * \param z The number of steps to advance the internal state by.
   */
  void discard(unsigned long long z){
    for (unsigned long long i = 0; i < z; ++i)
      (*this)();
  }

  /// Advance the engine by a large number of steps.
  /**
   * Alias for discard(2^(size of the state / 2)), but with constant complexity.
   */
  void jump(){
    switch (state_size){
      case 2:
        shift_state(std::vector<uint64_t>({0xdf900294d8f554a5, 0x170865df4b3201fc}));
        break;
      case 4:
        shift_state(std::vector<uint64_t>({0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c}));
        break;
      case 8:
        shift_state(std::vector<uint64_t>({0x33ed89b6e7a353f9, 0x760083d7955323be, 0x2837f2fbb5f22fae, 0x4b8c5674d309511c, 0xb11ac47a7ba28c25, 0xf1be7667092bcc1c, 0x53851efdb6df0aaf, 0x1ebbc8b23eaf25db}));
        break;
    }
  }

  /// Advance the engine by a huge number of steps.
  /**
   * Alias for discard(2^(size of the state * 3 / 4)), but with constant complexity.
   * Only available for state sizes 128 and 256.
   */
   void longjump(){
     switch (state_size){
      case 2:
        shift_state(std::vector<uint64_t>({0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1}));
        break;
      case 4:
        shift_state(std::vector<uint64_t>({0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635}));
        break;
      case 8:
        throw std::logic_error("longjump is not supported for state size 512.");
        break;
     }
   }

  /// Compares two engines
  /**
   * Two engines are equal if the they have the same internal state.
   * Engines with the same internal state produce the same sequences of numbers.
   * \param lhs,rhs The engines to compare.
   * \return true if engines are equal, false if not.
   */
  friend bool operator==(const xoshiro& lhs, const xoshiro& rhs){
    return (lhs.state == rhs.state);
  }

  /// Compares two engines
  /**
   * Two engines are equal if the they have the same internal state.
   * Engines with the same internal state produce the same sequences of numbers.
   * \param lhs,rhs The engines to compare.
   * \return false if engines are equal, true if not.
   */
  friend bool operator!=(const xoshiro& lhs, const xoshiro& rhs){
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
  friend std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, const xoshiro& e){
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
  friend std::basic_istream<CharT, Traits>& operator>> (std::basic_istream<CharT, Traits>& is, xoshiro& e){
    for (auto it = e.state.begin(); it != e.state.end(); ++it)
      is >> *it;
    return is;
  }

protected:

  size_t mode;
  size_t state_size;
  std::vector<uint64_t> state;

  static const uint64_t default_seed = 0x4f9db546677ecc8d;
  static const size_t default_mode = 0;
  static const size_t default_size = 256;

  static uint64_t rotl (uint64_t x, uint_fast8_t k) {
    return (x << k) | (x >> (64 - k));
  }

  void shift_state(const std::vector<uint64_t>& jump){
    std::vector<uint64_t> s (state_size,0);
    for (auto it = jump.begin(); it != jump.end(); ++it){
      for (uint_fast8_t b = 0; b < 64; ++b){
        if (*it & (1ULL << b))
          for (uint_fast8_t i = 0; i < state_size; ++i)
            s[i] ^= state[i];
        (*this)();
      }
    }
    for (uint_fast8_t i = 0; i < 4; ++i)
      state[i] = s[i];
  }

};
