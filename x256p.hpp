/**
 * \file x256p.hpp
 *
 * \brief Definition of the class xoshiro256plus_engine
 */
#pragma once

#include "xoshiro.hpp"

/**
 * \class xoshiro256plus_engine
 *
 * \brief A RandomNumberEngine implementation of xoshiro256+.
 *
 * Generates equidistributed bit vectors of length 64.
 * Might fail linearity tests on the lowest bits.
 * Should only be used by extracting the upper bits.
 * NOT CRYPTO SECURE!
 *
 * Fulfils the standard lib named requirements for a RandomNumberEngine.
 * Additionally offers a jump() and a longjump() function.
 *
 * \see https://en.cppreference.com/w/cpp/named_req/RandomNumberEngine
 * \see http://xoshiro.di.unimi.it/
 *
 * \version 1.0
 *
 * \date 2018/07/12
 */
class xoshiro256plus_engine : public xoshiro256 {

public:

  using xoshiro256::xoshiro256;

  /// Generates a pseudorandom number.
  /**
   * The engines internal state is advanced by one step.
   *
   * Has constant complexity.
   * \return A pseudorandom number of the result_type from the range [min(), max()].
   */
  uint64_t operator()() final {
    uint64_t ret = state[0] + state[3];
    step();
    return ret;
  }
};
