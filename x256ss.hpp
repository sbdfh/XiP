/**
 * \file x256ss.hpp
 *
 * \brief Definition of the class xoshiro256starsstar_engine
 */
#pragma once

#include "xoshiro.hpp"

/**
 * \class xoshiro256starstar_engine
 *
 * \brief A RandomNumberEngine implementation of xoshiro256**.
 *
 * Generates equidistributed bit vectors of length 64.
 * Is exceptionally fast and has a relatively large state of 256 bits.
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
class xoshiro256starstar_engine : public xoshiro256 {

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
    uint64_t ret = rotl(state[1] * 5 , 7) * 9;
    step();
    return ret;
  }
};
