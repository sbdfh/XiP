/**
 * \file splitmix.hpp
 *
 * \brief Implementation of splitmix_seq
 *
 * splitmix_seq is a SeedSequence based on the SplitMix64 PRNG.
 *
 * \see https://en.cppreference.com/w/cpp/named_req/SeedSequence
 * \see http://xoshiro.di.unimi.it/
 */
#pragma once

#include <cstdint>
#include <vector>
#include <initializer_list>

/**
 * \class splitmix_seq
 *
 * \brief A SeedSequence based on SplitMix64.
 *
 * Consumes a sequence of integer-valued data und produces a requested number of unsigned 64-bit integers based on the consumed data.
 * Can be used as a randomness extender, i.e. to produce evenly distributed values from a small or biased source.
 *
 * Fulfils the standard lib named requirements for a SeedSequence.
 *
 * \see https://en.cppreference.com/w/cpp/named_req/SeedSequence
 *
 * \version 1.0
 *
 * \date 2018/07/12
 */
class splitmix_seq {

public:

  /// Default constructor.
  /**
   * Creates a splitmix_seq with an initial sequence of length zero.
   */
  splitmix_seq() : state() {}

  /// Construct from a given sequence.
  /**
   * Creates a splitmix_seq obtained by consuming the values in the range [begin,end).
   * \tparam InputIterator must meet the standard lib named requirements of an InputIterator.
   *         std::interator_traits<>::value_type must be an integer type.
   * \param begin starting point of the iterator.
   * \param end end point of the iterator.
   * \see https://en.cppreference.com/w/cpp/named_req/InputIterator
   */
  template<typename InputIterator>
  splitmix_seq(InputIterator begin, InputIterator end) : state(begin, end) {}

  /// Construct from an initializer_list.
  /**
   * Alias for splitmix_seq(li.begin(), li.end()). Enable list-initialization.
   * \tparam T must be an integer type.
   * \param li std::initializer_list of integer type objects providing an initial sequence.
   */
  template<typename T>
  splitmix_seq(std::initializer_list<T> li) : splitmix_seq(li.begin(), li.end()) {}

  /// Provides evenly distributed 64-bit values.
  /**
   * Fills the range [begin,end) with evenly distributed unsigned 64-bit integers.
   * The produced sequence depends on the values provided during construction.
   * The algorithm is based on the SplitMix64 pseudorandom number generator.
   * \tparam RandomAccessIterator must meet the standard lib named requirements of a RandomAccessIterator.
   *         std::interator_traits<>::value_type must be an unsigned integer type of width at least 64.
   * \param begin starting point of the iterator.
   * \param end end point of the iterator.
   * \see https://en.cppreference.com/w/cpp/named_req/RandomAccessIterator
   * \see http://xoshiro.di.unimi.it/
   */
  template<typename RandomAccessIterator>
  void generate(RandomAccessIterator begin, RandomAccessIterator end){
    uint64_t x = 0;
    if (state.empty())
      x = 0xd609899c48c95c4c;
    else
      for (auto it = state.begin(); it != state.end(); ++it){
        x += 0x9e3779b97f4a7c15;
        x ^= *it;
      }
    for (RandomAccessIterator it = begin; it != end; ++it){
      uint64_t z = (x += 0x9e3779b97f4a7c15);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
      z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
      *it = z ^ (z >> 31);
    }
  }

  /// The size of the initially provided sequence.
  /**
   * \return The number of values stored during construction.
   */
  std::size_t size() const {return state.size();}

  /// The initially provided sequence
  /**
   * \tparam OutputIterator must meet the standard lib named requirements of an OutputIterator.
   *         std::interator_traits<>::value_type must be an unsigned integer type of width at least 64.
   * \param dest Start of an iterator with enough space to store size() many values of unsigned 64 bit ints.
   * \see https://en.cppreference.com/w/cpp/named_req/OutputIterator
   */
  template<typename OutputIterator>
  void param(OutputIterator dest) const {std::copy(state.begin(), state.end(), dest);}

private:

  std::vector<uint64_t> state;

};
