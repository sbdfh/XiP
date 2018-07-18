/**
 * \file xip.cpp
 *
 * \brief C++ side of XiP
 */
#include "splitmix.hpp"
#include "x256ss.hpp"
#include "x256p.hpp"
#include "reshaper.hpp"
#include "Python.h"
#include <sstream>

#include <iostream>

using namespace std;

/**
 * \addtogroup XiPInterface
 * \{
 *
 * \brief Functions interfacing the C++ implementation to python.
 *
 * This is the C++ side of XiP.
 * Provides an interface python can use to talk to a reshaped xoshiro generator.
 * Uses a SplitMix64 based seeding sequence to populate the initial state.
 *
 * A shared object can only be opened once per process.
 * Thus we hold multiple streams of pseudorandomness here.
 * When requesting a new stream the python side calls create_new_stream.
 * The different streams can then be referenced via the streamID parameter which is simply an increasing counter.
 *
 * \version 1.0
 *
 * \date 2018/07/12
 */

/// Helper function converting std::vector to a python list.
template <typename T, typename PyConvFunc>
PyObject* vector_to_py(std::vector<T> variates, PyConvFunc conv_f){
  if (variates.size() > 1){
    PyObject* ret = PyList_New(variates.size());
    for (size_t i = 0; i < variates.size(); ++i)
      PyList_SET_ITEM(ret, i, conv_f(variates[i]));
    return ret;
  } else
    return conv_f(variates[0]);
}

/// Helper function converting a python list to std::vector.
template <typename T, typename PyConvFunc>
std::vector<T> py_to_vector(PyObject* list, PyConvFunc conv_f){
  if (!PyList_CheckExact(list))
    return std::vector<T>(0);
  std::vector<T> ret(PyList_Size(list));
  for (size_t i = 0; i < ret.size(); ++i)
    ret[i] = conv_f(PyList_GetItem(list, i));
  return ret;
}

extern "C"{

  /// Currently opened streams of pseudorandomness.
  vector<xoshiro256*> streams;

  /**
   * \param genmode 0 for xoshiro256**, 1 for xoshiro256+.
   */
  void create_new_stream(uint64_t seed, size_t genmode){
    splitmix_seq seedgen({seed});
    switch (genmode){
      case 0: streams.push_back(new xoshiro256starstar_engine(seedgen));
              break;
      case 1:
              streams.push_back(new xoshiro256plus_engine(seedgen));
              break;
    }
  }

  PyObject* get_state(size_t streamID){
    stringstream s;
    s << *streams.at(streamID);
    return Py_BuildValue("s", s.str().c_str());
  }

  void set_state(size_t streamID, char* state){
    stringstream s;
    s << state;
    s >> *streams.at(streamID);
  }

  uint64_t get_rand_bits(size_t streamID){
    return (*streams.at(streamID))();
  }

  PyObject* uniform(size_t streamID, size_t n, double lb, double ub){
    return vector_to_py(reshaper::uniform(streams.at(streamID), n, lb, ub), &PyFloat_FromDouble);
  }

  PyObject* triangular(size_t streamID, size_t n, double low, double high, double mode){
    return vector_to_py(reshaper::triangular(streams.at(streamID), n, low, high, mode), &PyFloat_FromDouble);
  }

  PyObject* gammavariate(size_t streamID, size_t n, double alpha, double theta){
    return vector_to_py(reshaper::gamma(streams.at(streamID), n, alpha, theta), &PyFloat_FromDouble);
  }

  PyObject* betavariate(size_t streamID, size_t n, double alpha, double beta){
    return vector_to_py(reshaper::beta(streams.at(streamID), n, alpha, beta), &PyFloat_FromDouble);
  }

  PyObject* exponential(size_t streamID, size_t n, double lambda){
    return vector_to_py(reshaper::exponential(streams.at(streamID), n, lambda), &PyFloat_FromDouble);
  }

  PyObject* gaussian(size_t streamID, size_t n, double mu, double sigma){
    return vector_to_py(reshaper::gaussian(streams.at(streamID), n, mu, sigma), &PyFloat_FromDouble);
  }

  PyObject* lognormal(size_t streamID, size_t n, double mu, double sigma){
    return vector_to_py(reshaper::lognormal(streams.at(streamID), n, mu, sigma), &PyFloat_FromDouble);
  }

  PyObject* vonmises(size_t streamID, size_t n, double mu, double kappa){
    return vector_to_py(reshaper::vonmises(streams.at(streamID), n, mu, kappa), &PyFloat_FromDouble);
  }

  PyObject* pareto(size_t streamID, size_t n, double alpha, double xm){
    return vector_to_py(reshaper::pareto(streams.at(streamID), n, alpha, xm), &PyFloat_FromDouble);
  }

  PyObject* weibull(size_t streamID, size_t n, double lambda, double k){
    return vector_to_py(reshaper::weibull(streams.at(streamID), n, lambda, k), &PyFloat_FromDouble);
  }

  PyObject* uniform_int(size_t streamID, size_t n, int64_t lb, int64_t step, int64_t width){
    return vector_to_py(reshaper::uniform_int(streams.at(streamID), n, lb, step, width), &PyLong_FromLongLong);
  }

  /**
   * \param mode If mode is 1, then the weights are assumed to be cummulative.
   */
  PyObject* weighted_with_replacement(size_t streamID, size_t n, PyObject* weights, size_t mode){
    std::vector<double> c_weights = py_to_vector<double>(weights, &PyFloat_AsDouble);
    if (mode == 1)
      for (size_t i = c_weights.size() - 1; i > 0; --i)
        c_weights[i] -= c_weights[i-1];
    double sum = 0;
    for (auto it = c_weights.begin(); it != c_weights.end(); ++it)
      sum += *it;
    for (auto it = c_weights.begin(); it != c_weights.end(); ++it)
      *it /= sum;
    return vector_to_py(reshaper::alias_sampling(streams.at(streamID), n, c_weights), &PyLong_FromLongLong);
  }

  PyObject* uniform_without_replacement(size_t streamID, size_t ub, size_t n){
    return vector_to_py(reshaper::uniform_without_replacement(streams.at(streamID), ub, n), &PyLong_FromLongLong);
  }
}

/** \} */
