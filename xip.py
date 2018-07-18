##
# \file xip.py
#
# \brief Python side of XiP
##

##
# \package xip
#
# \brief Contains XiP
#
# Contains the class XiP, a python interface for pseudorandom number generation.
##
from os import urandom, path
import struct
from ctypes import *
from math import pi

##
# \class XiP
#
# \brief Generates pseudorandom numbers using a xoshiro generator.
#
# This is the python side of XiP.
# Binds to a reshaped xoshiro pseudorandom number generator implemented in C++.
# Compliant with the naming scheme used by pythons random module.
# However, most functions have additional optional parameters (mostly for generating more than one sample per call).
# Check with the python documentation for details on each function.
#
# Each object gets an ID on creation and uses that ID to talk to the binding.
#
# By default, XiP uses a xoshiro256** generator, which is the best all-purpose choice.
# XiP also offers a xoshiro256+ generator.
# However, this should only be used for generation of floating point numbers.
# Always use the default generator for bit strings and integers.
#
# This class provides an extended interface of python's random module.
# This means it easily be substituted into python projects using stock randomness.
#
# \see http://xoshiro.di.unimi.it/
# \see https://docs.python.org/3/library/random.html
#
# \version 1.0
#
# \date 2018/07/12
##
class XiP:
    _streamCount = 0
    _binding = CDLL(path.dirname(path.abspath(__file__))+'/libxip.so')

    x256ss = 0
    x256p = 1

    ## \cond
    # We need doxygen to ignore this part, otherwise it thinks argtypes and restype are members of XiP.
    _binding.create_new_stream.argtypes = [c_ulonglong, c_uint]
    _binding.get_state.argtypes = [c_uint]
    _binding.get_state.restype = py_object
    _binding.set_state.argtypes = [c_uint, c_char_p]
    _binding.get_rand_bits.argtypes = [c_uint]
    _binding.get_rand_bits.restype = c_ulonglong
    _binding.uniform.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.uniform.restype = py_object
    _binding.triangular.argtypes = [c_uint, c_uint, c_double, c_double, c_double]
    _binding.triangular.restype = py_object
    _binding.gammavariate.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.gammavariate.restype = py_object
    _binding.betavariate.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.betavariate.restype = py_object
    _binding.exponential.argtypes = [c_uint, c_uint, c_double]
    _binding.exponential.restype = py_object
    _binding.gaussian.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.gaussian.restype = py_object
    _binding.lognormal.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.lognormal.restype = py_object
    _binding.vonmises.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.vonmises.restype = py_object
    _binding.pareto.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.pareto.restype = py_object
    _binding.weibull.argtypes = [c_uint, c_uint, c_double, c_double]
    _binding.weibull.restype = py_object
    _binding.uniform_int.argtypes = [c_uint, c_uint, c_longlong, c_longlong, c_longlong]
    _binding.uniform_int.restype = py_object
    _binding.weighted_with_replacement.argtypes = [c_uint, c_uint, py_object, c_uint]
    _binding.weighted_with_replacement.restype = py_object
    _binding.uniform_without_replacement.argtypes = [c_uint, c_uint, c_uint]
    _binding.uniform_without_replacement.restype = py_object
    ## \endcond

    ##
    # \param generator XiP.x256ss or XiP.x256p for xoshiro256** and xoshiro256+, respectively.
    ##
    def __init__(self, seed = None, generator = 0):
        if not (generator in range(2)):
            raise ValueError('XiP: unknown generator mode')
        self._streamID = XiP._streamCount
        XiP._streamCount += 1
        if not seed:
            XiP._binding.create_new_stream(struct.unpack('@L',urandom(8))[0], generator)
        else:
            XiP._binding.create_new_stream(hash(seed), generator)

    def getstate(self):
        return XiP._binding.get_state(self._streamID)

    def setstate(self, state):
        XiP._binding.set_state(self._streamID, state.encode('utf8'))

    def getrandbits(self, k):
        cur = 0
        while k > 0:
            if k >= 64:
                cur = (cur << 64) + XiP._binding.get_rand_bits(self._streamID)
            else:
                cur = (cur << k) + (XiP._binding.get_rand_bits(self._streamID) >> (64 - k))
            k -= 64
        return cur

    def random(self, n=1, lb=0, ub = 1):
        if n <= 0:
            raise ValueError('random: can only generate a positive number of samples')
        return XiP._binding.uniform(self._streamID, n, lb, ub)

    def uniform(self, a, b, n=1):
        if n <= 0:
            raise ValueError('uniform: can only generate a positive number of samples')
        return self.random(n, a, b)

    def triangular (self, low = 0, high = 1, mode = None, n=1):
        if n <= 0:
            raise ValueError('triangular: can only generate a positive number of samples')
        if not mode:
            mode = (high-low)/2
        if low >= high:
            raise ValueError('triangular: low has to be strictly smaller than high')
        if mode > high or mode < low:
            raise ValueError('triangular: mode has to be between low and high')
        return XiP._binding.triangular(self._streamID, n, low, high, mode)

    def gammavariate(self, alpha, beta, n=1):
        if n <= 0:
            raise ValueError('gammavariate: can only generate a positive number of samples')
        if alpha <= 0 or beta <= 0:
            raise ValueError('gammavariate: alpha and beta must be strictly larger than 0')
        return XiP._binding.gammavariate(self._streamID, n, alpha, beta)

    def betavariate(self, alpha, beta, n=1):
        if n <= 0:
            raise ValueError('betavariate: can only generate a positive number of samples')
        if alpha <= 0 or beta <= 0:
            raise ValueError('betavariate: alpha and beta must be strictly larger than 0')
        return XiP._binding.betavariate(self._streamID, n, alpha, beta)

    def expovariate(self, lambd, n=1):
        if n <= 0:
            raise ValueError('expovariate: can only generate a positive number of samples')
        if lambd == 0:
            raise ValueError('expovariate: lambda has to ne nonzero')
        return XiP._binding.exponential(self._streamID, n, lambd)

    def gauss(self, mu, sigma, n=1):
        if n <= 0:
            raise ValueError('gauss: can only generate a positive number of samples')
        if sigma <= 0:
            raise ValueError('gauss: sigma has to be larger than zero')
        return XiP._binding.gaussian(self._streamID, n, mu, sigma)

    def lognormvariate(self, mu, sigma, n=1):
        if n <= 0:
            raise ValueError('lognormvariate: can only generate a positive number of samples')
        if sigma <= 0:
            raise ValueError('lognormvariate: sigma has to be larger than zero')
        return XiP._binding.lognormal(self._streamID, n, mu, sigma)

    def vonmisesvariate(self, mu, kappa, n=1):
        if n <= 0:
            raise ValueError('vonmisesvariate: can only generate a positive number of samples')
        if abs(mu) > pi:
            raise ValueError('vonmisesvariate: mu has to be between -pi and pi')
        if kappa < 0:
            raise ValueError('vonmisesvariate: kappa has to be larger of equal to 0')
        return XiP._binding.vonmises(self._streamID, n, mu, kappa)

    def paretovariate(self, alpha, n=1, xm=1):
        if n <= 0:
            raise ValueError('paretovariate: can only generate a positive number of samples')
        if alpha <= 0:
            raise ValueError('paretovariate: alpha has to be larger than 0')
        if xm <= 0:
            raise ValueError('paretovariate: xm has to be larger than 0')
        return XiP._binding.pareto(self._streamID, n, alpha, xm)

    def weibullvariate(self, alpha, beta, n=1):
        if n <= 0:
            raise ValueError('weibullvariate: can only generate a positive number of samples')
        if alpha <= 0 or beta <= 0:
            raise ValueError('weibullvariate: alpha and beta have to be larger than 0')
        return XiP._binding.weibull(self._streamID, n, alpha, beta)

    def randint(self, a, b, n=1):
        if n <= 0:
            raise ValueError('randint: can only generate a positive number of samples')
        if b < a:
            raise ValueError('randint: a needs to be smaller or equal to b')
        return XiP._binding.uniform_int(self._streamID, n, a, 1, b-a+1)

    def randrange(self, start, stop=None, step=1, n=1):
        if n <= 0:
            raise ValueError('randint: can only generate a positive number of samples')
        if stop is None:
            if start > 0:
                return XiP._binding.uniform_int(self._streamID, n, 0, 1, start)
            raise ValueError('randint: stop needs to be larger than 0')
        width = stop - start
        if step == 1 and width > 0:
            return XiP._binding.uniform_int(self._streamID, n, start, 1, width)
        if step == 1:
            raise ValueError('randint: stop needs to be larger than start')
        if step > 0:
            w = (width + step - 1) // step
        elif step < 0:
            w = (width + step + 1) // step
        else:
            raise ValueError('randint: step cannot be zero')
        if w <= 0:
            raise ValueError('randint: empty range')
        return XiP._binding.uniform_int(self._streamID, n, start, step, w)

    def shuffle(self, x):
        if len(x) <= 1:
            return
        seq = XiP._binding.uniform_int(self._streamID, len(x), 0, 1, len(x))
        for i in range(len(x)):
            x[i], x[seq[i]] = x[seq[i]], x[i]

    def choice(self, seq, n=1):
        if len(seq) == 0:
            raise IndexError('choice: cannot choose from an empty sequence')
        I = XiP._binding.uniform_int(self._streamID,n,0,1,len(seq))
        if n == 1:
            return seq[I]
        else:
            return [seq[i] for i in I]

    def choices(self, population, weights=None,*,cum_weights=None, k=1):
        n = len(population)
        if cum_weights is None:
            if weights is None:
                if k == 1:
                    return [choice(population, 1)]
                else:
                    return choice(population,k)
        elif weights is not None:
            raise TypeError('choices: cannot specify both weights and cumulative weights')
        else:
            if len(cum_weights) != n:
                raise ValueError('choices: number of weights does not match population size')
            I = XiP._binding.weighted_with_replacement(self._streamID, k, cum_weights, 1)
        if weights is not None:
            I = XiP._binding.weighted_with_replacement(self._streamID, k, weights, 0)
        if k == 1:
            return [population[I]]
        else:
            return [population[i] for i in I]

    def sample(self, population, k):
        if not 0 <= k <= len(population):
            raise ValueError('sample: sample must be between 0 and population size')
        I = XiP._binding.uniform_without_replacement(self._streamID, len(population), k)
        if k == 1:
            return [population[I]]
        else:
            return [population[i] for i in I]
