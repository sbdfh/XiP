/**
 * \file reshaper.hpp
 *
 * \brief Definition of several randomness reshaper functions
 */
#pragma once

#include <cmath>
#include <set>

/**
 * \addtogroup Reshaper
 * \{
 *
 * \brief Contains the namespace reshaper
 */

/**
 * \namespace reshaper
 *
 * \brief A collection of randomness reshaper functions.
 *
 * General purpose pseudorandom number reshaping.
 * These functions use a RandomNumberEngine to generate values distributed according to different probability distributions.
 *
 * All template parameter types have to meet the standard lib requirements of a RandomNumberEngine.
 *
 * \see https://en.cppreference.com/w/cpp/named_req/RandomNumberEngine
 *
 * \version 1.0
 *
 * \date 2019/07/11
 *
 * \}
 */
namespace reshaper {

  /// Generates a uniform unit double.
  /**
   * Uses the top 52 bits of a pseudorandom bitstring as the mantissa to obtain a value in [1,2), and then substracts one.
   * \param e The RandomNumberEngine to generate bits from.
   * \return A floating point number chosen uniformly from [0,1)
   */
  template<typename RandomNumberEngine>
  double uniform(RandomNumberEngine& e){
    union{uint64_t i; double d;} u;
    u.i = 0x3FFULL << 52 | e() >> 12;
    return u.d - 1.0;
  }

  /// Generates uniform doubles in a given interval.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param low Smallest values to possibly generate.
   * \param high Upper bound on the generated values.
   * \return Uniform floating point numbers from the interval [\p low, \p high).
   */
  template<typename RandomNumberEngine>
  std::vector<double> uniform(RandomNumberEngine& e, size_t n, double low, double high){
    std::vector<double> ret(n);
    for (auto it = ret.begin(); it != ret.end(); ++it){
      *it = low + (high-low)*uniform(e);
    }
    return ret;
  }

  /// Generates values according to a triangular distributions.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param low Lower limit for the generated values.
   * \param high Upper limit for the generated values.
   * \param mode Value with the highest probability.
   * \return Triangularly distributed values.
   * \see https://en.wikipedia.org/wiki/Triangular_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> triangular(RandomNumberEngine& e, size_t n, double low, double high, double mode){
    std::vector<double> ret(n);
    double f = (mode-low)/(high-low);
    double lowshift = (high-low)*(mode-low);
    double highshift = (high-low)*(high-mode);
    double u;
    for (auto it = ret.begin(); it != ret.end(); ++it){
      u = uniform(e);
      *it = u < f ? low + sqrt(u*lowshift) : high - sqrt((1-u)*highshift);
    }
    return ret;
  }

  /// Generates values according to an exponential distribution.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param lambda The rate parameter.
   * \return Exponentially distributed values.
   * \see https://en.wikipedia.org/wiki/Exponential_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> exponential(RandomNumberEngine& e, size_t n, double lambda){
    std::vector<double> ret(n);
    for (auto it = ret.begin(); it != ret.end(); ++it)
      *it = -log(uniform(e))/lambda;
    return ret;
  }

  /// Generates two normally distributed values.
  /**
   * The most efficient way to generate Gaussian variates is to generate them in pairs.
   * This helper function generates two independent variates and writes them to the location of the \p o and \p t parameters.
   *
   * \param e The RandomNumberEngine to generate bits from.
   * \param o Reference to the first value.
   * \param t Reference to the second value.
   * \param mu The location parameter.
   * \param sigma the scale parameter.
   * \see https://en.wikipedia.org/wiki/Normal_distribution
   */
  template<typename RandomNumberEngine>
  void gen_pair_of_gaussians(RandomNumberEngine& e, double& o, double& t, double mu, double sigma){
    double u,v,r;
    do{
      u = -1 + 2*uniform(e);
      v = -1 + 2*uniform(e);
      r = u*u + v*v;
    } while (r > 1 || r == 0);

    o = mu+sigma*u*sqrt(-2*log(r)/r);
    t = mu+sigma*v*sqrt(-2*log(r)/r);
  }

  /**
   * Gaussians are generated in pairs. Hence, if the user requests an uneven number, we can save one generated variate for later.
   * Is true if new variates need to be generated and false if there is a normal variate in store.
   */
  bool _gen_gaussian = true;
  /**
   * Gaussians are generated in pairs. Hence, if the user requests an uneven number, we can save one generated variate for later.
   * If there is a variate left over from the last generation, then this holds a standard normal variate.
   */
  double _saved_gaussian;
  /// Generates values according to a Gaussian normal distribution.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param mu The location parameter.
   * \param sigma the scale parameter.
   * \return Normally distributed values.
   * \see https://en.wikipedia.org/wiki/Normal_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> gaussian(RandomNumberEngine& e, size_t n, double mu, double sigma){
    std::vector<double> ret(n);
    size_t i = 0;

    if (!_gen_gaussian){
      ret[i++] = mu+sigma*_saved_gaussian;
      _gen_gaussian = true;
    }

    for(; i < n-1 ; i += 2)
      gen_pair_of_gaussians(e, ret[i], ret[i+1], mu, sigma);
    if (i < n){
      gen_pair_of_gaussians(e, ret[i], _saved_gaussian, mu, sigma);
      _saved_gaussian = (_saved_gaussian-mu)/sigma;
      _gen_gaussian = false;
    }
    return ret;
  }

  /// Generates values according to a gamma distribution.
  /**
   * Variates are generated by using Marsaglia and Tsang's method.
   *
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param alpha The shape parameter.
   * \param beta The rate parameter.
   * \return Gamma-distributed values.
   * \see https://en.wikipedia.org/wiki/Gamma_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> gamma(RandomNumberEngine& e, size_t n, double alpha, double beta){
    std::vector<double> ret(n);
    if (alpha >= 1.0){
      double d = alpha - 1.0/3;
      double c = 1.0/sqrt(9*d);
      for (auto it = ret.begin(); it != ret.end(); ++it){
        double V,Z, t;
        while(1){
          if (_gen_gaussian){
            gen_pair_of_gaussians(e, Z, _saved_gaussian, 0, 1);
            _gen_gaussian = false;
          } else {
            Z = _saved_gaussian;
            _gen_gaussian = true;
          }
          V = 1+c*Z;
          if (V >= 0){
            V *= V*V;
            t = d*V;
            double U = uniform(e);
            if (log(U) < 0.5*Z*Z+d-t+d*log(V))
              break;
          }
        }
        *it = t/beta;
      }
    } else {
      ret = gamma(e,n,alpha+1,beta);
      std::vector<double> u = uniform(e,n,0,1);
      for(auto it = ret.begin(), uit=u.begin(); it != ret.end();++it, ++uit)
        *it *= pow(*uit,1/alpha);
    }
    return ret;
  }

  /// Generates values according to a beta distribution.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param alpha The first shape parameter.
   * \param beta The second shape parameter.
   * \return Beta-distributed values.
   * \see https://en.wikipedia.org/wiki/Beta_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> beta(RandomNumberEngine& e, size_t n, double alpha, double beta){
    std::vector<double> ret = gamma(e, n, alpha, 1);
    std::vector<double> u = gamma(e, n, beta, 1);
    for (auto it = ret.begin(), uit = u.begin(); it != ret.end(); ++it, ++uit)
      *it /= *it + *uit;
    return ret;
  }

  /// Generates values according to a lognormal distribution.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param mu The location parameter.
   * \param sigma the scale parameter.
   * \return Log-Normally distributed values.
   * \see https://en.wikipedia.org/wiki/Log-normal_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> lognormal(RandomNumberEngine& e, size_t n, double mu, double sigma){
    std::vector<double> ret = gaussian(e, n, mu, sigma);
    for (auto it = ret.begin(); it != ret.end(); ++it)
      *it = exp(*it);
    return ret;
  }

  /// Generates values according to a vonmises distribution.
  /**
   * Variates are generated using the method proposed by Barabesi.
   *
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param mu The location parameter.
   * \param kappa The measure of concentration.
   * \return Von-Mises-distributed values.
   * \see https://en.wikipedia.org/wiki/Von_Mises_distribution
   * \see http://sa-ijas.stat.unipd.it/sites/sa-ijas.stat.unipd.it/files/417-426.pdf
   */
  template<typename RandomNumberEngine>
  std::vector<double> vonmises (RandomNumberEngine& e, size_t n, double mu, double kappa){
    std::vector<double> ret(n);
    double s = kappa > 1.3 ? 1/sqrt(kappa) : M_PI*exp(-kappa);
    for (auto it = ret.begin(); it != ret.end(); ++it){
      double t;
      while (1){
        double u1 = uniform(e), u2 = uniform(e);
        t = s*(2*u2 - 1) / u1;
        if (fabs(t) <= M_PI && (kappa*t*t < 4 - 4*u1 || kappa*cos(t) >= 2*log(u1) + kappa))
          break;
      }
      *it = atan2(sin(t + mu), cos(t + mu));
    }
    return ret;
  }

  /// Generates values according to a pareto distribution.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param alpha The shape parameter.
   * \param xm The scale parameter.
   * \return Pareto-distributed values.
   * \see https://en.wikipedia.org/wiki/Pareto_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> pareto(RandomNumberEngine& e, size_t n, double alpha, double xm){
    std::vector<double> ret(n);
    double alphainv = 1 / alpha;
    for (auto it = ret.begin(); it != ret.end(); ++it)
      *it = xm / pow(uniform(e),alphainv);
    return ret;
  }

  /// Generates values according to a weibull distribution.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param lambda The scale parameter.
   * \param k The shape parameter.
   * \return Weibull-distributed values.
   * \see https://en.wikipedia.org/wiki/Weibull_distribution
   */
  template<typename RandomNumberEngine>
  std::vector<double> weibull(RandomNumberEngine& e, size_t n, double lambda, double k){
    std::vector<double> ret(n);
    double kinv = 1 / k;
    for (auto it = ret.begin(); it != ret.end(); ++it)
      *it = lambda * pow(-log(uniform(e)),kinv);
    return ret;
  }

  /// Generates a uniformly distributed integer.
  /**
   * \param e The RandomNumberEngine to generate bits from.
   * \param ub Upper bound on the generated value.
   * \return A uniform integer from the interval [0,\p ub).
   */
  template<typename RandomNumberEngine>
  uint64_t uniform_uint(RandomNumberEngine& e, uint64_t ub){
    uint_fast8_t bits = floor(log2(ub)) + 1;
    uint64_t ret;
    do {
      ret = e() >> (64 - bits);
    } while (ret >= ub);
    return ret;
  }

  /// Generates uniformly distributed integers.
  /**
   * Values are chosen from [\p lb, \p lb + \p step, \p lb + 2*\p step, ... , \p lb + \p width).
   *
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param lb Lower bound on the generated values.
   * \param step Distance between two possible values.
   * \param width Largest distance a value can have to \p lb.
   * \return Uniformly distributed integers.
   */
  template<typename RandomNumberEngine>
  std::vector<int64_t> uniform_int(RandomNumberEngine& e, size_t n, int64_t lb, int64_t step, int64_t width){
    std::vector<int64_t> ret(n);
    for (auto it = ret.begin(); it != ret.end(); ++it)
      *it = lb + step*uniform_uint(e,width);
    return ret;
  }

  /// Generates integers according to a given distribution.
  /**
   * Given some weighting, generates integers from the range [0, |\p weights|) such that the probability for each i is \p weights[i]/sum(\p weights).
   *
   * \param e The RandomNumberEngine to generate bits from.
   * \param n The number of variates to generate.
   * \param weights Individual weighting of the events.
   * \see https://en.wikipedia.org/wiki/Alias_method
   * \see http://www.keithschwarz.com/darts-dice-coins/
   */
  template<typename RandomNumberEngine>
  std::vector<uint64_t> alias_sampling(RandomNumberEngine& e, size_t n, std::vector<double> weights){
    size_t len = weights.size();
    std::vector<double> prob_table(len);
    std::vector<size_t> alias_table(len);
    std::vector<double> probabilities(len);
    std::vector<size_t> small,large;

    for (size_t i = 0; i < len; ++i){
      probabilities[i] = weights[i] * len;
      if (probabilities[i] >= 1)
        large.push_back(i);
      else
        small.push_back(i);
    }

    while (!small.empty() && !large.empty()){
      size_t cursmall = small.back(), curlarge = large.back();
      small.pop_back();
      large.pop_back();
      prob_table[cursmall] = probabilities[cursmall];
      alias_table[cursmall] = curlarge;
      probabilities[curlarge] += probabilities[cursmall] - 1;
      if (probabilities[curlarge] >= 1)
        large.push_back(curlarge);
      else
        small.push_back(curlarge);
    }

    while (!small.empty()){
      prob_table[small.back()] = 1;
      small.pop_back();
    }
    while (!large.empty()){
      prob_table[large.back()] = 1;
      large.pop_back();
    }

    std::vector<uint64_t> ret(n);
    for (auto it = ret.begin(); it != ret.end(); ++it){
      uint64_t sample = uniform_uint(e,len);
      *it = uniform(e) < prob_table[sample] ? sample : alias_table[sample];
    }
    return ret;
  }

  /// Generates a sequence of unique integers.
  /**
   * Draws \p n values without replacement from [0,ub).
   *
   * \param e The RandomNumberEngine to generate bits from.
   * \param ub Upper bound on the generated values.
   * \param n Number of values to generate.
   * \return Sequence of integers drawn without replacement.
   */
  template<typename RandomNumberEngine>
  std::vector<uint64_t> uniform_without_replacement(RandomNumberEngine& e, size_t ub, size_t n){
    std::vector<uint64_t> ret(n);
    std::set<size_t> seen;
    auto it = ret.begin();
    while (seen.size() < n){
      uint64_t cur = uniform_uint(e,ub);
      if (seen.insert(cur).second)
        *it++ = cur;
    }
    return ret;
  }
}
