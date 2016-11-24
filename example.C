#include "rng-state.h"

#include <iostream>
#include <cmath>
#include <complex>

typedef std::complex<double> Complex;

const double PI = 3.141592653589793;

template <class T>
T sqr(const T& x)
{
  return x * x;
}

void testRngField()
{
  static const char* fname = "testRngField";
  displayln(fname);
  const RngState rs(getGlobalRngState(), fname);
  const long volume = 16 * 1024;
  double sum = 0.0;
  double sumsq = 0.0;
  for (long gindex = 0; gindex < volume; ++gindex) {
    RngState rsi = rs.newtype(gindex);
    double value = 0.0;
    for (long i = 0; i < 32; ++i) {
      value += uRandGen(rsi, sqrt(2.0)/32 + 0.05, sqrt(2.0)/32 - 0.05);
      value += gRandGen(rsi, sqrt(3.0)/32, std::sqrt(0.2));
    }
    sum += value;
    sumsq += sqr(value);
  }
  displayln(ssprintf("Expected : %24.16f", std::sqrt(2.0) + std::sqrt(3.0)));
  displayln(ssprintf("Mean     : %24.16f", sum / volume));
  displayln(ssprintf("Var      : %24.16f", std::sqrt(sumsq / volume - sqr(sum / volume)) / sqrt(volume-1)));
}

void testRngBlock()
{
  static const char* fname = "testRngBlock";
  displayln(fname);
  RngState rs(getGlobalRngState(), fname);
  const int Nb = 128;
  const int Ni = 8;
  const int Ndrop = 2;
  const int Ntake = 8;
  double sum = 0;
  double sigma2 = 0;
  for (long block = 0; block < Nb; ++block) {
    Complex a = 0;
    for (long id = 0; id < Ni; ++id) {
      long index = block * Ni + id;
      RngState rsi = rs.newtype(index);
      for (long i = 0; i < Ndrop; ++i) {
        uRandGen(rsi);
      }
      for (long i = 0; i < Ntake; ++i) {
        a += std::polar(1.0, uRandGen(rsi, PI, -PI));
      }
    }
    sum += norm(a);
    sigma2 += sqr(norm(a));
  }
  displayln(ssprintf("Expected : %24.16f", (double)Ni * Ntake));
  displayln(ssprintf("Mean     : %24.16f", sum / Nb));
  displayln(ssprintf("Var      : %24.16f", sqrt(sigma2 / Nb - sqr(sum / Nb)) / sqrt(Nb-1)));
}

void testGaussianRngProb()
{
  static const char* fname = "testGaussianRngProb";
  displayln(fname);
  RngState rs(getGlobalRngState(), fname);
  displayln(ssprintf("Expected: %24.16f", 0.682689492137));
  const long limit = 1024 * 1024;
  long count = 0;
  for (long i = 0; i < limit; ++i) {
    double x = gRandGen(rs);
    if (std::abs(x) <= 1.0) {
      count += 1;
    }
  }
  displayln(ssprintf("Find-out: %24.16f", (double)count / limit));
  count = 0;
  for (int i = 0; i < limit; ++i) {
    double x = gRandGen(rs, 1.0, 2.0);
    if (std::abs(x - 1.0) <= 2.0) {
      count += 1;
    }
  }
  displayln(ssprintf("Find-out: %24.16f", (double)count / limit));
}

void testRngApi()
{
  static const char* fname = "testRngApi";
  displayln(fname);
  RngState rs(getGlobalRngState(), fname);
  // test equal
  assert(rs == rs);
  RngState rsi;
  // test not equal
  assert(rs != rsi);
  // test split
  splitRngState(rsi, rs, "test-split");
  assert(rsi == RngState(rs, "test-split"));
  assert(rsi == rs.split("test-split"));
  splitRngState(rsi, rs, 314);
  assert(rsi == RngState(rs, 314));
  assert(rsi == RngState(rs, "314"));
  assert(rsi == rs.split(314));
  assert(rsi == rs.split("314"));
  splitRngState(rsi, rs, -314);
  assert(rsi == RngState(rs, -314));
  assert(rsi == RngState(rs, "-314"));
  assert(rsi == rs.split(-314));
  assert(rsi == rs.split("-314"));
  // test reset
  reset(rsi);
  assert(rsi == RngState());
  reset(rsi, "test-reset");
  assert(rsi == RngState("test-reset"));
  reset(rsi, 273);
  assert(rsi == RngState(273));
  assert(rsi == RngState("273"));
  reset(rsi, "273");
  assert(rsi == RngState(273));
  // test type
  rsi = rs;
  setType(rsi, 0);
  assert(rsi == rs.newtype(0));
  assert(rsi != rs);
  assert(rsi.split("type-split") != rs.split("type-split"));
  rsi = rs;
  setType(rsi, 23);
  assert(rsi == rs.newtype(23));
  assert(rsi != rs);
  assert(rsi.split("type-split") != rs.split("type-split"));
}

int main()
{
  assert(getGlobalRngState() == RngState());
  reset(getGlobalRngState(), "exampleRngState");
  assert(getGlobalRngState() == RngState("exampleRngState"));
  testRngField();
  testRngBlock();
  testGaussianRngProb();
  testRngApi();
  assert(getGlobalRngState() == RngState("exampleRngState"));
  return 0;
}
