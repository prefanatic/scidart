/*
 *     Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, Gamma();
 *
 * y = Gamma( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Gamma function of the argument.  The result is
 * correctly signed.
 *
 * Arguments |x| <= 34 are reduced by recurrence and the function
 * approximated by a rational function of degree 6/7 in the
 * interval (2,3).  Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE    -170,-33      20000       2.3e-15     3.3e-16
 *    IEEE     -33,  33     20000       9.4e-16     2.2e-16
 *    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 */

/*                                                     lgam()
 *
 *     Natural logarithm of Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgam();
 *
 * y = lgam( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the Gamma function of the argument.
 *
 * For arguments greater than 13, the logarithm of the Gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGM return NPY_INFINITY and an error
 * message.  MAXLGM = 2.556348e305 for IEEE arithmetic.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic      domain        # trials     peak         rms
 *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
 *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
 *
 */

/*
 * Cephes Math Library Release 2.2:  July, 1992
 * Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 * Converted to Dart by @Javierd
 */

import 'dart:math';

import './pointers.dart';
import './polevl.dart';

const List<double> P = [
  1.60119522476751861407E-4,
  1.19135147006586384913E-3,
  1.04213797561761569935E-2,
  4.76367800457137231464E-2,
  2.07448227648435975150E-1,
  4.94214826801497100753E-1,
  9.99999999999999996796E-1
];

const List<double> Q = [
  -2.31581873324120129819E-5,
  5.39605580493303397842E-4,
  -4.45641913851797240494E-3,
  1.18139785222060435552E-2,
  3.58236398605498653373E-2,
  -2.34591795718243348568E-1,
  7.14304917030273074085E-2,
  1.00000000000000000320E0
];

const double MAXGAM = 171.624376956302725;
const double LOGPI = 1.14472988584940017414;

const List<double> STIR = [
  7.87311395793093628397E-4,
  -2.29549961613378126380E-4,
  -2.68132617805781232825E-3,
  3.47222221605458667310E-3,
  8.33333333333482257126E-2,
];

const double MAXSTIR = 143.01608;
const double SQTPI = 2.50662827463100050242E0;

const double MAXLOG = 7.09782712893383996732E2;
const double MAXLGM = 2.556348e305;

/* Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
double stirf(double x){
  double y, w, v;

  if (x >= MAXGAM) {
    return double.infinity;
  }
  w = 1.0 / x;
  w = 1.0 + w * polevl(w, STIR, 4);
  y = exp(x);
  if (x > MAXSTIR) {		/* Avoid overflow in pow() */
    v = pow(x, 0.5 * x - 0.25) as double;
    y = v * (v / y);
  }
  else {
    y = pow(x, x - 0.5) / y;
  }
  y = SQTPI * y * w;
  return (y);
}

double gamma(double x){
  double p, q, z;
  int i;
  int sgngam = 1;

  if (!x.isFinite) {
    return x;
  }
  q = x.abs();

  if (q > 33.0) {
    if (x < 0.0) {
      p = q.floor().toDouble();
      if (p == q) {
        return _gamnan();
      }
      i = p.round();
      if ((i & 1) == 0)
        sgngam = -1;
      z = q - p;
      if (z > 0.5) {
        p += 1.0;
        z = q - p;
      }
      z = q * sin(pi * z);
      if (z == 0.0) {
        return (sgngam * double.infinity);
      }
      z = z.abs();
      z = pi / (z * stirf(q));
    }
    else {
      z = stirf(x);
    }
    return (sgngam * z);
  }

  z = 1.0;
  while (x >= 3.0) {
    x -= 1.0;
    z *= x;
  }

  while (x < 0.0) {
    if (x > -1E-9){
      return _small(x, z);
    }
    z /= x;
    x += 1.0;
  }

  while (x < 2.0) {
    if (x < 1e-9){
      return _small(x, z);
    }

    z /= x;
    x += 1.0;
  }

  if (x == 2.0)
    return (z);

  x -= 2.0;
  p = polevl(x, P, 6);
  q = polevl(x, Q, 7);
  return (z * p / q);
}

double _gamnan(){
  /* sf_error("Gamma", SF_ERROR_OVERFLOW, NULL); */
  return double.infinity;
}

double _small(double x, double z){
  if (x == 0.0) {
    return _gamnan();
  }
  else
    return (z / ((1.0 + 0.5772156649015329 * x) * x));
}

/* A[]: Stirling's formula expansion of log Gamma
* B[], C[]: log Gamma function between 2 and 3
*/
const List<double> A = [
  8.11614167470508450300E-4,
  -5.95061904284301438324E-4,
  7.93650340457716943945E-4,
  -2.77777777730099687205E-3,
  8.33333333333331927722E-2
];

const List<double> B = [
  -1.37825152569120859100E3,
  -3.88016315134637840924E4,
  -3.31612992738871184744E5,
  -1.16237097492762307383E6,
  -1.72173700820839662146E6,
  -8.53555664245765465627E5
];

const List<double> C = [
/* 1.00000000000000000000E0, */
  -3.51815701436523470549E2,
  -1.70642106651881159223E4,
  -2.20528590553854454839E5,
  -1.13933444367982507207E6,
  -2.53252307177582951285E6,
  -2.01889141433532773231E6
];

/* log(sqrt(2*pi)) */
const double LS2PI = 0.91893853320467274178;

/* Logarithm of Gamma function */
double lgam(double x) {
  IntPointer sign = IntPointer();
  return lgamSgn(x, sign);
}

double lgamSgn(double x, IntPointer sign){
  double p, q, u, w, z;
  int i;

  sign.value = 1;

  if (!x.isFinite){
    return x;
  }

  if (x < -34.0) {
    q = -x;
    w = lgamSgn(q, sign);
    p = q.floor().toDouble();
    if (p == q) {
      return double.infinity;
    }

    i = p.round();
    if ((i & 1) == 0){
      sign.value = -1;
    }else{
      sign.value = 1;
    }

    z = q - p;
    if (z > 0.5) {
      p += 1.0;
      z = p - q;
    }

    z = q * sin(pi * z);
    if (z == 0.0){
      return double.infinity;
    }

    /*     z = log(NPY_PI) - log( z ) - w; */
    z = LOGPI - log(z) - w;
    return (z);
  }

  if (x < 13.0) {
    z = 1.0;
    p = 0.0;
    u = x;
    while (u >= 3.0) {
      p -= 1.0;
      u = x + p;
      z *= u;
    }
    while (u < 2.0) {
      if (u == 0.0){
        return double.infinity;
      }

      z /= u;
      p += 1.0;
      u = x + p;
    }
    if (z < 0.0) {
      sign.value = -1;
      z = -z;
    }else{
      sign.value = 1;
    }

    if (u == 2.0){
      return (log(z));
    }

    p -= 2.0;
    x = x + p;
    p = x * polevl(x, B, 5) / p1evl(x, C, 6);
    return (log(z) + p);
  }

  if (x > MAXLGM) {
    return (sign.value * double.infinity);
  }

  q = (x - 0.5) * log(x) - x + LS2PI;
  if (x > 1.0e8)
    return (q);

  p = 1.0 / (x * x);
  if (x >= 1000.0){
    q += ((7.9365079365079365079365e-4 * p
        - 2.7777777777777777777778e-3) * p
        + 0.0833333333333333333333) / x;
  }else{
    q += polevl(p, A, 4) / x;
  }

  return (q);

}