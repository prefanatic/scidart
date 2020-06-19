/*                                                     beta.c
 *
 *     Beta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, y, beta();
 *
 * y = beta( a, b );
 *
 *
 *
 * DESCRIPTION:
 *
 *                   -     -
 *                  | (a) | (b)
 * beta( a, b )  =  -----------.
 *                     -
 *                    | (a+b)
 *
 * For large arguments the logarithm of the function is
 * evaluated using lgam(), then exponentiated.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,30       30000       8.1e-14     1.1e-14
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * beta overflow    log(beta) > MAXLOG       0.0
 *                  a or b <0 integer        0.0
 *
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

import 'dart:math';

import './gamma.dart';
import './pointers.dart';

const double ASYMP_FACTOR = 1e6;

double beta(double a, double b) {
  double y;
  IntPointer sign = IntPointer(1);

  if (a <= 0.0) {
    if (a == a.floor()) {
      if (a == a.toInt()) {
        return betaNegint(a.toInt(), b);
      } else {
        return _overflow(sign.value);
      }
    }
  }

  if (b <= 0.0) {
    if (b == b.floor()) {
      if (b == b.toInt()) {
        return betaNegint(b.toInt(), a);
      } else {
        return _overflow(sign.value);
      }
    }
  }

  if (a.abs() < b.abs()) {
    y = a;
    a = b;
    b = y;
  }

  if (a.abs() > ASYMP_FACTOR * b.abs() && a > ASYMP_FACTOR) {
    /* Avoid loss of precision in lgam(a + b) - lgam(a) */
    y = lbetaAsymp(a, b, sign);
    return sign.value * exp(y);
  }

  y = a + b;
  if (y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM) {
    IntPointer sgngam = IntPointer();
    y = lgamSgn(y, sgngam);
    sign.value *= sgngam.value; /* keep track of the sign */
    y = lgamSgn(b, sgngam) - y;
    sign.value *= sgngam.value;
    y = lgamSgn(a, sgngam) + y;
    sign.value *= sgngam.value;
    if (y > MAXLOG) {
      return _overflow(sign.value);
    }
    return (sign.value * exp(y));
  }

  y = gamma(y);
  a = gamma(a);
  b = gamma(b);
  if (y == 0.0) return _overflow(sign.value);

  if ((a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs()) {
    y = b / y;
    y *= a;
  } else {
    y = a / y;
    y *= b;
  }

  return (y);
}

double _overflow(int sign) {
  return sign * double.infinity;
}

/* Natural log of |beta|. */

double lbeta(double a, double b) {
  double y;
  IntPointer sign = IntPointer(1);

  if (a <= 0.0) {
    if (a == a.floor()) {
      if (a == a.toInt()) {
        return lbetaNegint(a.toInt(), b);
      } else {
        return _over(sign);
      }
    }
  }

  if (b <= 0.0) {
    if (b == b.floor()) {
      if (b == b.toInt()) {
        return lbetaNegint(b.toInt(), a);
      } else {
        return _over(sign);
      }
    }
  }

  if (a.abs() < b.abs()) {
    y = a;
    a = b;
    b = y;
  }

  if (a.abs() > ASYMP_FACTOR * b.abs() && a > ASYMP_FACTOR) {
    /* Avoid loss of precision in lgam(a + b) - lgam(a) */
    y = lbetaAsymp(a, b, sign);
    return y;
  }

  y = a + b;
  if (y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM) {
    IntPointer sgngam = IntPointer();
    y = lgamSgn(y, sgngam);
    sign.value *= sgngam.value; /* keep track of the sign */
    y = lgamSgn(b, sgngam) - y;
    sign.value *= sgngam.value;
    y = lgamSgn(a, sgngam) + y;
    sign.value *= sgngam.value;
    return (y);
  }

  y = gamma(y);
  a = gamma(a);
  b = gamma(b);
  if (y == 0.0) {
    return _over(sign);
  }

  if ((a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs()) {
    y = b / y;
    y *= a;
  } else {
    y = a / y;
    y *= b;
  }

  if (y < 0) {
    y = -y;
  }

  return (log(y));
}

double _over(IntPointer sign) {
  return sign.value * double.infinity;
}

/*
 * Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
 */
double lbetaAsymp(double a, double b, IntPointer sgn) {
  double r = lgamSgn(b, sgn);
  r -= b * log(a);

  r += b * (1 - b) / (2 * a);
  r += b * (1 - b) * (1 - 2 * b) / (12 * a * a);
  r += -b * b * (1 - b) * (1 - b) / (12 * a * a * a);

  return r;
}

/*
 * Special case for a negative integer argument
 */

double betaNegint(int a, double b) {
  int sgn;
  if (b == b.toInt() && 1 - a - b > 0) {
    sgn = (b.toInt() % 2 == 0) ? 1 : -1;
    return sgn * beta(1 - a - b, b);
  } else {
    return double.infinity;
  }
}

double lbetaNegint(int a, double b) {
  double r;
  if (b == b.toInt() && 1 - a - b > 0) {
    r = lbeta(1 - a - b, b);
    return r;
  } else {
    return double.infinity;
  }
}
