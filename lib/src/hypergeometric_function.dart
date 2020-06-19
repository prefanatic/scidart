/*                                                      hyp2f1.c
 *
 *      Gauss hypergeometric function   F
 *                                     2 1
 *
 *
 * SYNOPSIS:
 *
 * double a, b, c, x, y, hyp2f1();
 *
 * y = hyp2f1( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 *  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
 *                           2 1
 *
 *           inf.
 *            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
 *   =  1 +   >   -----------------------------  x   .
 *            -         c(c+1)...(c+k) (k+1)!
 *          k = 0
 *
 *  Cases addressed are
 *      Tests and escapes for negative integer a, b, or c
 *      Linear transformation if c - a or c - b negative integer
 *      Special case c = a or c = b
 *      Linear transformation for  x near +1
 *      Transformation for x < -0.5
 *      Psi function expansion if x > 0.5 and c - a - b integer
 *      Conditionally, a recurrence on c to make c-a-b > 0
 *
 *      x < -1  AMS 15.3.7 transformation applied (Travis Oliphant)
 *         valid for b,a,c,(b-a) != integer and (c-a),(c-b) != negative integer
 *
 * x >= 1 is rejected (unless special cases are present)
 *
 * The parameters a, b, c are considered to be integer
 * valued if they are within 1.0e-14 of the nearest integer
 * (1.0e-13 for IEEE arithmetic).
 *
 * ACCURACY:
 *
 *
 *               Relative error (-1 < x < 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1,7        230000      1.2e-11     5.2e-14
 *
 * Several special cases also tested with a, b, c in
 * the range -7 to 7.
 *
 * ERROR MESSAGES:
 *
 * A "partial loss of precision" message is printed if
 * the internally estimated relative error exceeds 1^-12.
 * A "singularity" message is printed on overflow or
 * in cases not addressed (such as x < -1).
 */

/*                                                      hyp2f1  */


/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
 * Converted to Dart by @Javierd
 */

import 'dart:math';

import './gamma.dart';
import './pointers.dart';
import './psi.dart';

const double EPS = 1.0e-13;
const double EPS2 = 1.0e-10;
const double ETHRESH = 1.0e-12;

const int MAX_ITERATIONS = 10000;

const double MACHEP = 1.11022302462515654042E-16;

double hyp2f1(double a, double b, double c, double x){
  double d, d1, d2, e;
  double p, q, r, s, y, ax;
  double ia, ib, ic, id;
  double t1;
  int i, aid;
  int neg_int_a = 0, neg_int_b = 0;
  int neg_int_ca_or_cb = 0;

  DoublePointer err = DoublePointer(0.0);

  ax = x.abs();
  s = 1.0 - x;
  ia = a.round().toDouble();		/* nearest integer to a */
  ib = b.round().toDouble();

  if (x == 0.0) {
    return 1.0;
  }

  d = c - a - b;
  id = d.round().toDouble();

  if ((a == 0 || b == 0) && c != 0) {
    return 1.0;
  }

  if (a <= 0 && (a - ia).abs() < EPS) {	/* a is a negative integer */
    neg_int_a = 1;
  }

  if (b <= 0 && (b - ib).abs() < EPS) {	/* b is a negative integer */
    neg_int_b = 1;
  }

  if (d <= -1 && !((d - id).abs() > EPS && s < 0)
      && !(neg_int_a != 0 || neg_int_b != 0)) {
    return pow(s, d) * hyp2f1(c - a, c - b, c, x);
  }

  if (d <= 0 && x == 1 && !(neg_int_a != 0 || neg_int_b != 0)){
    /* Return infinity*/
    return _hypdiv();
  }

  if (ax < 1.0 || x == -1.0) {
    /* 2F1(a,b;b;x) = (1-x)**(-a) */
    if ((b - c).abs() < EPS) {	/* b = c */
      if (neg_int_b != 0) {
        y = hyp2f1_neg_c_equal_bc(a, b, x);
      } else {
        y = pow(s, -a);	/* s to the -a power */
      }
      return _hypdon(y, err);
    }
    if ((a - c).abs() < EPS) {	/* a = c */
      y = pow(s, -b);	/* s to the -b power */
      return _hypdon(y, err);
    }
  }

  if (c <= 0.0) {
    ic = c.round().toDouble();		/* nearest integer to c */
    if ((c - ic).abs() < EPS) {	/* c is a negative integer */
      /* check if termination before explosion */
      if (neg_int_a != 0 && (ia > ic)){
       return _hypok(a, b, c, x, err);
      }
      if (neg_int_b != 0 && (ib > ic)){
        return _hypok(a, b, c, x, err);
      }
      return _hypdiv();
    }
  }

  if (neg_int_a != 0 || neg_int_b != 0){ 	/* function is a polynomial */
    return _hypok(a, b, c, x, err);
  }

  t1 = (b - a).abs();
  if (x < -2.0 && (t1 - t1.round()).abs() > EPS) {
    /* This transform has a pole for b-a integer, and
   * may produce large cancellation errors for |1/x| close 1
   */
    p = hyp2f1(a, 1 - c + a, 1 - b + a, 1.0 / x);
    q = hyp2f1(b, 1 - c + b, 1 - a + b, 1.0 / x);
    p *= pow(-x, -a);
    q *= pow(-x, -b);
    t1 = Gamma(c);
    s = t1 * Gamma(b - a) / (Gamma(b) * Gamma(c - a));
    y = t1 * Gamma(a - b) / (Gamma(a) * Gamma(c - b));
    return s * p + y * q;
  }else if (x < -1.0) {
    if (a.abs() < b.abs()) {
      return pow(s, -a) * hyp2f1(a, c - b, c, x / (x - 1));
    }
    else {
      return pow(s, -b) * hyp2f1(b, c - a, c, x / (x - 1));
    }
  }

  if (ax > 1.0){ /* series diverges  */
    return _hypdiv();
  }

  p = c - a;
  ia = p.round().toDouble();		/* nearest integer to c-a */
  if ((ia <= 0.0) && ((p - ia).abs() < EPS))	/* negative int c - a */
    neg_int_ca_or_cb = 1;

  r = c - b;
  ib = r.round().toDouble();		/* nearest integer to c-b */
  if ((ib <= 0.0) && ((r - ib).abs() < EPS))	/* negative int c - b */
    neg_int_ca_or_cb = 1;

  id = d.round().toDouble();		/* nearest integer to d */
  q = (d - id).abs();

  /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
   * for reporting a bug here.  */
  if ((ax - 1.0).abs() < EPS) {	/* |x| == 1.0   */
    if (x > 0.0) {
      if (neg_int_ca_or_cb != 0) {
        if (d >= 0.0)
          return _hypf(a, b, c, x, s, d, err);
        else{
          return _hypdiv();
        }
      }
      if (d <= 0.0){
        return _hypdiv();
      }
      y = Gamma(c) * Gamma(d) / (Gamma(p) * Gamma(r));
      return _hypdon(y, err);
    }
    if (d <= -1.0){
      return _hypdiv();
    }
  }


  /* Conditionally make d > 0 by recurrence on c
   * AMS55 #15.2.27
   */
  if (d < 0.0) {
    /* Try the power series first */
    y = hyt2f1(a, b, c, x, err);
    if (err.value < ETHRESH){
      return _hypdon(y, err);
    }

    /* Apply the recurrence if power series fails */
    err.value = 0.0;
    aid = 2 - id.round();
    e = c + aid;
    d2 = hyp2f1(a, b, e, x);
    d1 = hyp2f1(a, b, e + 1.0, x);
    q = a + b + 1.0;
    for (i = 0; i < aid; i++) {
      r = e - 1.0;
      y = (e * (r - (2.0 * e - q) * x) * d2 +
          (e - a) * (e - b) * x * d1) / (e * r * s);
      e = r;
      d1 = d2;
      d2 = y;
    }
    return _hypdon(y, err);
  }

  if (neg_int_ca_or_cb != 0) {
    /* negative integer c-a or c-b */
    /* goto hypf; */
    return _hypf(a, b, c, x, s, d, err);
  }

  return _hypok(a, b, c, x, err);

}

double _hypok(double a, double b, double c, double x, DoublePointer err){
  double y = hyt2f1(a, b, c, x, err);
  return _hypdon(y, err);
}

double _hypf(double a, double b, double c, double x, double s, double d, DoublePointer err){
  /* The transformation for c-a or c-b negative integer
   * AMS55 #15.3.3
   */
  double y = pow(s, d) * hys2f1(c - a, c - b, c, x, err);
  return _hypdon(y, err);
}

double _hypdon(double y, DoublePointer err){
  if (err.value > ETHRESH) {
    /* sf_error("hyp2f1", SF_ERROR_LOSS, NULL); */
    /* print( "Estimated err = ${err.value}\n"); */
  }
  return (y);
}

double _hypdiv(){
  /*print("Error overflow");*/
  return double.infinity;
}

/* Apply transformations for |x| near 1
* then call the power series
*/
double hyt2f1(double a, double b, double c, double x, DoublePointer loss){
  double p, q, r, s, t, y, w, d;
  double ax, id, d1, d2, e, y1;
  int i, aid, sign;
  DoublePointer err, err1;

  int ia, ib, neg_int_a = 0, neg_int_b = 0;

  ia = a.round();
  ib = b.round();

  if (a <= 0 && (a - ia).abs() < EPS) {	/* a is a negative integer */
    neg_int_a = 1;
  }

  if (b <= 0 && (b - ib).abs() < EPS) {	/* b is a negative integer */
    neg_int_b = 1;
  }

  err = DoublePointer(0.0);
  s = 1.0 - x;
  if (x < -0.5 && !(neg_int_a != 0 || neg_int_b != 0)) {
    if (b > a)
      y = pow(s, -a) * hys2f1(a, c - b, c, -x / s, err);

    else
      y = pow(s, -b) * hys2f1(c - a, b, c, -x / s, err);

    loss.value = err.value;
    return y;
  }

  d = c - a - b;
  id = d.round().toDouble();		/* nearest integer to d */

  if (x > 0.9 && !(neg_int_a != 0 || neg_int_b != 0)) {
    if ((d - id).abs() > EPS) {
      IntPointer sgngam = IntPointer();

      /* test for integer c-a-b */
      /* Try the power series first */
      y = hys2f1(a, b, c, x, err);
      if (err.value < ETHRESH){
        loss.value = err.value;
        return y;
      }

      /* If power series fails, then apply AMS55 #15.3.6 */
      q = hys2f1(a, b, 1.0 - d, s, err);
            sign = 1;
            w = lgam_sgn(d, sgngam);
            sign *= sgngam.value;
            w -= lgam_sgn(c-a, sgngam);
            sign *= sgngam.value;
            w -= lgam_sgn(c-b, sgngam);
            sign *= sgngam.value;
      q *= sign * exp(w);
      err1 = DoublePointer();
      r = pow(s, d) * hys2f1(c - a, c - b, d + 1.0, s, err1);
            sign = 1;
            w = lgam_sgn(-d, sgngam);
            sign *= sgngam.value;
            w -= lgam_sgn(a, sgngam);
            sign *= sgngam.value;
            w -= lgam_sgn(b, sgngam);
            sign *= sgngam.value;
      r *= sign * exp(w);
      y = q + r;

      q = q.abs();	/* estimate cancellation error */
      r = r.abs();
      if (q > r)
    r = q;
      err.value += err1.value + (MACHEP * r) / y;

      y *= Gamma(c);
      loss.value = err.value;
      return y;
    }else {
      /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12
       *
       * Although AMS55 does not explicitly state it, this expansion fails
       * for negative integer a or b, since the psi and Gamma functions
       * involved have poles.
       */

      if (id >= 0.0) {
        e = d;
        d1 = d;
        d2 = 0.0;
        aid = id.round();
      }else {
        e = -d;
        d1 = 0.0;
        d2 = d;
        aid = -id.round();
      }

      ax = log(s);

      /* sum for t = 0 */
      y = psi(1.0) + psi(1.0 + e) - psi(a + d1) - psi(b + d1) - ax;
      y /= Gamma(e + 1.0);

      p = (a + d1) * (b + d1) * s / Gamma(e + 2.0);	/* Poch for t=1 */
      t = 1.0;
      do {
        r = psi(1.0 + t) + psi(1.0 + t + e) - psi(a + t + d1)
            - psi(b + t + d1) - ax;
        q = p * r;
        y += q;
        p *= s * (a + t + d1) / (t + 1.0);
        p *= (b + t + d1) / (t + 1.0 + e);
        t += 1.0;
          if (t > MAX_ITERATIONS) {	/* should never happen */
              /* sf_error("hyp2f1", SF_ERROR_SLOW, NULL);*/
              loss.value = 1.0;
              return double.nan;
          }
      }while (y == 0 || (q / y).abs() > EPS);

      if (id == 0.0) {
        y *= Gamma(c) / (Gamma(a) * Gamma(b));
        loss.value = err.value;
        return y;
      }

      y1 = 1.0;

      if (aid == 1){
        /*goto nosum*/
        p = Gamma(c);
        y1 *= Gamma(e) * p / (Gamma(a + d1) * Gamma(b + d1));

        y *= p / (Gamma(a + d2) * Gamma(b + d2));
        if ((aid & 1) != 0){
          y = -y;
        }

        q = pow(s, id);	/* s to the id power */
        if (id > 0.0){
          y *= q;
        }else{
          y1 *= q;
        }

        y += y1;

        /* goto done */
        loss.value = err.value;
        return y;
      }

      t = 0.0;
      p = 1.0;
      for (i = 1; i < aid; i++) {
        r = 1.0 - e + t;
        p *= s * (a + t + d2) * (b + t + d2) / r;
        t += 1.0;
        p /= t;
        y1 += p;
      }

      /*no sum*/
      p = Gamma(c);
      y1 *= Gamma(e) * p / (Gamma(a + d1) * Gamma(b + d1));

      y *= p / (Gamma(a + d2) * Gamma(b + d2));
      if ((aid & 1) != 0){
        y = -y;
      }

      q = pow(s, id);	/* s to the id power */
      if (id > 0.0){
        y *= q;
      }else{
        y1 *= q;
      }

      y += y1;
      /* done */
      loss.value = err.value;
      return y;
    }

  }

  /* Use defining power series if no special cases */
  y = hys2f1(a, b, c, x, err);
  loss.value = err.value;
  return y;

}

/* Defining power series expansion of Gauss hypergeometric function */
double hys2f1(double a, double b, double c, double x, DoublePointer loss){
  double f, g, h, k, m, s, u, umax;
  int i;
  int ib, intflag = 0;

  if (b.abs() > a.abs()) {
    /* Ensure that |a| > |b| ... */
    f = b;
    b = a;
    a = f;
  }

  ib = b.round();

  if ((b - ib).abs() < EPS && ib <= 0 && b.abs() < a.abs()) {
    /* .. except when `b` is a smaller negative integer */
    f = b;
    b = a;
    a = f;
    intflag = 1;
  }

  if ((a.abs() > c.abs() + 1 || intflag != 0) && (c - a).abs() > 2
      && a.abs() > 2) {
    /* |a| >> |c| implies that large cancellation error is to be expected.
   *
   * We try to reduce it with the recurrence relations
   */
    return hyp2f1ra(a, b, c, x, loss);
  }

  i = 0;
  umax = 0.0;
  f = a;
  g = b;
  h = c;
  s = 1.0;
  u = 1.0;
  k = 0.0;

  do {
    if (h.abs() < EPS) {
      loss.value = 1.0;
      return double.infinity;
    }
    m = k + 1.0;
    u = u * ((f + k) * (g + k) * x / ((h + k) * m));
    s += u;
    k = u.abs();		/* remember largest term summed */
    if (k > umax)
    umax = k;
    k = m;
    if (++i > MAX_ITERATIONS) {	/* should never happen */
      loss.value = 1.0;
      return (s);
    }
  }while (s == 0 || (u / s).abs() > MACHEP);

  /* return estimated relative error */
  loss.value = (MACHEP * umax) / s.abs() + (MACHEP * i);
  return (s);

}

/*
 * Evaluate hypergeometric function by two-term recurrence in `a`.
 *
 * This avoids some of the loss of precision in the strongly alternating
 * hypergeometric series, and can be used to reduce the `a` and `b` parameters
 * to smaller values.
 *
 * AMS55 #15.2.10
 */
double hyp2f1ra(double a, double b, double c, double x, DoublePointer loss){
  double f2, f1, f0;
  int n;
  double t, da;
  DoublePointer err;

  /* Don't cross c or zero */
  if ((c < 0 && a <= c) || (c >= 0 && a >= c)) {
    da = (a - c).round().toDouble();
  }
  else {
    da = a.round().toDouble();
  }

  t = a - da;

  err = DoublePointer(0.0);

  assert(da != 0);

  if ((da).abs() > MAX_ITERATIONS) {
    /* Too expensive to compute this value, so give up */
    print("Too expensive to compute for this value");
    loss.value = 1.0;
    return double.nan;
  }


  if (da < 0) {
    /* Recurse down */
    f2 = 0;
    f1 = hys2f1(t, b, c, x, err);
    loss.value += err.value;
    f0 = hys2f1(t - 1, b, c, x, err);
    loss.value += err.value;
    t -= 1;
    for (n = 1; n < -da; ++n) {
      f2 = f1;
      f1 = f0;
      f0 = -(2 * t - c - t * x + b * x) / (c - t) * f1 - t * (x - 1) / (c - t) * f2;
      t -= 1;
    }
  }else {
    /* Recurse up */
    f2 = 0;
    f1 = hys2f1(t, b, c, x, err);
    loss.value += err.value;
    f0 = hys2f1(t + 1, b, c, x, err);
    loss.value += err.value;
    t += 1;
    for (n = 1; n < da; ++n) {
      f2 = f1;
      f1 = f0;
      f0 = -((2 * t - c - t * x + b * x) * f1 + (c - t) * f2) / (t * (x - 1));
      t += 1;
    }
  }

  return f0;
}

/*
  15.4.2 Abramowitz & Stegun.
*/
double hyp2f1_neg_c_equal_bc(double a, double b, double x) {
  double k;
  double collector = 1;
  double sum = 1;
  double collector_max = 1;

  if (!(b.abs() < 1e5)) {
    return double.nan;
  }

  for (k = 1; k <= -b; k++) {
    collector *= (a + k - 1)*x/k;
    collector_max = max<double>(collector.abs(), collector_max);
    sum += collector;
  }

  if (1e-16 * (1 + collector_max/sum.abs()) > 1e-7) {
    return double.nan;
  }

  return sum;
}
