import 'package:scipy_port/scipy_port.dart';

main() {
  /* Check the scipy documentation for full examples. Here you just have
   * some functions. */

  /* P(T<t) where T is a Student's t-distribution with n degrees of freedom */
  double t = 1.2341;
  int n = 3;
  print("${stdtr(n, t)}");

  /* Hypergeometric function */
  double res = hyp2f1(1, 1/2, 1, 3);
  print(res);

  /* Gamma */
  res = gamma(1/2);
  print(res);
}