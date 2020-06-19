# scipy_port

This Dart package include some scipy functionalities which are directly ported from its source code.

Some minors modifications have been made, including the use of IntegerPointer and DoublePointer objects to replace C pointer's where they were used. You can check the oficial documentation 

[here]:https://docs.scipy.org/doc/scipy/reference/index.html although most of them are not ported (yet?).

The included utilities are mostly from `scipy.special`and some from `scipy.stats`:

- Beta function

- Gamma function

- Hypergeometric function

- Incomplete beta integral

- Inverse of incomplete beta integral

- Inverse of Normal distribution function

- Evaluate polynomial

- Psi (digamma) function

- Student's t distribution CDF and it's inverse.

Feel free to contribute adding more functionality by making Pull requests. I will check and add them.