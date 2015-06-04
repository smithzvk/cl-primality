
(defpackage :cl-primality
  (:use :cl :iterate)
  (:export #:primep
           #:gen-prime
           #:expt-mod))

(in-package :cl-primality)

;; @\section{Introduction}

;; @CL-Primality is a small library that can test whether integers are prime or
;; not, and perform common actions that require that knowledge.  As of now, the
;; implementation is based on the Miller-Rabin probabilistic primality
;; algorithm.  It is written with some speed considerations in mind.

;; @\section{Utilities}

;; The sort of number theoretical calculations involved in primality testing
;; typically need some support for modular arithmetic.  We introduce the
;; functions <<*-mod>> and <<expt-mod>>, which perform multiplication and
;; exponentiation modulo some number.

;;<<>>=
(defun *-mod (n m md)
  "Multiply N by M, modulo MD."
  (mod (* n m) md))

(defun expt-mod-ref (b e md)
  "A reference implementation, do not use except for testing purposes."
  (mod (expt b e) md))

;;<<>>=
(defun expt-mod (b e md &optional (tot 1))
  "Raise B to the power of E, modulo MD \(leave TOT as 1)."
  (declare (type integer e))
  (cond ((= e 0) tot)
        ((oddp e)
         (expt-mod (mod (* b b) md)
                   (ash e -1)
                   md
                   (mod (* tot b) md)))
        (t (expt-mod (mod (* b b) md)
                     (ash e -1)
                     md
                     tot))))

;; @According to Euler's theorem, $a^b \mod m = a^{b \mod \phi(m)} \mod m$,
;; where $phi(m)$ is Euler's totient function equal to the number of integers
;; coprime to $m$ in the range
;; $[1,m]$.\footnote{https://stackoverflow.com/questions/11448658/modular-exponentiation}
;; For a prime number, $\phi(m) = m$.  Euler's totient product formula says
;; that\footnote{https://en.wikipedia.org/wiki/Euler%27s_totient_function#Euler.27s_product_formula}:

;; \[
;; \phi(n) = n \prod_{p|n}\left(1 - \frac 1 p \right)
;; \]

;; ... where the product is over all integers that divide $n$ including $1$ and $p$. This requires 

;; (defun eulers-totient (n)
;;   (if (and nil (primep n))
;;       n
;;       ;; Fall back to a slow method of calculating
;;       (let ((tot 0))
;;         (iter (for i :below n)
;;           (unless (> (gcd i (mod n (if (> i 0) i n))) 1)
;;             (incf tot))
;;           (finally (return tot))))))

;; (defun %expt-mod (b e md)
;;   "Raise B to the power of E, modulo MD \(leave TOT as 1)."
;;   (declare (type integer e))
;;   (let ((e (mod e (eulers-totient md))))
;;     (expt-mod b e md)))

;; @\section{Primality Algorithms}

;; @\section{Simple Trial Division}

;; The function <<trial-division>> is just a stupid implementation that is here
;; for the sake of completion.  Trial division can be done in a much more
;; efficient way, but this seems like a waste of effort that will only serve to
;; muddy the implementation and have no impact on the asymptotic time
;; complexity.  This should never be used for any purpose, really.

;;<<>>=
(defun trial-division (n)
  "Test for primality by effectively attempting to divide N by every integer
between 2 and \(/ N 2).  This should not actually be used."
  (when (> n (expt 2 32))
    (warn "~A is such a large number that trial division will likely not finish ~
           in a timely manner" n))
  (and (not (integerp (/ n 2)))
       (iter (for i from 3 to (floor n 2) by 2)
         (never (integerp (/ n i))))
       n))

;; @\subsection{The Miller-Rabin Algorithm}

;; The Miller-Rabin algorithm is a common implementation for primality testing.
;; It performs a probabilistic check for primality that gives guarantees on the
;; maximum likelihood of a false positive (a number identified as prime when it
;; is actually composite) and never gives false negatives (a number identified
;; as composite whin it is in fact prime).  This value is set via the optional
;; parameter <chance-of-error>.  By default, <chance-of-error> is set to a very
;; small number which can slow things down if you don't need that strong of a
;; guarantee.

;;<<,2>>=
(defun miller-rabin (n &optional (chance-of-error 1d-300))
  "Miller-Rabin probabilistic primality test:

Checks if N is prime with the chance of a false positive less than
CHANCE-OF-ERROR.  This algorithm never gives false negatives."
  (declare (optimize (speed 3) (debug 0)))
  (cond ((= n 1) nil)
        ((= n 2) n)
        (t (let ((n-iter (ceiling (log chance-of-error 1/4))))
             (labels
                 ((rec (n n-iter)
                    (cond ((= n-iter 0) n)
                          (t (and (miller-rabin-pass n (1+ (random (- n 1))))
                                  (rec n (- n-iter 1)))))))
               (rec n n-iter))))))

(defun miller-rabin-pass (n a)
  "Performs one 'pass' of the Miller-Rabin primality algorithm."
  (declare (optimize (speed 3) (debug 0))
           (inline miller-rabin-pass))
  (labels ((decompose-val (n s)
             (cond ((or (= n 0) (oddp n)) (values n s))
                   (t (decompose-val (/ n 2) (1+ s))))))
    (multiple-value-bind (d s) (decompose-val (- n 1) 0)
         (cond ((= 1 (expt-mod a d n)) n)
               ((do* ((a-loc (expt-mod a d n) (expt-mod a-loc 2 n))
                      (i 0 (1+ i))
                      (ret (= (- n 1) a-loc) (= (- n 1) a-loc)))
                     ((or ret (= i s)) (if (/= i s) t))) n)
               (t nil)))))

;; @Even though the Miller-Rabin test is "probabilitic," this shouldn't scare
;; you away from using this.  With the default probability of a false positive
;; of one in $10^300$, this is, so be sure, an extremely unlikely occurrence.
;; As an illustration, as of 2012, a back of the envelope calculation tells us
;; that there have been less much less than $10^21$ CPU intructions processed in
;; the history of the world.\footnote{This is estimated assuming that computing
;; is 50 years old (age of Lisp), that the population is 7 billion people and
;; they each own a PC, and that computers run at around 2 GHz on each of its
;; four cores, and that these current trends were historically constant.  This
;; gives a good, loose upper bound.}  So, even if the Miller-Rabin algorithm was
;; a built-in instruction on the CPU, and every CPU cycle of human history has
;; been devoted to calculating prime numbers, you can be extremely certain that
;; every positive given by the Miller-Rabin algorithm (with <chance-or-error>
;; equal to even one in $10^30$ much less the default value) is prime.

;; @\section{General Interface}

;; A simple interface, <<primep>> is provided for times when you don't care
;; about how your primes are found.  This interface uses a very small chance of
;; false positive in order to be more bullet/idiot proof.

;;<<>>=
(defun primep (n)
  "Determine if N is prime."
  (declare (inline miller-rabin)
           (type integer n))
  (miller-rabin n 1d-300))

;; @\subsection{Generating primes}

;; If you need prime numbers, you can use the <<gen-prime>> function which will
;; generate a random prime number of a specified number of bits.  This uses your
;; Common Lisp's random number generator.  This means that if you really want
;; numbers that someone cannot guess, make sure that you initiallize the random
;; state on your implementation and that you trust your implementation's RNG.

;; You can specify the actual primality algorithm as the second argument.
;; However, as of now, there is only one algorithm available.

;;<<>>=
(defun gen-prime (n-bits &optional (primep-fn #'miller-rabin))
  "Generate a prime that is N-BITS long (less than 2^N-BITS).  Just try random
numbers of the right length until we find one that is prime \(we use
MILLER-RABIN for the test by default but it can be specified via PRIMEP-FN)."
  (let* ((try (random (expt 2 n-bits)))
         (odd-try (logior 1 try)))
    (or (funcall primep-fn odd-try)
        (gen-prime n-bits primep-fn))))

(defpackage :cl-primality-algorithms
  (:import-from :cl-primality
                #:miller-rabin
                #:trial-division)
  (:export #:miller-rabin
           #:trial-division))

