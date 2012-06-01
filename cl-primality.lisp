
(defpackage :cl-primality
  (:use :cl)
  (:export #:primep
           #:miller-rabin
           #:gen-prime))

(in-package :cl-primality)

(defun *-mod (n m md)
  (mod (* n m) md))

(defun expt-mod (b e md &optional (tot 1))
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

;;; Miller-Rabin algorithm

(defun miller-rabin (n &optional (chance-of-error 1d-10))
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
                          (t (and (miller-rabin-pass n (1+ (random (1- n))))
                                  (rec n (1- n-iter)))))))
               (rec n n-iter))))))

(defun miller-rabin-pass (n a)
  (declare (optimize (speed 3) (debug 0))
           (inline miller-rabin-pass))
  (labels ((decompose-val (n s)
             (cond ((or (= n 0) (oddp n)) (values n s))
                   (t (decompose-val (/ n 2) (1+ s))))))
    (multiple-value-bind (d s) (decompose-val (1- n) 0)
         (cond ((= 1 (expt-mod a d n)) n)
               ((do* ((a-loc (expt-mod a d n) (expt-mod a-loc 2 n))
                      (i 0 (1+ i))
                      (ret (= (1- n) a-loc) (= (1- n) a-loc)))
                     ((or ret (= i s)) (if (/= i s) t))) n)
               (t nil)))))

(defun gen-prime (n-bits &optional (primep-fn #'miller-rabin))
  "Generate a prime that is N-BITS long (less than 2^N-BITS).  Just try random
number of the right length until we find one that is prime (we use MILLER-RABIN
for the test here)."
  (declare (optimize (speed 3) (debug 0)))
  (let ((max (1- (expt 2 n-bits))))
    (or (funcall primep-fn (1+ (* 2 (random max))))
        (gen-prime n-bits primep-fn))))

(defun primep (n)
  "Determine if N is prime."
  (declare (inline miller-rabin)
           (type integer n))
  (miller-rabin n 1d-300))
