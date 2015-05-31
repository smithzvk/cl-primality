
(asdf:defsystem #:cl-primality-test
  :author "Zach Kost-Smith <zachkostsmith@gmail.com>"
  :license "LLGPL (http://opensource.franz.com/preamble.html)"
  :description "CL-Primality test suite"
  :components ((:file "cl-primality-tests"))
  :serial t
  :depends-on (:cl-primality :stefil :iterate))

