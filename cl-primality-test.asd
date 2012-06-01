
(asdf:defsystem #:cl-primality-test
  :name "Primality testing, testing"
  :author "Zach Kost-Smith <zachkostsmith@gmail.com>"
  :license "LLGPL"
  :components ((:file "cl-primality-tests"))
  :serial t
  :depends-on (cl-primality stefil iterate))

