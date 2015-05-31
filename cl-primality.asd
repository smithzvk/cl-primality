
(asdf:defsystem #:cl-primality
  :author "Zach Kost-Smith <zachkostsmith@gmail.com>"
  :license "LLGPL (http://opensource.franz.com/preamble.html)"
  :description "Primality testing"
  :components ((:file "cl-primality"))
  :serial t
  :depends-on (:iterate))

