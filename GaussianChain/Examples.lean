import GaussianChain.Discriminant

namespace GaussianChain

def fibSevenCertificate : Certificate where
  x := 233
  y := 1098
  z := 1597
  p := 865
  q := 499
  r := 1364
  U := 6
  V := 6
  w := 10
  c := 2
  J := 72
  h := 5
  N := 22698161
  hp := by norm_num
  hq := by norm_num
  hr := by norm_num
  hw := by norm_num
  hc := by norm_num
  hJ := by norm_num
  hN := by norm_num

def fibSevenGeomCertificate : GeomCertificate where
  toCertificate := fibSevenCertificate
  n2 := 32736
  n3num := 27144
  hn2_left := by norm_num [fibSevenCertificate, Certificate.n2Left]
  hn2_right := by norm_num [fibSevenCertificate, Certificate.n2Right]
  hn3 := by norm_num [fibSevenCertificate, Certificate.n3]
  hsquare3 := by norm_num [fibSevenCertificate]
  hJ_ne_zero := by norm_num [fibSevenCertificate]

example : fibSevenCertificate.endpointValid := by
  norm_num [Certificate.endpointValid, fibSevenCertificate]

example :
    discLeading fibSevenCertificate =
      -4 * fibSevenCertificate.y ^ 2 * fibSevenCertificate.z ^ 2 *
        fibSevenCertificate.V ^ 2 := by
  exact fibSevenGeomCertificate.discLeading_eq

end GaussianChain
