import unittest
import duphold

suite "test duphold":

  test "test bnd position":
     check get_bnd_mate_pos("N]1:81660351]", "1") == 81660351
     check get_bnd_mate_pos("N]1:81660351]", "2") == -1
     check get_bnd_mate_pos("N]2:81660351]", "1") == -1

     check get_bnd_mate_pos("N]2:81660351]", "2") == 81660351

     check get_bnd_mate_pos("[3:81660357[N", "3") == 81660357

     check get_bnd_mate_pos("]3:110413393]N", "3") == 110413393
     check get_bnd_mate_pos("[1:64839545[N", "1") == 64839545

  test "median":

      var m = MedianStats()
      check m.median == 0
      m.add(10, true)
      check m.median == 10
      m.add(1, true)
      m.add(5, true)
      check m.median == 5

      for k in 0..100:
          m.add(2000, true)
      check m.median == m.counts.len - 1
      for k in 0..100:
          m.drop(2000, true)
      check m.median == 5

      m.clear
      check m.median == 0

      for v in m.counts:
          check v == 0
