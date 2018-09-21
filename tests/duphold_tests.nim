import unittest
import duphold

suite "test duphold":
  test "test start near 0":
      var values = newSeq[int32](2000)
      for i, v in values:
          values[i] = i.int32
      var r = check_rapid_depth_change[int32](0, 2, values, w=7)
      check r >= 0

  test "test end near chrom length":
      var values = newSeq[int32](2000)
      for i, v in values:
          values[i] = i.int32
      for w in 1..70:
          var r = check_rapid_depth_change[int32](values.len - 3, values.len - 1, values, w=w)
          check r >= 0

  test "test bnd position":
     check get_bnd_mate_pos("N]1:81660351]", "1") == 81660351
     check get_bnd_mate_pos("N]1:81660351]", "2") == -1
     check get_bnd_mate_pos("N]2:81660351]", "1") == -1

     check get_bnd_mate_pos("N]2:81660351]", "2") == 81660351

     check get_bnd_mate_pos("[3:81660357[N", "3") == 81660357

     check get_bnd_mate_pos("]3:110413393]N", "3") == 110413393
     check get_bnd_mate_pos("[1:64839545[N", "1") == 64839545
