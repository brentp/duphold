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
