import scipy.stats as ss
import numpy as np

MAX = 256

def tofloat(v, i):
    if i > 0:
      return "%.3g" % v

    if v == 0:
        return "0'f64"
    if v == 1:
        return "1'f64"
    return "%.3g" % v

def toint8(v, i):
    if i > 0:
        return str(v)
    return str(v) + "'i8"

row = np.zeros(MAX, dtype=np.int8)

print "# 0 if a diploid HET is more likely, 1 if a triploid (0.33/0.67) het is more likely and -1 if same"
print "# lookup is dp -> row, alt -> column"
print "const binHET: array[%d, array[%d, int8]] = [" % (MAX, MAX)
for dp in range(MAX):
    row[:] = -1
    for nalts in range(0, dp + 1):
        dh = ss.binom_test(nalts, dp, p=0.48)
        th = max(ss.binom_test(nalts, dp, p=0.32), ss.binom_test(nalts, dp, p=0.68))
        if abs(dh - th) < 0.2:
          row[nalts] = -1
        elif dh > th:
            row[nalts] = 0
        else:
            row[nalts] = 1

    print "  [" + ",".join(toint8(v, i) for i, v in enumerate(row)) + "],"

print "]"

print """
func het_lookup*(refs:int32, alts:int32, min_depth:int): int8 {.inline.} =
  var c = refs + alts
  if c > %d or c.int < min_depth:
    return -1
  return binHET[c][alts]
""" % MAX

"""
print "const binp32: array[%d, array[%d, float64]] = [" % (MAX, MAX)

print "const binp48: array[%d, array[%d, float64]] = [" % (MAX, MAX)
for dp in range(MAX):
    row[:] = 0
    for nalts in range(0, dp + 1):
        row[nalts] = ss.binom_test(nalts, dp, p=0.48)
    print dp, "  [" + ", ".join(tofloat(v, i) for i, v in enumerate(row)) + "],"
print "]"

print "const binp32: array[%d, array[%d, float64]] = [" % (MAX, MAX)
for dp in range(MAX):
    row[:] = 0
    for nalts in range(0, dp):
        if nalts / float(dp) > 0.5:
            nalts = dp - nalts
        row[nalts] = ss.binom_test(nalts, dp, p=0.32)
    print "  [" + ", ".join(tofloat(v, i) for i, v in enumerate(row)) + "],"
print "]"

print "const binp0: array[%d, array[%d, float64]] = [" % (MAX, MAX)
for dp in range(MAX):
    row[:] = 0
    for nalts in range(0, dp):
        if nalts / float(dp) > 0.5:
            nalts = dp - nalts
        row[nalts] = ss.binom_test(nalts, dp, p=0.01)
    print "  [" + ", ".join(tofloat(v, i) for i, v in enumerate(row)) + "],"
print "]"
"""
