import hts/vcf
import hts/private/hts_concat
import algorithm
import strutils
import os
import random

# make random variants of similar size distribution within the give bed file.
# all that changes is the chrom and position.

type br = object
  chrom: string
  start: int
  stop: int

proc read_bed(path: string): seq[br] =
  result = newSeqOfCap[br](2000)
  for l in path.lines:
    var toks = l.strip().split('\t')
    if toks[0] == "Y": continue
    result.add(br(chrom: toks[0], start: parseInt(toks[1]), stop: parseInt(toks[2])))

proc choose(v:Variant, regions: var seq[br]): br =
  var L = v.stop - v.start
  for i in 0..100:
    var b = rand(regions)
    if b.stop - b.start > L or i == 99:
      var bL = b.stop - b.start
      var diff = int((bL - L).float64.abs / 2.0)
      return br(chrom: b.chrom, start:b.start + diff, stop: b.start + L + diff)

proc make_homref_variant(v:Variant, regions: var seq[br]): Variant =
  result = v.copy()
  var b = v.choose(regions)
  stderr.write_line "before:", v.start, "...", v.stop
  result.c.pos = b.start.int32
  result.c.rid = bcf_hdr_name2id(v.vcf.header.hdr, v.CHROM)
  var s = result.tostring.replace("\t0/1:", "\t0/0:").replace("\t1/1:", "\t0/0:").strip()
  var vals = @[1'i32]
  result.from_string(result.vcf.header, s)
  doAssert result.info.set("fake", vals) == Status.OK
  var stop = @[(result.start + b.stop - b.start).int32]
  doAssert result.info.set("END", stop) == Status.OK
  stderr.write_line "after:", result.start, "...", result.stop

proc main() =
  var
    vcf:VCF
    wtr:VCF
  if not open(vcf, paramStr(1)):
    quit "couldn't open vcf"
  var bed = read_bed(paramStr(2))
  var vs = newSeqOfCap[Variant](8192)

  doAssert vcf.header.add_info("fake", "1", "Integer", "added as a fake variant") == Status.OK

  if not open(wtr, paramStr(3), mode="w"):
    quit "bad"

  wtr.copy_header(vcf.header)

  for v in vcf:
    if v.FILTER != "PASS": continue
    doAssert v.format.fields[0].name == "GT"
    vs.add(v.copy())
    vs.add(v.make_homref_variant(bed))
  stderr.write_line $vs.len

  vs.sort(proc(a, b: Variant): int =
    if a.rid != b.rid: return cmp(a.rid, b.rid)
    if a.start != b.start: return cmp(a.start, b.start)
    return cmp(a.stop, b.stop)
  )

  doAssert wtr.write_header
  for v in vs:
    doAssert wtr.write_variant(v)
  wtr.close


when isMainModule:
  main()

