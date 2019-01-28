import hts/vcf
import hts/fai
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

proc choose(v:Variant, regions: var seq[br], fai:Fai): br =
  var L = v.stop - v.start
  for i in 0..100:
    var b = random(regions)
    if b.stop - b.start > L or i == 99:
      var bL = b.stop - b.start
      var diff = int((bL - L).float64.abs / 3.0)
      if b.start + L + diff < fai.chrom_len(b.chrom):
        result = br(chrom: b.chrom, start:b.start + diff, stop: b.start + L + diff)
      else:
        result = br(chrom: b.chrom, start:b.start, stop: b.start + L)
      var fseq: string
      doAssert result.stop <= fai.chrom_len(b.chrom)
      if result.stop - result.start > 10:
        fseq = fai.get(result.chrom, result.start, result.stop)
      else:
        fseq = fai.get(result.chrom, result.start - 50, result.stop + 50)
      # keep trying if 'N' count is too high.
      if fseq.count('G').float64 != 0 and fseq.count('N').float64 / fseq.len.float64 < 0.1: return
  stderr.write_line "returning after 100 tries"

proc make_homref_variant(v:Variant, regions: var seq[br], fai:Fai): Variant =
  result = v.copy()
  var b = v.choose(regions, fai)
  #stderr.write_line "before:", v.start, "...", v.stop
  result.c.pos = b.start.int32
  result.c.rid = bcf_hdr_name2id(v.vcf.header.hdr, b.chrom)
  var s = result.tostring.replace("\t0/1:", "\t0/0:").replace("\t1/1:", "\t0/0:").strip()
  var vals = @[1'i32]
  result.from_string(result.vcf.header, s)
  doAssert result.info.set("fake", vals) == Status.OK
  var stop = @[(result.start + b.stop - b.start).int32]
  doAssert result.info.set("END", stop) == Status.OK
  #stderr.write_line "after:", result.start, "...", result.stop

proc main() =
  var
    vcf:VCF
    wtr:VCF
    fai:Fai
  if not open(vcf, paramStr(1)):
    quit "couldn't open vcf"
  var bed = read_bed(paramStr(2))
  var vs = newSeqOfCap[Variant](8192)

  doAssert vcf.header.add_info("fake", "1", "Integer", "added as a fake variant") == Status.OK

  if not open(wtr, paramStr(3), mode="w"):
    quit "bad"

  if not open(fai, paramStr(4)):
    quit "couldn't open fai at:" & paramStr(4)

  wtr.copy_header(vcf.header)

  for v in vcf:
    if v.FILTER != "PASS": continue
    doAssert v.format.fields[0].name == "GT"
    var stops = @[v.stop.int32]
    doAssert v.info.set("END", stops) == Status.OK
    vs.add(v.copy())
    vs.add(v.make_homref_variant(bed, fai))
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

