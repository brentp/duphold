import hts
import times
import algorithm
import random
import strutils
import logging
import strformat
import os
import math
import docopt
import genoiser
import ./dupholdpkg/version

var STEP = 250 ## size of GC window
var SIZE = 1_000 # flank size
const GC_BINS = 20 ## number of bins to partition gc-content

proc setenvvars() =
  try:
    SIZE = parseInt(getEnv("DUPHOLD_FLANK"))
    stderr.write_line "[duphold] set flank size to " & $SIZE
  except:
    discard
  try:
    STEP = parseInt(getEnv("DUPHOLD_GC_STEP"))
    stderr.write_line "[duphold] set gc-step to " & $STEP
  except:
    discard


type MedianStats* = object
  ## since we expect nearly all values to be < 1000
  ## we can store counts of observed values in a fixed-size array
  ## and get median (or any percentile) from those.
  ## we can even keep a moving window of medians as we can drop
  ## values by decrementing the counts.
  counts*: array[16384, int]
  i_median: int
  n*: int
  i_ok: bool

proc add*(m:var MedianStats, d:int) {.inline.} =
    if d < m.counts.high:
      m.counts[d] += 1
    else:
      m.counts[m.counts.high] += 1
    m.n += 1
    m.i_ok = false

proc drop*(m:var MedianStats, d:int) {.inline.} =
    if d < m.counts.high:
      m.counts[d] -= 1
    else:
      m.counts[m.counts.high] -= 1
    m.n -= 1
    m.i_ok = false

proc median*(m:var MedianStats, skip_zeros:bool=false): int {.inline.} =
    if m.i_ok:
      return m.i_median

    var cum = 0
    var stop_n = (0.5 + m.n.float64 * 0.5).int
    if skip_zeros:
      stop_n = (0.5 + float64(m.n - m.counts[0]) * 0.5).int

    for i, cnt in m.counts:
        if skip_zeros and i == 0: continue
        cum += cnt
        if cum >= stop_n:
            m.i_median = i
            m.i_ok = true
            break
    return m.i_median

proc mean*(m:var MedianStats): int {.inline.} =
    for i, cnt in m.counts:
        result += (i * cnt)
    result = int(0.5 + result.float64 / m.n.float64)

proc clear*(m:var MedianStats) =
    zeroMem(m.counts[0].addr, sizeof(m.counts[0]) * m.counts.len)
    m.i_ok = false
    m.i_median = 0
    m.n = 0

proc percentile*(m:MedianStats, pct:float64): int {.inline.} =
  ## return the value at the requested percentile.
  var cum = 0
  var p = pct
  if p > 1: p /= 100.0
  if p > 1:
      quit "can't get value outside of distribution (> 1)"

  var stop_n = (m.n.float64 * p).int
  for i, cnt in m.counts:
        cum += cnt
        if cum >= stop_n:
            return i
  return -1

proc `$`(m:var MedianStats): string =
    return &"MedianStats(n:{m.n}, vals: {m.counts[0..1000]})"

type Discordant* = ref object
  left*: uint32
  right*: uint32

proc `$`*(d:Discordant): string =
    return &"Discordant(left: {d.left}, right: {d.right}"

proc find(targets:seq[Target], chrom:string): int =
    for i, t in targets:
        if t.name == chrom: return i
    return -1

proc count_gc(gc: var seq[bool], s:int64, e:int64): float32 {.inline.} =
    var tot = 0
    for i in s..<min(e, gc.len):
        tot += gc[i].int
    return tot.float32 / max(1, min(e, gc.len) - s).float32

proc gc_content*(fai:Fai, chrom:string, step:int, gc_bool: var seq[bool]): seq[float32] =

    var s = fai.cget(chrom)
    defer:
        free(s)
    if gc_bool.len != s.len:
        info("setting gc_bool length")
        gc_bool.set_len(s.len)

    result = newSeq[float32]((s.len/step+1).int)

    for i, c in s:
        if c == 'C' or c == 'G' or c == 'c' or c == 'g':
            result[(i/step).int] += 1
            gc_bool[i] = true
        elif c == 'N' or c == 'n':
            result[(i/step).int] -= 1
    for i, v in result:
        result[i] = v / step.float32

proc get_or_empty[T](variant:Variant, field:string, input:var seq[T], nper:int=1) {.inline.} =
  ## if, for example we've already annotated a sample in the VCF with duphold
  ## we dont want to overwite those values with nan so try to grab existing
  ## values but otherwise make an empty array.
  if variant.format.get(field, input) == Status.OK:
      return

  if input.len != variant.vcf.n_samples * nper:
    input.set_len(variant.vcf.n_samples * nper)
  for i, f in input:
    when T is float32:
        input[i] = cast[T](bcf_float_missing)
    elif T is int32:
        input[i] = int32.low
    else:
        quit "unknown type in get_or_empty"


proc get_bnd_mate_pos*(a:string, vchrom:string): int {.inline.} =
    if not (':' in a): return -1
    var tmp = a.split(':')
    var left = tmp[0]
    var i = 0
    while i < 3:
      if left[i] in {'[', ']'}:
        break
      i += 1

    var chrom = left[i+1..left.high]
    if chrom != vchrom: return -1

    var right = tmp[1]
    i = 0
    while right[i].isdigit:
        i += 1
    result = parseInt(right[0..<i])

proc count(discordants:var seq[Discordant], s:int, e:int, i99:int, slop:int=25): int =
    ## count the number of discordants that span s..e but are within d of both.
    var ilow: int
    if s <= i99:
      ilow = 1
    else:
      ilow = lowerBound(discordants, Discordant(left:uint32(s - i99), right: uint32(s - i99)), proc(a, b: Discordant):int =
        if a.left == b.left:
            return a.right.int - b.right.int
        return a.left.int - b.left.int
      )
    for i in ilow..discordants.high:
        var disc = discordants[i]
        if disc.left.int > s + i99: break
        if s - disc.left.int > i99: continue

        if disc.right.int - (e - slop) < i99:
          # after correcting for the event, the insert should fall in the expected distribution.
          #echo &"s-e:{s}-{e} supported by: {disc} with corrected span: {(disc.right.int - disc.left.int) - (e - s)}"
          var corrected = (disc.right.int - disc.left.int) - (e - s)
          if corrected < 2 * i99 and corrected > -2*slop:
            result += 1

proc get_bnd_mate_pos(variant:Variant): int {.inline.} =
    return get_bnd_mate_pos(variant.ALT[0], $variant.CHROM)

proc duphold*[T](variant:Variant, values:var seq[T], sample_i: int, stats:var MedianStats, gc_stats:var seq[MedianStats], gc_bool:var seq[bool], w:int=7, discordants:var seq[Discordant], i99:int): float64 =
    ## sets FORMAT fields for sample i in the variant and returns the DHBFC value
    var
      s = variant.start
      e = variant.stop
      isBnd = false


    if variant.ALT[0][0] != '<':
      var bnd = get_bnd_mate_pos(variant)
      var isBnd = false
      if bnd != -1:
        if bnd < s and s - bnd < 20000000:
            e = s
            s = bnd
        elif bnd > e and bnd - e < 20000000:
            s = e
            e = bnd - 1
        isBnd = true

      if bnd == -1:
          # skip distant BND's as this is not informative
          if ':' in variant.ALT[0]: return -1

    if e - s < 10:
      s = max(0, s - 50)
      e = min(e + 50, values.high)
    e = min(e, values.high)
    if s > e:
      stderr.write_line "weird variant:", $variant

    var gc = count_gc(gc_bool, s, e)
    var gci = ((GC_BINS - 1) * gc).int
    var gc_stat = gc_stats[gci]
    if gc_stat.n < 1000:
      gc_stat = stats

    var local_stats = MedianStats()
    for i in (s+1)..e:
        local_stats.add(values[i])

    # look left and right up to $size bases.
    # also offset by `off` bases to avoid sketchiness near event bounds.
    var
      flank_stats = MedianStats()
      off = 20
      size = SIZE
    for i in max(0, s - size - off)..(s - off):
      flank_stats.add(values[i])
    for i in (e + off)..min(e + size + off, values.high):
      flank_stats.add(values[i])

    var tmp = @[gc]
    if variant.info.set("GCF", tmp) != Status.OK:
        quit "couldn't set GCF"

    var floats = newSeq[float32](variant.vcf.n_samples)

    get_or_empty(variant, "DHFC", floats)
    var fc = (local_stats.median / stats.median).float32
    floats[sample_i] = fc
    if variant.format.set("DHFC", floats) != Status.OK:
        quit "error setting DHFC in VCF"

    get_or_empty(variant, "DHFFC", floats)
    var ffc = (local_stats.median / flank_stats.median(skip_zeros=true)).float32

    floats[sample_i] = ffc
    if variant.format.set("DHFFC", floats) != Status.OK:
        quit "error setting DHFFC in VCF"

    get_or_empty(variant, "DHBFC", floats)
    var gfc = (local_stats.median / gc_stat.median).float32
    result = gfc
    floats[sample_i] = gfc
    if variant.format.set("DHBFC", floats) != Status.OK:
        quit "error setting DHBFC in VCF"

    if discordants.len == 0: return

    if isBnd or variant.ALT[0].startswith("<DEL") or (variant.ALT[0] != "<" and len(variant.REF) > len(variant.ALT[0])):
      var ints = newSeq[int32](variant.vcf.n_samples)
      get_or_empty(variant, "DHSP", ints)
      ints[sample_i] = discordants.count(s.int, e.int, i99).int32
      if variant.format.set("DHSP", ints) != Status.OK:
          quit "error setting DHSP in VCF"

proc fill_stats*[T](depths: var seq[T], stats:var MedianStats, gc_stats:var seq[MedianStats], gc_count:var seq[float32], step:int, target_length:int) =
  stats.clear()
  for v in depths:
    if v == 0: continue
    stats.add(v)

  # for each window of length step, gc_count holds the proportion of bases that were G or C
  # now, for each window, we determine the gc bin (multiply by 20 to get the i) and update the
  # stats for that bin.
  var wi = -1
  for w0 in countup(0, target_length - step, step):
      wi += 1
      if gc_count[wi] < 0: continue
      var gci = ((GC_BINS - 1) * gc_count[wi]).int
      # get the correct stat for the gc in this window and update it.
      for i in w0..<(w0 + step):
          if depths[i] == 0: continue
          gc_stats[gci].add(depths[i])

proc get_isize_distribution*(bam:Bam, n:int=5_000_000, skip=1_500_000): MedianStats =
    result = MedianStats()
    var k = 0
    var mates = false
    for aln in bam:
        if aln.mate_pos != -1:
            mates = true

        k += 1
        if k > 10000 and not mates:
            info("no mates found. assuming single end reads")
            return
        if k < skip: continue

        if aln.mapping_quality == 0: continue
        #if not aln.flag.proper_pair: continue
        if aln.start > aln.mate_pos or aln.isize <= 0: continue
        result.add(max(0, aln.isize.int))
        if result.n == n:
            return

iterator duphold*(bam:Bam, vcf:VCF, fai:Fai, sample_i:int, step:int=STEP): Variant =
  var
      targets = bam.hdr.targets
      target: Target

  var
    last_chrom = ""
    stats:MedianStats
    gc_stats:seq[MedianStats]
    gc_count:seq[float32]
    w = 21

  var id = bam.get_isize_distribution()
  # we add a lot of stuff to the discs so that we don't miss anything
  # but we use a more stringent i99 to decide what's normal when we check
  # for read-pairs that a deletion--this actually ends up making more pairs
  # support the deletion.
  var i99 = id.percentile(99)
  var discs = newSeqOfCap[Discordant](16384)

  proc ifun(aln:Record, posns: var seq[mrange]) =
      if aln.mapping_quality == 0:
          return
      var f = aln.flag
      if f.unmapped or f.secondary or f.qcfail or f.dup: return
      var stop = aln.stop
      if stop - aln.start < 500:
        posns.add((aln.start.int, stop.int, 1))
      else: # for long reads, we actually parse the cigar.
        var pos = aln.start.int
        for op in aln.cigar:
          var c = op.consumes
          if not c.reference: continue

          if c.query:
            if len(posns) == 0 or pos != posns[len(posns)-1].stop:
              # for this function, we want the span of the event so we increment start .. stop.$
              posns.add((pos, pos+op.len, 1))
            else:
              posns[len(posns)-1].stop = pos + op.len
          pos += op.len

      if aln.isize < i99: return
      if stop >= aln.mate_pos: return
      if aln.tid != aln.mate_tid: return
      #var cig = aln.cigar
      #if cig[cig.len-1].op != CigarOp(soft_clip):
      discs.add(Discordant(left:stop.uint32, right: aln.mate_pos.uint32))
      #else:
      #  discs.add(Discordant(left:(aln.stop.uint32 + cig[cig.len-1].len.uint32), right: aln.mate_pos.uint32))

  var depths : Fun[int16] = Fun[int16](f:ifun)
  var gc_bool: seq[bool]

  for variant in vcf:
      if variant.CHROM == last_chrom:
          discard variant.duphold(depths.values, sample_i, stats, gc_stats, gc_bool, w=w, discordants=discs, i99=i99)

          yield variant
          continue

      last_chrom = $variant.CHROM
      stats = MedianStats()
      gc_stats = newSeq[MedianStats](GC_BINS)
      for i, g in gc_stats:
          gc_stats[i] = MedianStats()
      gc_count = gc_count[0..<0]
      var
        i = targets.find(last_chrom)

      if i == -1:
          yield variant
          continue

      target = targets[i]
      depths.values.setLen(target.length.int + 1)
      zeroMem(depths.values[0].addr, depths.values.len * depths.values[0].sizeof)
      gc_bool.setLen(target.length.int)
      zeroMem(gc_bool[0].addr, gc_bool.len * gc_bool[0].sizeof)
      discs.setLen(0)
      GC_fullCollect()

      # required for recent nim version
      var odepths:seq[Fun[int16]] = @[depths]

      discard genoiser[int16](bam, odepths, $target.name, 0, target.length.int)

      gc_count = fai.gc_content(last_chrom, step, gc_bool)
      depths.values.fill_stats(stats, gc_stats, gc_count, step, target.length.int)

      sort(discs, proc(a, b:Discordant): int =
          if a.left == b.left:
              return a.right.int - b.right.int
          return a.left.int - b.left.int
      )

      discard variant.duphold(depths.values, sample_i, stats, gc_stats, gc_bool, w=w, discordants=discs, i99=i99)
      yield variant

type inner = tuple[left: int, right:int]

proc innerCI(v:Variant): inner {.inline.} =
    result = (v.start.int, v.stop.int)
    var ci = newSeq[int32](2)
    if v.info.get("CIPOS95", ci) == Status.OK:
        result[0] += ci[1]
        if v.info.get("CIEND95", ci) == Status.OK:
            result[1] += ci[0]

    if result[1] - result[0] < 5000: return
    # for a longer variant, get a tighter bound
    if v.info.get("CIPOS", ci) == Status.OK:
        result[0] += ci[1]
    if v.info.get("CIEND", ci) == Status.OK:
        result[1] += ci[0]
    if result[1] < result[0]:
        result = (result[1], result[0])

type snpset = ref object
  starts: seq[int32]
  nalts: seq[int8] # 0 == HOM-REF, 1 == HET, 2 == HOM-ALT, 3 == UNKNOWN, 4 == Low-Q

template AB(sample_i:int, AD:seq[int32]): float32 =
  var a = AD[2*sample_i+1].float32
  a / max(1, a + AD[2*sample_i].float32)

template DP(sample_i: int, AD:seq[int32]): int32 =
  AD[2*sample_i+1] + AD[2*sample_i]

proc hq(sample_i: int, GQ: seq[int32], AD: seq[int32], alts: seq[int8], min_dp=9): bool =
  if sample_i.DP(AD) < min_dp: return false
  if GQ[sample_i] < 20: return false
  var ab = sample_i.AB(AD)
  case alts[sample_i]:
    of 0:
      return ab < 0.01
    of 1:
      return ab > 0.3 and ab < 0.7
    of 2:
      return ab > 0.99
    else:
      return false

proc read(snps:VCF, chrom:string): snpset =
  ## read all bi-allelic, high-quality snps into the snpset.
  info("starting to read snps for chrom: " & chrom)

  result = snpset(starts:newSeqOfCap[int32](32768), nalts:newSeqOfCap[int8](32768))
  var ad = newSeq[int32](2)
  var gt = newSeq[int32](1)
  var gq = newSeq[int32](1)
  for snv in snps.query(chrom):
    if len(snv.REF) > 1 or len(snv.ALT) != 1 or len(snv.ALT[0]) != 1 or snv.QUAL < 20: continue
    if snv.FILTER notin ["", ".", "PASS"]: continue
    var alts = snv.format.genotypes(gt).alts

    result.starts.add(snv.start.int32)

    var sample_alts = alts[0]
    if sample_alts == -1:
      result.nalts.add(3)
      continue

    if snv.format.get("AD", ad) != Status.OK:
      quit "expected AD field in snps VCF"
    var st = snv.format.get("GQ", gq)
    if st == Status.UnexpectedType:
      var o = cast[seq[float32]](gq)
      if snv.format.get("GQ", o) != Status.OK:
        quit "unable to extract GQ field in snps VCF"
      gq.setLen(o.len)
      for i, v in o:
        gq[i] = v.int32
    elif st != Status.OK:
      quit "expected GQ field in snps VCF"

    if not hq(0, gq, ad, alts):
      result.nalts.add(4)
    else:
      result.nalts.add(sample_alts)

  info("done reading " & $result.starts.len & " bi-allelic snps for chrom: " & chrom)


proc annotate*(snps:snpset, variant:Variant, sample_i:int) =
  ## annotate a structural variant with the snps.
  if len(snps.starts) == 0: return

  # dont annotate BNDs (or INVs)
  var alt = variant.ALT[0]
  if (':' in alt and ('[' in alt or ']' in alt)) or alt == "<INV>": return

  var dhgt = newSeq[int32](5 * variant.vcf.n_samples)
  get_or_empty(variant, "DHGT", dhgt, 5)
  dhgt[5*sample_i] = 0
  dhgt[5*sample_i+1] = 0
  dhgt[5*sample_i+2] = 0
  dhgt[5*sample_i+3] = 0
  dhgt[5*sample_i+4] = 0

  var ci = innerCI(variant)
  var i = lowerBound(snps.starts, ci.left.int32, system.cmp)
  var n = 0
  while i < snps.starts.len and snps.starts[i] < ci.right:
    dhgt[5 * sample_i + snps.nalts[i]] += 1
    n += 1
    i += 1

  if n == 0:
    return

  if variant.format.set("DHGT", dhgt) != Status.OK:
      quit "error setting DHGT in VCF"

proc sample_name(b:Bam): string =
  for line in ($b.hdr).split("\n"):
    if not line.startswith("@RG"): continue
    var tmp = line.split("SM:")[1].split("\t")[0].strip()
    if result != "" and tmp != result:
      raise newException(ValueError, "found multiple samples in bam:" & result & "," & tmp)
    result = tmp

proc main(argv: seq[string]) =
  var L = newConsoleLogger(fmtStr=verboseFmtStr)
  addHandler(L)

  let doc = format("""
  version: $version

  Usage: duphold [options]

Options:
  -v --vcf <path>           path to sorted SV VCF/BCF
  -b --bam <path>           path to indexed BAM/CRAM
  -f --fasta <path>         indexed fasta reference.
  -s --snp <path>           optional path to snp/indel VCF/BCF with which to annotate SVs. BCF is highly recommended as it's much faster to parse.
  -t --threads <int>        number of decompression threads. [default: 4]
  -o --output <string>      output VCF/BCF (default is VCF to stdout) [default: -]
  -d --drop                 drop all samples from a multi-sample --vcf *except* the sample in --bam. useful for parallelization by sample followed by merge.
  -h --help                 show help
  """, @["version", dupholdVersion])

  let args = docopt(doc, argv=argv)
  var
    fai:Fai
    vcf:VCF
    bam:Bam
    ovcf:VCF
    sample_i: int
    snps:VCF

  if $args["--fasta"] == "nil":
    echo doc
    quit "--fasta is required"
  if $args["--vcf"] == "nil":
    echo doc
    quit "--vcf is required"
  if $args["--bam"] == "nil":
    echo doc
    quit "--bam is required"

  if not open(fai, $args["--fasta"]):
    quit "invalid --fasta: " & $args["--fasta"]
  if not open(vcf, $args["--vcf"]):
    quit "invalid --vcf: " & $args["--vcf"]

  if not open(ovcf, $args["--output"], mode="w"):
    quit "unable to open output vcf"

  var threads = parseInt($args["--threads"])
  setenvvars()
  if $args["--snp"] != "nil" and not open(snps, $args["--snp"], threads=threads):
      quit "couldn't open snp VCF at:" & $args["--snp"]


  if vcf.header.add_info("GCF", "1", "Float", "GC-content fraction for the variant region between 0 and 1.") != Status.OK:
      quit "unable to add to header"
  #if vcf.header.add_format("DHZ", "1", "Float", "duphold z-score for depth") != Status.OK:
  #    quit "unable to add to header"
  if vcf.header.add_format("DHFC", "1", "Float", "duphold depth fold-change") != Status.OK:
      quit "unable to add to header"
  if vcf.header.add_format("DHBFC", "1", "Float", "duphold depth fold-change compared to bins with matching GC") != Status.OK:
      quit "unable to add to DHBFC header"
  if vcf.header.add_format("DHFFC", "1", "Float", &"duphold depth flank fold-change compared to {SIZE}bp left and right of event") != Status.OK:
      quit "unable to add to DHFFC header"
  if vcf.header.add_format("DHSP", "1", "Integer", "duphold count of spanning read-pairs") != Status.OK:
      quit "unable to add DHSP to header"
  if snps != nil and vcf.header.add_format("DHGT", "5", "Integer", "count of hi-quality hom-ref, het, hom-alt, unknown, low-quality SNP variant genotypes in the event") != Status.OK:
      quit "unable DHGT to add to header"

  open(bam, $args["--bam"], index=true, threads=threads, fai=($args["--fasta"]))
  if bam == nil:
      quit "could not open bam file"
  if bam.idx == nil:
      quit "could not open bam index"
  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 510)
  sample_i = vcf.samples.find(bam.sample_name)
  if sample_i == -1:
    if getEnv("DUPHOLD_SAMPLE_NAME") != "":
      stderr.write_line "[duphold] " & "trying sample name from env:" & getEnv("DUPHOLD_SAMPLE_NAME")
      sample_i = vcf.samples.find(getEnv("DUPHOLD_SAMPLE_NAME"))
    if sample_i == -1:
      quit "couldn't find sample from bam:" & bam.sample_name & " or ENV in vcf which had:" & join(vcf.samples, ",")

  if args["--drop"]:
    vcf.set_samples(@[bam.sample_name])
    sample_i = 0
  ovcf.header = vcf.header

  if snps != nil:
    if snps.samples.find(bam.sample_name) == -1:
      quit "[duphold] couldn't find sample:" & bam.sample_name & " in snps vcf which had:" & join(snps.samples, ",")
    snps.set_samples(@[bam.sample_name])

  if not ovcf.write_header():
      quit "couldn't write vcf header"

  var last_chrom = "".cstring
  var snpst: snpset
  for variant in bam.duphold(vcf, fai, sample_i):
      if snps != nil:
          if last_chrom != variant.CHROM:
            snpst = snps.read($variant.CHROM)
            last_chrom = variant.CHROM

          snpst.annotate(variant, sample_i)
      if not ovcf.write_variant(variant):
          quit "couldn't write variant"

  vcf.close()
  bam.close()
  if snps != nil:
    snps.close()
  ovcf.close()
  stderr.write_line "[duphold] finished"

when isMainModule:
    main(commandLineParams())

