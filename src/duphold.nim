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
import ./dupholdpkg/bin_table

const STEP = 200

type MedianStats* = object
  ## since we expect nearly all values to be < 1000
  ## we can store counts of observed values in a fixed-size array
  ## and get median (or any percentile) from those.
  ## we can even keep a moving window of medians as we can drop
  ## values by decrementing the counts.
  counts*: array[1000, int]
  i_median: int
  i_ok: bool
  n*: int

proc addm*(m:var MedianStats, d:int, include_zero:bool) {.inline.} =
    if not include_zero and d == 0:
        return
    if d > m.counts.high:
      m.counts[m.counts.high] += 1
    else:
      m.counts[d] += 1
    m.n += 1
    m.i_ok = false

proc dropm*(m:var MedianStats, d:int, include_zero:bool) {.inline.} =
    if d == 0 and (not include_zero):
        return
    m.i_ok = false
    if d > m.counts.high:
      m.counts[m.counts.high] -= 1
    else:
      m.counts[d] -= 1
    m.n -= 1

proc clear*(m:var MedianStats) =
    zeroMem(m.counts[0].addr, sizeof(m.counts[0]) * m.counts.len)
    m.i_ok = false
    m.i_median = 0
    m.n = 0

proc median*(m:var MedianStats): int {.inline.} =
    if m.i_ok:
      return m.i_median

    var cum = 0
    var stop_n = (0.5 + m.n.float32 * 0.5'f32).int
    for i, cnt in m.counts:
        cum += cnt
        if cum >= stop_n:
            m.i_median = i
            m.i_ok = true
            break
    return m.i_median

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
    return &"MedianStats(n:{m.n}, median:{m.median} vals: {m.counts[0..100]})"

type Stats* = ref object
    n*: int
    mean*: float64

proc addm*(s:var Stats, d:int, include_zero:bool) {.inline.} =
    ## streaming mean
    if (not include_zero) and d == 0:
        return
    s.n += 1
    s.mean += (d.float64 - s.mean) / s.n.float64

proc dropm*(s:var Stats, d:int, include_zero:bool) {.inline.} =
    ## streaming mean, sd
    ## https://lingpipe-blog.com/2009/07/07/welford-s-algorithm-delete-online-mean-variance-deviation/
    if (not include_zero) and d == 0:
        return
    var df = d.float64
    if s.n == 1:
        s.mean = 0
        s.n = 0
        return
    var mprev = (s.n.float64 * s.mean - df) / (s.n - 1).float64
    s.n -= 1
    #s.S -= (df - s.mean) * (df - mprev)
    s.mean = mprev
    if s.mean < 0:
        s.mean = 0

proc clear*(s:var Stats) =
    s.n = 0
    s.mean = 0

proc fc*(a:Stats, b:Stats): float64 {.inline.} =
    ## fold-change of a relative to b so that value is always > 1.
    if a.mean > b.mean:
      return a.mean / b.mean
    return -b.mean / a.mean

type Discordant* = ref object
  left*: uint32
  right*: uint32

proc `$`*(d:Discordant): string =
    return &"Discordant(left: {d.left}, right: {d.right}"

proc idepthfun*(aln:Record, posns:var seq[mrange]) =
  ## depthfun is an example of a `fun` that can be sent to `genoiser`.
  ## it sets reports the depth at each position
  if aln.mapping_quality == 0: return
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  #refposns(aln, posns)
  posns.add((aln.start, aln.stop, 1))

proc find(targets:seq[Target], chrom:string): int =
    for i, t in targets:
        if t.name == chrom: return i
    return -1

{.push checks: off.}
proc count_gc(gc: var seq[bool], s:int, e:int): float32 {.inline.} =
    var tot = 0
    for i in s..<e:
        tot += gc[i].int
    return tot.float32 / max(1, e - s).float32

proc gc_content*(fai:Fai, chrom:string, step:int, gc_bool: var seq[bool]): seq[float32] =
    var s = fai.cget(chrom)
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
    free(s)
    return result
{.push checks: on.}

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


type fc_position = tuple[fc:float64, position:int]

proc check_rapid_depth_change_end*[T](position:int, dist:int, values:var seq[T], w:int, cutoff:float64): fc_position =
    var
      cs = max(w, position - dist - w)
      left = Stats()
      right = Stats()

    for i in (cs-w)..<(cs):
      left.addm(values[i], true)
    for i in cs..<(cs + w):
      right.addm(values[i], true)

    result = (0'f64, -1.int)

    var gt1p1 = 0
    var lastgt = -10000

    for k in cs..min(position + dist + w, values.len - w - 1):

        left.dropm(values[k-w], true)
        left.addm(values[k], true)

        right.dropm(values[k], true)
        right.addm(values[k + w], true)

        if values[k] == 0: continue

        if left.mean < 7 and right.mean < 7:
            continue

        #when defined(bdebug):
        #  echo k, " ", left.fc(right), " rm:", right.mean, " lm:", left.mean

        if left.fc(right).abs > cutoff:
            if k - lastgt > w:
                gt1p1 += 1
            lastgt = k

        if left.fc(right).abs > result.fc.abs:
            result.fc = left.fc(right)
            result.position = k
    if gt1p1 >= 3 and result.fc < cutoff:
        result.position = -1
        result.fc = 0

proc toIndex(fc: float64, vhigh:int, mult:float64): int {.inline.} =
    if fc.classify == fcInf:
        return vhigh
    result = (fc * mult).int
    if result > vhigh:
        result = vhigh

proc find_fc_cutoffs[T](values: var seq[T], fdrs:seq[float64], w:int): seq[float64] =
    var mult = 3000'f64

    var
      cs = w
      left = Stats()
      right = Stats()

    var t = cpuTime()

    for i in (cs-w)..<(cs):
      left.addm(values[i], true)
    for i in cs..<(cs + w):
      right.addm(values[i], true)

    result = newSeq[float64](fdrs.len)

    var fcs = newSeq[int](10000)
    var added = 0

    for k in cs..(values.len - w - 1):

        left.dropm(values[k-w], true)
        left.addm(values[k], true)

        right.dropm(values[k], true)
        right.addm(values[k + w], true)

        if right.n != w:
            quit "bad"
        if left.n != w:
            quit "bad"

        #if k mod 100000 == 0:
        #    echo left.mean, " ", right.mean
        if values[k] == 0: continue

        if left.mean < 7 and right.mean < 7:
            continue

        if left.mean > 800 or right.mean > 800: continue

        var fc = left.fc(right)

        added += 1
        fcs[fc.abs.toIndex(fcs.high, mult)] += 1

    var ncum = 0
    var ni = 0
    for k, fdr in fdrs:
        while ncum < int((1 - fdr) * (added).float64):
            #ncum += negatives[ni]
            ncum += fcs[ni]
            ni += 1
        # ni is the index and now we calculate the fold-change requred
        result[k] = ni.float64 / mult


proc check_rapid_depth_change*[T](start:int, stop:int, values: var seq[T], cutoffs:seq[float64], w:int): int32 =
    ## if start and end indicate the bounds of a deletion, we can often expect to see a rapid change in
    ## depth at or near the break-point.
    ## cutoffs should have length 2 and the first value is less stringnet-- for looking near the
    ## break-point. the 2nd value is more stringent in case we have to look farther from the break.
    var
      # if we see too many changes then we can't trust the result.
      changes = 0
      d: float64
      last_change = 0
      dist = min(80, max(25, 0.08 * (stop - start).float64).int)

      left = check_rapid_depth_change_end(start, 4, values, w, cutoffs[0])
      right = check_rapid_depth_change_end(stop, 4, values, w, cutoffs[0])

    when defined(bdebug):
      echo "start:", start, " stop:", stop
      echo left
      echo right
      echo result

    if left.fc.abs < cutoffs[0] and left.position != -1:
      left = check_rapid_depth_change_end(start, dist, values, w, cutoffs[0])
    if right.fc.abs < cutoffs[0] and right.position != -1:
      right = check_rapid_depth_change_end(stop, dist, values, w, cutoffs[0])

    when defined(bdebug):
      echo "start:", start, " stop:", stop
      echo left
      echo right
      echo result


    # if position is -1, we had no data for 1 side.
    if (left.fc > 0) == (right.fc > 0) and (left.position != -1 and right.position != -1): return 0
    if left.position == -1 and right.position == -1: return 0

    if left.position != -1:
      if left.fc > cutoffs[1]: result -= 1
      if left.fc < -cutoffs[1]: result += 1

    if right.position != -1:
      if right.fc < -cutoffs[1]: result -= 1
      if right.fc > cutoffs[1]: result += 1

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
    var ilow = lowerBound(discordants, Discordant(left:uint32(s - i99)), proc(a, b: Discordant):int =
        if a.left == b.left:
            return a.right.int - b.right.int
        return a.left.int - b.left.int
    )
    if ilow > ilow: ilow -= 1
    for i in ilow..discordants.high:
        var disc = discordants[i]
        if disc.left.int > s + slop: break
        if s - disc.left.int > i99: continue

        if disc.right.int - (e - slop) < i99:
          # after correcting for the event, the insert should fall in the expected distribution.
          #echo &"s-e:{s}-{e} supported by: {disc} with corrected span: {(disc.right.int - disc.left.int) - (e - s)}"
          var corrected = (disc.right.int - disc.left.int) - (e - s)
          if corrected < i99 and corrected > -2*slop:
            result += 1

proc get_bnd_mate_pos(variant:Variant): int {.inline.} =
    return get_bnd_mate_pos(variant.ALT[0], $variant.CHROM)

proc duphold*[T](variant:Variant, values:var seq[T], sample_i: int, stats:var MedianStats, gc_stats:var seq[MedianStats], gc_bool:var seq[bool], w:int=7, cutoffs:seq[float64], discordants:var seq[Discordant], i99:int): float64 =
    ## sets FORMAT fields for sample i in the variant and returns the DHBFC value
    var
      s = variant.start
      e = variant.stop

    var bnd = get_bnd_mate_pos(variant)
    if bnd != -1:
      if bnd < s and s - bnd < 20000000:
          e = s
          s = bnd
      elif bnd > e and bnd - e < 20000000:
          s = e
          e = bnd - 1

    if bnd == -1:
        # skip distant BND's as this is not informative
        if ':' in variant.ALT[0]: return -1

    var gc = count_gc(gc_bool, s, e)
    var gci = (19 * gc).int
    var gc_stat = gc_stats[gci]

    var local_stats = MedianStats()
    for i in (s+1)..e:
        local_stats.addm(values[i], true)

    var tmp = @[gc]
    if variant.info.set("GCF", tmp) != Status.OK:
        quit "couldn't set GCF"

    var floats = newSeq[float32](variant.vcf.n_samples)

    get_or_empty(variant, "DHFC", floats)
    var fc = (local_stats.median / stats.median).float32

    floats[sample_i] = fc
    if variant.format.set("DHFC", floats) != Status.OK:
        quit "error setting DHFC in VCF"

    get_or_empty(variant, "DHBFC", floats)
    var gfc = (local_stats.median / gc_stat.median).float32
    result = gfc
    floats[sample_i] = gfc
    if variant.format.set("DHBFC", floats) != Status.OK:
        quit "error setting DHBFC in VCF"

    var ints = newSeq[int32](variant.vcf.n_samples)
    get_or_empty(variant, "DHD", ints)
    ints[sample_i] = check_rapid_depth_change(s, e, values, cutoffs=cutoffs, w=w)
    if variant.format.set("DHD", ints) != Status.OK:
        quit "error setting DHD in VCF"

    if variant.ALT[0] == "<DEL>" or (variant.ALT[0] != "<" and len(variant.REF) > len(variant.ALT[0])):
      get_or_empty(variant, "DHSP", ints)
      ints[sample_i] = discordants.count(s, e, i99).int32
      if variant.format.set("DHSP", ints) != Status.OK:
          quit "error setting DHSP in VCF"

proc fill_stats*[T](depths: var seq[T], stats:var MedianStats, gc_stats:var seq[MedianStats], gc_count:var seq[float32], step:int, target_length:int) =
  stats.clear()
  for v in depths:
    stats.addm(v, true)

  # for each window of length step, gc_count holds the proportion of bases that were G or C
  # now, for each window, we determine the gc bin (multiply by 20 to get the i) and update the
  # stats for that bin.
  var wi = -1
  for w0 in countup(0, target_length - step, step):
      wi += 1
      if gc_count[wi] < 0: continue
      var gci = (19 * gc_count[wi]).int
      # get the correct stat for the gc in this window and update it.
      for i in w0..<(w0 + step):
          gc_stats[gci].addm(depths[i], true)

proc get_isize_distribution(bam:Bam, n:int=4000000, skip=500000): MedianStats =
    result = MedianStats()
    var k = 0
    for aln in bam:
        k += 1
        if k < skip: continue

        if aln.mapping_quality == 0: continue
        #if not aln.flag.proper_pair: continue
        if aln.stop >= aln.mate_pos: continue
        result.addm(aln.isize.abs, true)
        if result.n == n:
            return

iterator duphold*(bam:Bam, vcf:VCF, fai:Fai, sample_i:int, step:int=STEP): Variant =
  var depths : Fun[int16] = Fun[int16](f:idepthfun)
  var
      targets = bam.hdr.targets
      target: Target

  var
    last_chrom = ""
    stats:MedianStats
    gc_stats:seq[MedianStats]
    gc_count:seq[float32]
    dhd_cutoffs: seq[float64]
    w = 1

  var id = bam.get_isize_distribution()
  # we add a lot of stuff to the discs so that we don't miss anything
  # but we use a more stringent i99 to decide what's normal when we check
  # for read-pairs that a deletion--this actually ends up making more pairs
  # support the deletion.
  var i95 = id.percentile(95)
  var i99 = id.percentile(99)
  var discs = newSeqOfCap[Discordant](16384)
  proc ifun(aln:Record, posns: var seq[mrange]) =
      if aln.mapping_quality == 0:
          return
      var f = aln.flag
      if f.unmapped or f.secondary or f.qcfail or f.dup: return
      if aln.stop >= aln.mate_pos: return
      if aln.isize < i95: return
      if aln.isize > 50000000: return
      #var cig = aln.cigar
      #if cig[cig.len-1].op != CigarOp(soft_clip):
      discs.add(Discordant(left:aln.stop.uint32, right: aln.mate_pos.uint32))
      #else:
      #  discs.add(Discordant(left:(aln.stop.uint32 + cig[cig.len-1].len.uint32), right: aln.mate_pos.uint32))

  var iempty : Fun[int16] = Fun[int16](values:newSeq[int16](), f:ifun)
  var gc_bool: seq[bool]

  for variant in vcf:
      if variant.CHROM == last_chrom:
          discard variant.duphold(depths.values, sample_i, stats, gc_stats, gc_bool, cutoffs=dhd_cutoffs, w=w, discordants=discs, i99=i99)

          yield variant
          continue

      target = nil
      last_chrom = $variant.CHROM
      stats = MedianStats()
      gc_stats = newSeq[MedianStats](20)
      for i, g in gc_stats:
          gc_stats[i] = MedianStats()
      gc_count = gc_count[0..<0]
      var
        start = 0
        i = targets.find(last_chrom)

      if i == -1:
          yield variant
          continue

      target = targets[i]
      depths.values = newSeq[int16](target.length.int + 1)
      gc_bool = newSeq[bool](target.length.int)

      info("starting read of bam values for chrom: " & target.name)
      discard genoiser[int16](bam, @[depths, iempty], target.name, 0, target.length.int)
      info("finished reading:" & target.name & ". starting gc-content, discordant sorting, and gc cutoffs")

      gc_count = fai.gc_content(last_chrom, step, gc_bool)
      depths.values.fill_stats(stats, gc_stats, gc_count, step, target.length.int)

      dhd_cutoffs = depths.values.find_fc_cutoffs(fdrs=(@[0.002, 0.0002]), w=w)
      sort(discs, proc(a, b:Discordant): int =
          if a.left == b.left:
              return a.right.int - b.right.int
          return a.left.int - b.left.int
      )
      info("finished gc-content, sorting, and cutoffs")

      discard variant.duphold(depths.values, sample_i, stats, gc_stats, gc_bool, cutoffs=dhd_cutoffs, w=w, discordants=discs, i99=i99)
      yield variant

type inner = tuple[left: int, right:int]

proc innerCI(v:Variant): inner {.inline.} =
    result = (v.start, v.stop)
    var ci = newSeq[int32](2)
    if v.info.get("CIPOS95", ci) == Status.OK:
        result[0] += ci[1]
        if v.info.get("CIEND95", ci) == Status.OK:
            result[0] += ci[0]

    if result[1] - result[0] < 10000: return
    # for a longer variant, get a tighter bound

    if v.info.get("CIPOS", ci) == Status.OK:
        result[0] += ci[1]
    if v.info.get("CIEND", ci) == Status.OK:
        result[1] += ci[0]
    if result[1] < result[0]:
        result = (result[1], result[0])

type snpset = ref object
  starts: seq[int32]
  nalts: seq[int8]
  hetc: seq[int8]

proc read(snps:VCF, chrom:string): snpset =
  ## read all bi-allelic, high-quality snps into the snpset.
  info("starting to read snps for chrom: " & chrom)
  if snps.n_samples != 1:
      quit "only works for 1 snp"
  result = snpset(starts:newSeqOfCap[int32](8192), nalts:newSeqOfCap[int8](8192), hetc:newSeqOfCap[int8](8192))
  var ad = newSeq[int32](2)
  var gt = newSeq[int32](2)
  for snv in snps.query(chrom):
    if len(snv.REF) > 1 or len(snv.ALT) > 1 or len(snv.ALT[0]) > 1 or snv.QUAL < 20: continue
    var gts = snv.format.genotypes(gt)[0]
    var nalts = gts.alts

    result.starts.add(snv.start.int32)
    result.nalts.add(nalts)
    if nalts == 1:
      if snv.format.get("AD", ad) != Status.OK:
        quit "expected AD field in snps VCF"
      result.hetc.add(het_lookup(ad[0], ad[1], 13))
    else:
      result.hetc.add(-1)

  info("done reading " & $result.starts.len & " bi-allelic snps for chrom: " & chrom)


proc annotate*(snps:snpset, variant:Variant, sample_i:int) =
    ## annotate a structural variant with the snps.
    if len(snps.starts) == 0: return

    # dont annotate BNDs (or INVs)
    if ':' in variant.ALT[0] or variant.ALT[0] == "<INV>": return
    var dhet = newSeq[int32](2 * variant.vcf.n_samples)
    get_or_empty(variant, "DHET", dhet, 2)
    dhet[2 * sample_i] = 0
    dhet[2 * sample_i + 1] = 0

    var dhhu = newSeq[int32](3 * variant.vcf.n_samples)
    get_or_empty(variant, "DHHU", dhhu, 3)
    dhhu[3 * sample_i] = 0
    dhhu[3 * sample_i + 1] = 0
    dhhu[3 * sample_i + 2] = 0

    var ci = innerCI(variant)
    var i = lowerBound(snps.starts, ci.left.int32, system.cmp)
    var n = 0
    while snps.starts[i] < ci.right:
        # only compare diploid to triploid allele balance if this was called a HET.
        if snps.nalts[i] == 1:
          if snps.hetc[i] == 0: # diploid
              dhet[2 * sample_i] += 1
          elif snps.hetc[i] == 1: # triploid
              dhet[2 * sample_i + 1] += 1

        elif snps.nalts[i] == 0: # hom-ref
            dhhu[3 * sample_i] += 1
        elif snps.nalts[i] == 2: # hom-alt
            dhhu[3 * sample_i + 1] += 1
        elif snps.nalts[i] == -1: # unknown
            dhhu[3 * sample_i + 2] += 1
        i += 1
        n += 1
    if n == 0:
        return

    if variant.format.set("DHET", dhet) != Status.OK:
        quit "error setting DHET in VCF"
    if variant.format.set("DHHU", dhhu) != Status.OK:
        quit "error setting DHHU in VCF"

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
  -s --snp <path>           optional path to snp/indel VCF with which to annotate SVs.
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
  if $args["--snp"] != "nil" and not open(snps, $args["--snp"], threads=threads):
      quit "couldn't open snp VCF at:" & $args["--snp"]

  if vcf.header.add_info("GCF", "1", "Float", "GC-content fraction for the variant region betwee 0 and 1.") != Status.OK:
      quit "unable to add to header"
  #if vcf.header.add_format("DHZ", "1", "Float", "duphold z-score for depth") != Status.OK:
  #    quit "unable to add to header"
  if vcf.header.add_format("DHFC", "1", "Float", "duphold depth fold-change") != Status.OK:
      quit "unable to add to header"
  if vcf.header.add_format("DHBFC", "1", "Float", "duphold depth fold-change compared to bins with matching GC") != Status.OK:
      quit "unable to add to DHBFC header"
  if vcf.header.add_format("DHD", "1", "Integer", "duphold rapid change in depth at one of the break-points (1 for higher. 0 for no or conflicting changes. -1 for drop, 2 for both break points)") != Status.OK:
      quit "unable to add DHD to header"
  if vcf.header.add_format("DHSP", "1", "Integer", "duphold count of spanning read-pairs") != Status.OK:
      quit "unable to add DHSP to header"

  if snps != nil and vcf.header.add_format("DHET", "2", "Integer", "duphold count of SNP heterozygotes in the event supporting: a normal heterozygote, a triploid heterozygote") != Status.OK:
      quit "unable DHET to add to header"
  if snps != nil and vcf.header.add_format("DHHU", "3", "Integer", "count of hom-ref, hom-alt, unknown SNP variants in the event") != Status.OK:
      quit "unable DHUU to add to header"

  open(bam, $args["--bam"], index=true, threads=threads, fai=($args["--fasta"]))
  if bam == nil:
      quit "could not open bam file"
  if bam.idx == nil:
      quit "could not open bam index"
  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 510)
  sample_i = vcf.samples.find(bam.sample_name)
  if sample_i == -1:
      quit "couldn't find sample:" & bam.sample_name & " in vcf which had:" & join(vcf.samples, ",")
  if args["--drop"]:
    vcf.set_samples(@[bam.sample_name])
    sample_i = 0
  ovcf.header = vcf.header

  if snps != nil:
    if snps.samples.find(bam.sample_name) == -1:
      quit "couldn't find sample:" & bam.sample_name & " in snps vcf which had:" & join(snps.samples, ",")
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

  ovcf.close()
  vcf.close()
  bam.close()

when isMainModule:
    main(commandLineParams())

