import hts
import strutils
import os
import math
import docopt
import genoiser
import ./dupholdpkg/version

const STEP = 200


type Stats* = ref object
    n*: int
    S*: float64
    m*: float64
    t*: int # total used for DHD

proc update*(s:var Stats, d:int, include_zero:bool) {.inline.} =
    ## streaming mean, sd
    ## https://dsp.stackexchange.com/questions/811/determining-the-mean-and-standard-deviation-in-real-time
    if (not include_zero) and d == 0:
        return
    s.n += 1
    var mprev = s.m
    var df = d.float64
    s.m += (df - s.m) / s.n.float64
    s.S += (df - s.m) * (df - mprev)

proc clear*(s:var Stats) =
    s.n = 0
    s.S = 0
    s.m = 0
    s.t = 0

proc stddev*(s:Stats): float64 {.inline.} =
   sqrt(s.S/s.n.float64)

proc addm*(s:Stats, d:int) {.inline.} =
  # update only the t and n.
  # not that this method is never used on a struct where the update method is also used.
  s.n += 1
  s.t += d

proc dropm*(s:Stats, d:int) {.inline.} =
  # drp the t and the n.
  s.n -= 1
  s.t -= d

proc mean*(s:Stats):float64 {.inline.} =
    return s.t.float64 / s.n.float64

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

proc count_gc(s:string): float32 =
    for c in s:
        if c == 'C' or c == 'G':
            result += 1
    result /= s.len.float32

proc gc_content*(fai:Fai, chrom:string, step:int): seq[float32] =
    var s = fai.get(chrom).toUpperAscii()
    result = newSeq[float32]((s.len/step+1).int)

    for i, c in s:
        if c == 'C' or c == 'G':
            result[(i/step).int] += 1
        elif c == 'N':
            result[(i/step).int] -= 1
    for i, v in result:
        result[i] = v / step.float32
    return result

proc get_or_empty[T](variant:Variant, field:string, input:var seq[T]) =
  ## if, for example we've already annotated a sample in the VCF with duphold
  ## we dont want to overwite those values with nan so try to grab existing
  ## values but otherwise make an empty array.
  if variant.format.get(field, input) == Status.OK:
      return

  if input.len != variant.vcf.n_samples:
    input.set_len(variant.vcf.n_samples)
  for i, f in input:
    when T is float32:
        input[i] = cast[T](bcf_float_missing)
    elif T is int32:
        input[i] = int32.low
    else:
        quit "unknown type in get_or_empty"

proc check_rapid_depth_change*[T](start:int, stop:int, values: var seq[T], w:int=7): int32 =
    ## if start and end indicate the bounds of a deletion, we can often expect to see a rapid change in
    ## depth at or near the break-point.
    var
      # could use CI for this, but larger CI == less confident anyway.
      dist = min(80, max(25, 0.05 * (stop - start).float64).int)
      # if we see too many changes then we can't trust the result.
      changes = 0
      d: float64
      last_change = 0

    for bi, bp in @[start, stop]:
      var
        close_changes = 0
        cs = max(w, bp - dist - w)
        left = Stats()
        right = Stats()

      for i in (cs-w)..<(cs):
          left.addm(values[i])
      for i in cs..<(cs + w):
          right.addm(values[i])

      for k in cs..min(bp + dist + w, values.len - w - 1):
        if (k - last_change) > w:
            var
                lm = left.mean
                rm = right.mean
            if lm < 8 and rm < 8:
              left.dropm(values[k-w])
              left.addm(values[k])

              right.dropm(values[k])
              right.addm(values[k + w])
              continue

            if bi == 0:
              d = rm / lm
            else:
              d = lm / rm
            #if (k > start - 10 and k < start + 10) or (k > stop - 10 and k < stop + 10):
            #    echo k, " ",d

            # in addition to normal change (1.25), we have special-case here for when we're very close to either break-point.
            # the changes <= bi makes sure we dont mess up an existing, better change that was
            # already found.
            var close = (((k - start).abs <= 3) or (k - stop).abs <= 3)
            if close_changes == bi and d > 1.25 or (changes <= bi and d > 1.17 and close):
              # if we are right near the break-point and we already found a conflicting change
              # but this one is very good, we over-ride
              if close and close_changes == bi and result < 0 and d > 1.25:
                  changes = bi
                  close_changes = bi + 1
                  result = 0
              changes += 1
              #echo "DUP position:", k, " fc:", d, " changes:", changes, " left:", left.mean, " right:", right.mean , "close:", close
              last_change = k
              result += 1
            if d < 0.7 or (changes <= bi and d < 0.8 and close):
              changes += 1
              last_change = k
              #echo "DEL position:", k, " fc:", d, " changes:", changes, " left:", left.mean, " right:", right.mean
              result -= 1
        # now update left and right
        left.dropm(values[k-w])
        left.addm(values[k])

        right.dropm(values[k])
        right.addm(values[k + w])

    #echo "FINISHED:", "changes:", changes, "result:", result
    if changes > 2:
        result = 0

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

proc get_bnd_mate_pos(variant:Variant): int {.inline.} =
    return get_bnd_mate_pos(variant.ALT[0], $variant.CHROM)

proc duphold*[T](variant:Variant, values:var seq[T], sample_i: int, stats:var Stats, gc_stats:var seq[Stats], fai:Fai): float64 =
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

    var ss = fai.get($variant.CHROM, s, e).toUpperAscii()
    var gc = count_gc(ss)
    var gci = (19 * gc).int
    var gc_stat = gc_stats[gci]
    #echo gci, " ", gc, " ", gc_stat.n, " ", gc_stat.S

    var local_stats = Stats()
    for i in (s+1)..e:
        local_stats.update(values[i], true)

    var tmp = @[gc]
    if variant.info.set("GCF", tmp) != Status.OK:
        quit "couldn't set GCF"

    var floats = newSeq[float32](variant.vcf.n_samples)

    get_or_empty(variant, "DHFC", floats)
    var fc = local_stats.m / stats.m
    floats[sample_i] = fc.float32
    if variant.format.set("DHFC", floats) != Status.OK:
        quit "error setting DHFC in VCF"

    get_or_empty(variant, "DHBZ", floats)
    var gcz = (local_stats.m - gc_stat.m) / sqrt(gc_stat.S/gc_stat.n.float64)
    floats[sample_i] = gcz.float32
    if variant.format.set("DHBZ", floats) != Status.OK:
        quit "error setting DHBZ in VCF"

    get_or_empty(variant, "DHBFC", floats)
    var gfc = local_stats.m / gc_stat.m
    result = gfc
    floats[sample_i] = gfc.float32
    if variant.format.set("DHBFC", floats) != Status.OK:
        quit "error setting DHBFC in VCF"

    var ints = newSeq[int32](variant.vcf.n_samples)
    get_or_empty(variant, "DHD", ints)
    ints[sample_i] = check_rapid_depth_change(s, e, values)
    if variant.format.set("DHD", ints) != Status.OK:
        quit "error setting DHD in VCF"

proc fill_stats*[T](depths: var seq[T], stats:var Stats, gc_stats:var seq[Stats], gc_count:var seq[float32], step:int, target_length:int) =
  stats.clear()
  for v in depths:
    stats.update(v, false)

  var dmax = int(stats.m + 4 * sqrt(stats.m))
  stats.clear()
  for v in depths:
    if v < dmax:
      stats.update(v, false)

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
          gc_stats[gci].update(depths[i], false)


iterator duphold*(bam:Bam, vcf:VCF, fai:Fai, sample_i:int, step:int=STEP): Variant =
  var depths : Fun[int16] = Fun[int16](values: newSeq[int16](), f:idepthfun)
  var
      targets = bam.hdr.targets
      target: Target

  var
    last_chrom = ""
    stats:Stats
    gc_stats:seq[Stats]
    gc_count:seq[float32]


  for variant in vcf:
      if variant.CHROM == last_chrom:
          discard variant.duphold(depths.values, sample_i, stats, gc_stats, fai)

          yield variant
          continue

      target = nil
      last_chrom = $variant.CHROM
      stats = Stats()
      gc_stats = newSeq[Stats](20)
      for i, g in gc_stats:
          gc_stats[i] = Stats()
      gc_count = gc_count[0..<0]
      var
        start = 0
        i = targets.find(last_chrom)

      depths.values.set_len(targets[i].length.int + 1)
      zeroMem(depths.values[0].addr, depths.values.len * sizeof(depths.values[0]))

      if i == -1:
          yield variant
          continue

      target = targets[i]
      discard genoiser[int16](bam, @[depths], target.name, 0, target.length.int)

      gc_count = fai.gc_content(last_chrom, step)
      depths.values.fill_stats(stats, gc_stats, gc_count, step, target.length.int)

      discard variant.duphold(depths.values, sample_i, stats, gc_stats, fai)
      yield variant

proc sample_name(b:Bam): string =
    for line in ($b.hdr).split("\n"):
        if not line.startswith("@RG"): continue
        var tmp = line.split("SM:")[1].split("\t")[0].strip()
        if result != "" and tmp != result:
            raise newException(ValueError, "found multiple samples in bam:" & result & "," & tmp)
        result = tmp


proc main(argv: seq[string]) =

  let doc = format("""
  version: $version

  Usage: duphold [options]

Options:
  -v --vcf <path>           path to sorted VCF/BCF
  -b --bam <path>           path to indexed BAM/CRAM
  -f --fasta <path>         indexed fasta reference.
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

  if $args["--fasta"] == "nil":
    quit "--fasta is required"
  if $args["--vcf"] == "nil":
    quit "--vcf is required"
  if $args["--bam"] == "nil":
    quit "--bam is required"

  if not open(fai, $args["--fasta"]):
    quit "invalid --fasta: " & $args["--fasta"]
  if not open(vcf, $args["--vcf"]):
    quit "invalid --vcf: " & $args["--vcf"]

  if not open(ovcf, $args["--output"], mode="w"):
    quit "unable to open output vcf"

  if vcf.header.add_info("GCF", "1", "Float", "GC-content fraction for the variant region betwee 0 and 1.") != Status.OK:
      quit "unable to add to header"
  #if vcf.header.add_format("DHZ", "1", "Float", "duphold z-score for depth") != Status.OK:
  #    quit "unable to add to header"
  if vcf.header.add_format("DHFC", "1", "Float", "duphold depth fold-change") != Status.OK:
      quit "unable to add to header"
  if vcf.header.add_format("DHBZ", "1", "Float", "duphold z-score for depth compared to bins with matching GC") != Status.OK:
      quit "unable to add to header"
  if vcf.header.add_format("DHBFC", "1", "Float", "duphold depth fold-change compared to bins with matching GC") != Status.OK:
      quit "unable to add to header"
  if vcf.header.add_format("DHD", "1", "Integer", "duphold rapid change in depth at one of the break-points (1 for higher. 0 for no or conflicting changes. -1 for drop, 2 for both break points)") != Status.OK:
      quit "unable to add to header"


  open(bam, $args["--bam"], index=true, threads=parseInt($args["--threads"]), fai=($args["--fasta"]))
  if bam == nil:
      quit "could not open bam file"
  if bam.idx == nil:
      quit "could not open bam index"
  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 511)
  sample_i = vcf.samples.find(bam.sample_name)
  if sample_i == -1:
      quit "couldn't find sample:" & bam.sample_name & " in vcf which had:" & join(vcf.samples, ",")
  if args["--drop"]:
    vcf.set_samples(@[bam.sample_name])
    sample_i = 0
  ovcf.header = vcf.header

  if not ovcf.write_header():
      quit "couldn't write vcf header"

  for variant in bam.duphold(vcf, fai, sample_i):
      if not ovcf.write_variant(variant):
          quit "couldn't write variant"

  ovcf.close()
  vcf.close()
  bam.close()

when isMainModule:
    main(commandLineParams())
