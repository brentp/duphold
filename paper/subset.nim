import random
import strutils
import algorithm
import hts

var path = "ALL.wgs.merged-svs.vcf.gz"

var vcf:VCF

var variants = newSeq[Variant]()

if not open(vcf, path, threads=3):
 quit "couldn't open vcf"

for v in vcf:
  if v.ALT[0].startswith("<INS"): continue
  variants.add(v.copy())
echo variants.len

discard vcf.header.add_info("CIPOS", "2", "Integer", "")
discard vcf.header.add_info("CIEND", "2", "Integer", "")

var wtr:VCF
var svt = ""
var ci = @[2'i32, 2]
for n in @[100, 1000, 5000, 10000, 20000, 35000, 50000]:
  shuffle(variants)

  var subset = variants[0..<n]

  sort(subset, proc(a, b:Variant): int = 
      if a.CHROM == b.CHROM:
        return a.start - b.start
      return cmp(a.CHROM, b.CHROM)
  )

  if not open(wtr, path & "." & $n & ".vcf.gz", mode="w"):
    quit "couldnt open for writing"
  wtr.copy_header(vcf.header)
  if not wtr.write_header():
      quit "coudln't write header"

  for ss in subset:
    var s = ss.copy()
    discard s.info.set("CIPOS", ci)
    discard s.info.set("CIEND", ci)
    discard s.info.get("SVTYPE", svt)
    if svt != "DUP" and svt != "DEL" and svt != "INS" and svt != "INV":
        var ch: string
        if svt.startswith("DEL"):
            ch = "DEL"
            discard s.info.set("SVTYPE", ch)
        elif svt == "ALU" or svt == "LINE1" or svt == "SVA" or svt == "SVA":
            ch = "INS"
            discard s.info.set("SVTYPE", ch)
        elif svt.startswith("CNV"):
            ch = "DEL"
            discard s.info.set("SVTYPE", ch)


    if not wtr.write_variant(s):
      quit "couldn't write variant"

  wtr.close()

