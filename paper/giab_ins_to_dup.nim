import hts
import strformat
import os
import osproc
import tables

var
  vcf_path = commandLineParams()[0]
  fasta_path = commandLineParams()[1]
  vcf: VCF
  ovcf: VCF
  cmds: File
  ofa: File
  reads = "bwa-hg002-dups.fa"

if not open(vcf, vcf_path):
  quit "couldnt open vcf:" & vcf_path

if not open(ovcf, "-", mode="w"):
  quit "couldnt open vcf to stdout"

if not open(ofa, reads, fmWrite):
  quit "couldnt open reads file"

ovcf.header= vcf.header

var svtype: string
for v in vcf:
  if v.info.get("SVTYPE", svtype) != Status.OK:
    quit "couldn't get SVTYPE"
  if svtype != "INS": continue

  #assert ovcf.write_variant(v)
  ofa.write(&">{v.CHROM}_{v.start}_{v.stop}_{v.ID}\n{v.ALT[0]}\n")

ofa.close()

var prefix = "bwa-hg002-dups"
if not open(cmds, prefix & ".sh", fmWrite):
  quit "bad"

cmds.write(&"""set -euo pipefail; bwa mem -R "@RG\tID:hg002\tSM:hg002" -t 3 {fasta_path} {reads} | samtools sort -o {prefix}.bam && samtools index {prefix}.bam""")
cmds.close
let (outp, errc) = execCmdEx(&"""bash {prefix}.sh""")
if outp.len > 0:
  stderr.write_line outp
doAssert errc == 0

var bam:Bam
open(bam, prefix & ".bam")

var mappings = initTable[string, Record](32768)

for aln in bam:
    # TODO: handle SA.
    if aln.flag.secondary or aln.flag.supplementary or aln.flag.unmapped: continue
    if aln.mapping_quality == 0: continue
    var nm = tag[int](aln, "NM")
    if not nm.isNone and nm.get.int > 5: continue

    mappings[aln.qname] = aln.copy()

vcf.close()
if not open(vcf, vcf_path):
  quit "couldnt reopen vcf:" & vcf_path
discard ovcf.write_header

var dup = "DUP"
var svt = ""
for v in vcf:
  var key = &"{v.CHROM}_{v.start}_{v.stop}_{v.ID}"
  if v.info.get("SVTYPE", svt) != Status.OK:
      quit "couldn't get SVTYPE"
  if svt != "DEL" and svt != "INS": continue
  if v.FILTER != "PASS": continue
  if key in mappings and svt == "INS":
    var aln = mappings[key]
    #stderr.write_line $aln.start, " ", $(v.start + 1)
    var L = v.ALT[0].len.float64
    var alnL = aln.stop - aln.start
    if abs(aln.start - v.start) < 3 and abs(alnL - L.int) < 3:
      if v.info.set("SVTYPE", dup) != Status.OK:
        quit "couldn't set type to dup"
      var stop = @[int32(v.start + v.ALT[0].len + 1)]
      if v.info.set("END", stop) != Status.OK:
        quit "couldn't set end of dup"

  doAssert ovcf.write_variant(v)

ovcf.close()
