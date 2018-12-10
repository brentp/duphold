v0.1.2 (dev)
============
+ remove PCRE dependency.
+ reduce logging output.
+ support un-indexed VCFs without contig definitions.

v0.1.1
======
+ fix bug when later chroms had larger values (see: 9f88400 and thanks @raul-w for reporting and providing a test-case).
+ adjust calculation of global and per-gc bin coverage to use only non-zero bases. this improves the ratio when comparing
  to the depth inside of events (which will still count zero-coverage bases) in chromosomes with sparse coverage. (thanks @raul-w
  for suggesting).

v0.1.0
======
+ reduce memory usage in discordant calculation.

v0.0.9
======
+ drop DHD in favor of DHFFC for (F)lank fold change. 
+ make changes to support long, single-end reads.
+ fix bug in hts-nim that made duphold grab the wrong sample info from a SNP BCF (not vcf)

v0.0.8
======
+ skip duphold calculation for distant or interchromosomal BND's
+ don't decode CRAM qname. This makes duphold use ~15% less CPU for CRAM
+ allow annotating SVs with SNP calls as a additional validation of each event.

v0.0.7
======
+ remove DHBZ for the z-score. this was less useful than the fold-change values.
+ greatly improve DHD. the cutoffs for this are now data-driven so that we expect a low false-positive rate
  of ends that are called to have a rapid change in depth. this requires an extra pass over each chromosome that
  will add about 90 seconds of runtime per (human or similar-sized) genome.
+ `DHBFC` is now based on the median instead of mean so it is less susceptible to outliers and off-center calls
+ fix bug that in rare cases resulted in exiting with "how"
+ also report `DHSP` for number of supporting read-pairs for the event. This is usually redundant with SVTyper info, but might
  be more accurate in some cases.

v0.0.6
======
+ `duphold` will now annotate the space between BND's that occur within 20MB on the same chromosome. Sometimes events that are obvious
   deletions by their change in coverage will be called as a pair of BNDs. This will help prioritize those that are accompanied by a
   depth change.

v0.0.5
======
+ fix bug that occurred with SVs near either end of chromosome (brentp/duphold#2). thanks Brad for reporting and providing a test-case.

v0.0.4
======
+ small bug-fixes
+ duphold now uses about half as much memory
+ add -d/--drop flag which will drop all samples from a VCF except the
  one matching the sample in --bam. this simplifies per-sample 
  parallelization followed by merge.


v0.0.3
======
+ fix bug when annotatign multi-sample vcf (brentp/duphold#2)
+ get sample name from bam read group info. **NOTE** that this changes the command line parameters by removing the `--sample` argument.
