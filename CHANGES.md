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
