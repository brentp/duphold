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
