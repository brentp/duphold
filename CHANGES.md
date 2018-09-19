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
