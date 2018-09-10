# duphold: uphold your DUP and DEL calls

SV callers like [lumpy](https://github.com/arq5x/lumpy) look at split-reads and pair distances to find structural variants.
This tool is a fast way to add depth information to those calls. This can be used as additional
information for filtering variants; for example **we will be skeptical of deletion calls that
do not have lower than average coverage** compared to regions with similar gc-content.


`duphold` takes a **bam/cram**, a **VCF/BCF** of SV calls, and a **fasta** reference and it updates the FORMAT field for a
single sample with:

+ **DHZ**: z-score for the variant depth *relative to the rest of the chromosome* the variant was found on.
+ **DHFC**: fold-change for the variant depth *relative to the rest of the chromosome* the variant was found on
+ **DHBZ**: z-score for the variant depth *relative to bins in the genome with similar GC-content*.
+ **DHBFC**: fold-change for the variant depth *relative to bins in the genome with similar GC-content*.

It also adds **GCF** to the INFO field indicating the fraction of G or C bases in the variant.

## Performance

### Speed

`duphold` runtime depends almost entirely on how long it takes to parse the BAM/CRAM files; it is relatively independent of the number of variants evaluated. It will also run quite a bit faster on CRAM than on BAM. It can be < 20 minutes of CPU time for a 30X CRAM.

### Accuracy

coming soon.

## Install

`duphold` is distributed as a binary [here] and requires libhts.so in standard locations or indicated with `LD_LIBRARY_PATH`.



## Usage

```
duphold -t 4 -v $svvcf -b $cram -f $fasta -o $output.bcf
duphold --threads 4 --vcf $svvcf --bam $cram --fasta $fasta --output $output.bcf
```

the threads are decompression threads so increasing up to about 4 works.

## Acknowledgements

I stole this idea from Ira Hall.
