# duphold: uphold your DUP and DEL calls

SV callers like [lumpy](https://github.com/arq5x/lumpy) look at split-reads and pair distances to find structural variants.
This tool is a fast way to add depth information to those calls. This can be used as additional
information for filtering variants; for example **we will be skeptical of deletion calls that
do not have lower than average coverage** compared to regions with similar gc-content.


`duphold` takes a **bam/cram**, a **VCF/BCF** of SV calls, and a **fasta** reference and it updates the FORMAT field for a
single sample with:

+ **DHFC**: fold-change for the variant depth *relative to the rest of the chromosome* the variant was found on
+ **DHBZ**: z-score for the variant depth *relative to bins in the genome with similar GC-content*.
+ **DHBFC**: fold-change for the variant depth *relative to bins in the genome with similar GC-content*.
+ **DHD**: rapid change in depth at one of the break-points (1 for higher (DUP). 0 for no or conflicting changes. -1 for drop (DUP), 2 or -2 for both break points)

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

## Examples

#### Duplication

Here is a duplication with clear break-points (`DHD`) and clear
change in depth (`DHBFC`)

![image](https://user-images.githubusercontent.com/1739/45895409-5a224080-bd8e-11e8-844f-e7ffc13c7972.png "example IGV screenshot")

`duphold` annotated this with

+ **DHD**: 2
+ **DHBFC**: 1.79

where together these indicate rapid (DUP-like) change in depth at the break-points and a coverage that 1.79 times higher than the mean for the genome--again indicative of a DUP. Together, these recapitulate (or anticipate) what we see on visual inspection.

#### Deletion

A clear deletion will have rapid drop in depth at the left and increase in depth at the right and a lower mean coverage.

![image](https://user-images.githubusercontent.com/1739/45895721-2dbaf400-bd8f-11e8-88b3-9fd5a90ef39e.png)

`duphold` annotated this with:

+ **DHD**: -2
+ **DHBFC**: 0.6

These indicate that both break-points are consistent with a deletion and that the coverage is ~60% of expected. So this is a clear deletion.

#### BND

when lumpy decides that a cluster of evidence does not match a DUP or DEL or INV, it creates a BND with 2 lines in the VCF. Sometimes these
are actual deletions. For example:

![image](https://user-images.githubusercontent.com/1739/45906495-987d2700-bdb1-11e8-8ba5-eacdf8221f68.png)

shows where a deletion is bounded by 2 DEL calls. `duphold` annotates this with:

+ **DHD**: -2
+ **DHBFC**: 0.01

indicating a homozygous deletion with clear break-points.

## Acknowledgements

I stole the idea of annotating SVs with depth-change from Ira Hall.
