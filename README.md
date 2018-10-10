[![Build Status](https://travis-ci.org/brentp/duphold.svg?branch=master)](https://travis-ci.org/brentp/duphold)

# duphold: uphold your DUP and DEL calls

SV callers like [lumpy](https://github.com/arq5x/lumpy) look at split-reads and pair distances to find structural variants.
This tool is a fast way to add depth information to those calls. This can be used as additional
information for filtering variants; for example **we will be skeptical of deletion calls that
do not have lower than average coverage** compared to regions with similar gc-content.

In addition, `duphold` will annotate the SV vcf with information from a SNP/Indel VCF. For example, **we will not
believe a large deletion that has many heterozygote SNP calls**.


`duphold` takes a **bam/cram**, a **VCF/BCF** of SV calls, and a **fasta** reference and it updates the FORMAT field for a
single sample with:

+ **DHFC**: fold-change for the variant depth *relative to the rest of the chromosome* the variant was found on
+ **DHBFC**: fold-change for the variant depth *relative to bins in the genome with similar GC-content*.
+ **DHFFC**: fold-change for the variant depth *relative to **F**lanking regions*.

If a SNP/Indel VCF/BCF is given, `duphold` will annotate each DEL/DUP call with (see below for more detail on what it does):

+ **DHET**: counts of SNP heterozygotes in the SV supporting: [0] a normal heterozygote, [1] a triploid heterozygote.
            for a DUP, we expect most hets to have an allele balance closer to 0.33 or 0.67 than to 0.5. A good heterozygous
            DUP will have larger values of [1], than [0], though it's possible there are no HETs in small events.

+ **DHHU**: counts of [0] Hom-ref, [1] Hom-alt, [2] Unknown variants in the event. A heterozygous deletion may have more hom-alt SNP calls.
            A homozygous deletion may have only unknown SNP calls.

It also adds **GCF** to the INFO field indicating the fraction of G or C bases in the variant.


## SNP/Indel annotation

**NOTE** it is strongly recommended to use BCF for the `--snp` argument as otherwise VCF parsing will be a bottleneck.

+ A DEL call with many HETs is unlikely to be valid.
+ A DUP call that has many HETs that have a 0.5 allele balance is unlikely to be valid.

When the user specifies a `--snp` VCF, `duphold` finds the appropriate sample in that file and extracts high (> 20) quality, bi-allelic
SNP calls. For each chromosome, it will store a minimal (low-memory representation) in a sorted data-structure for fast access. It will
then query this data structure for each SV and count the number of heterozygotes supporting a diploid HET (allele balance close to 0.5)
or a triploid HET (allele balance close to 0.33 or 0.67) into `DHET`. It will store the number of Hom-Ref, Hom-Alt, Unnkown calls in
`DHHU`.

## Performance

### Speed

`duphold` runtime depends almost entirely on how long it takes to parse the BAM/CRAM files; it is relatively independent of the number of variants evaluated. It will also run quite a bit faster on CRAM than on BAM. It can be < 20 minutes of CPU time for a 30X CRAM.

### Accuracy

coming soon.

## Install

`duphold` is distributed as a binary [here](https://github.com/brentp/duphold/releases/latest) and requires libhts.so in standard locations or indicated with `LD_LIBRARY_PATH`.



## Usage

```
duphold -s $gatk_vcf -t 4 -v $svvcf -b $cram -f $fasta -o $output.bcf
duphold --snp $gatk_bcf --threads 4 --vcf $svvcf --bam $cram --fasta $fasta --output $output.bcf
```

`--snp` can be a multi-sample VCF/BCF. `duphold` will be much faster with a BCF, especially if
the snp/indel file contains many (>20 or so) samples.

the threads are decompression threads so increasing up to about 4 works.

## Examples

#### Duplication

Here is a duplication with clear change in depth (`DHBFC`)

![image](https://user-images.githubusercontent.com/1739/45895409-5a224080-bd8e-11e8-844f-e7ffc13c7972.png "example IGV screenshot")

`duphold` annotated this with

+ **DHBFC**: 1.79

where together these indicate rapid (DUP-like) change in depth at the break-points and a coverage that 1.79 times higher than the mean for the genome--again indicative of a DUP. Together, these recapitulate (or anticipate) what we see on visual inspection.

#### Deletion

A clear deletion will have rapid drop in depth at the left and increase in depth at the right and a lower mean coverage.

![image](https://user-images.githubusercontent.com/1739/45895721-2dbaf400-bd8f-11e8-88b3-9fd5a90ef39e.png)

`duphold` annotated this with:

+ **DHBFC**: 0.6

These indicate that both break-points are consistent with a deletion and that the coverage is ~60% of expected. So this is a clear deletion.

#### BND

when lumpy decides that a cluster of evidence does not match a DUP or DEL or INV, it creates a BND with 2 lines in the VCF. Sometimes these
are actual deletions. For example:

![image](https://user-images.githubusercontent.com/1739/45906495-987d2700-bdb1-11e8-8ba5-eacdf8221f68.png)

shows where a deletion is bounded by 2 BND calls. `duphold` annotates this with:

+ **DHBFC**: 0.01

indicating a homozygous deletion with clear break-points.

## Acknowledgements

I stole the idea of annotating SVs with depth-change from Ira Hall.
