set -euo pipefail

#ncftp -R ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/ .

#for f in $(find . -type f -name "*_R1_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 14 -c > hg002_R1.fastq.gz &
#for f in $(find . -type f -name "*_R2_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 14 -c > hg002_R2.fastq.gz &

#wait

ref=~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa
<<DONE
~/Projects/src/bwa/bwa mem -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tPU:HG002\tLB:HG002' -t 32 $ref hg002_R1.fastq.gz hg002_R2.fastq.gz \
	| samblaster \
	| samtools sort -@ 4 -m 15G --output-fmt CRAM --reference $ref -o hg002.cram

bcftools view -i SVTYPE="DEL" -O z -o HG002_SVs_Tier1_v0.6.DEL.vcf.gz HG002_SVs_Tier1_v0.6.vcf.gz
tabix HG002_SVs_Tier1_v0.6.DEL.vcf.gz
bcftools view -i 'REPTYPE="DUP"' -O z -o HG002_SVs_Tier1_v0.6.DUP.vcf.gz HG002_SVs_Tier1_v0.6.vcf.gz
tabix HG002_SVs_Tier1_v0.6.DUP.vcf.gz


samtools index hg002.cram

smoove call -F --genotype -o lf-dev2/ -x -p 5 -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -n HG002 hg002.cram
duphold -v lf-dev2/HG002-smoove.genotyped.vcf.gz -o lf-dev2/HG002-smoove.genotyped.duphold.vcf.gz -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t 4  -b hg002.cram

bcftools view -i 'SVTYPE="DEL"' lf-dev2/HG002-smoove.genotyped.duphold.vcf.gz -O z -o lf-dev2/HG002-smoove.genotyped.DEL.vcf.gz
tabix lf-dev2/HG002-smoove.genotyped.DEL.vcf.gz
DONE

ev=DEL
filt="< 0.7"

eval=eval-$ev
mkdir -p $eval/
rm -rf $eval/*

sizemax=15000000
sizemin=300


python ~/Projects/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b HG002_SVs_Tier1_v0.6.$ev.vcf.gz -c lf-dev2/HG002-smoove.genotyped.$ev.vcf.gz -o $eval/unfiltered/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.6

bcftools view -i "(FMT/DHBFC[0] $filt)" lf-dev2/HG002-smoove.genotyped.$ev.vcf.gz -O z -o lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHBFC.vcf.gz
tabix lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHBFC.vcf.gz
python ~/Projects/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b HG002_SVs_Tier1_v0.6.$ev.vcf.gz -c lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHBFC.vcf.gz -o $eval/dhbfc/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.6

bcftools view -i "(FMT/DHFFC[0] $filt)" lf-dev2/HG002-smoove.genotyped.$ev.vcf.gz -O z -o lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHFFC.vcf.gz
tabix lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHFFC.vcf.gz
python ~/Projects/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b HG002_SVs_Tier1_v0.6.$ev.vcf.gz -c lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHFFC.vcf.gz -o $eval/dhffc/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.6

bcftools view -i "(FMT/DHFC[0] $filt)" lf-dev2/HG002-smoove.genotyped.$ev.vcf.gz -O z -o lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHFC.vcf.gz
tabix lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHFC.vcf.gz
python ~/Projects/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b HG002_SVs_Tier1_v0.6.$ev.vcf.gz -c lf-dev2/HG002-smoove.genotyped.$ev.duphold-DHFC.vcf.gz -o $eval/dhfc/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.6

python figure1.py eval-DEL/{unfiltered,dhfc,dhbfc,dhffc}/summary.txt
exit

#python figure1.py eval-DUP/{unfiltered,dhfc,dhbfc,dhffc}/summary.txt
export SAMPLOT_COVERAGE_ONLY=FALSE
rm -rf ~/public_html/samplot-fps/
~/Projects/src/samplot/src/samplot_vcf.sh -O png -o ~/public_html/samplot-fps/ -v $eval/dhffc/fp.vcf -r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -S ~/Projects/src/samplot/src/samplot.py hg002.cram
export SAMPLOT_COVERAGE_ONLY=FALSE
rm -rf ~/public_html/samplot-tp/
~/Projects/src/samplot/src/samplot_vcf.sh -O pdf -o ~/public_html/samplot-tp/ -v $eval/dhffc/tp-call.vcf -r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -S ~/Projects/src/samplot/src/samplot.py hg002.cram

### DUPS
#wget http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam.sniffles1kb_auto_noalts.vcf.gz

bcftools view -h all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam.sniffles1kb_auto_noalts.vcf.gz | grep -v '#CHROM' > sniffles.dups.vcf
zgrep $'INFO\tFORMAT' all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam.sniffles1kb_auto_noalts.vcf.gz | perl -pe 's/(.*)\t.+$/\1\tHG002/' >> sniffles.dups.vcf
bcftools view --apply-filters PASS all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam.sniffles1kb_auto_noalts.vcf.gz | awk '$5 == "<DUP>"' >> sniffles.dups.vcf
gsort sniffles.dups.vcf  ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa.fai | bgzip -c > sniffles.dups.vcf.gz
tabix -f sniffles.dups.vcf.gz

duphold -v sniffles.dups.vcf.gz -o duphold.annotated.dups.vcf.gz -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t 6 -b hg002.cram

unset SAMPLOT_COVERAGE_ONLY
~/Projects/src/samplot/src/samplot_vcf.sh -O png -o ~/public_html/samplot-dups/ -v duphold.annotated.dups.vcf.gz -r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -S ~/Projects/src/samplot/src/samplot.py hg002.cram



### MISSED (False Negative)
eval=missed
bcftools view -i '(FMT/DHFFC[0] > 0.85)' lf-dev2/HG002-smoove.genotyped.DEL.vcf.gz -O z -o lf-dev2/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz
tabix lf-dev2/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz
rm -rf $eval-no-support
python ~/Projects/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c \
	lf-dev2/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz -o $eval-no-support --passonly --pctsim=0 \
   	-r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref \
	--includebed HG002_SVs_Tier1_v0.6.bed -O 0.95

rm -rf ~/public_html/samplot-tp-missed
~/Projects/src/samplot/src/samplot_vcf.sh -o ~/public_html/samplot-tp-missed/ -v $eval-no-support/tp-call.vcf \
	-r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa \
	-S ~/Projects/src/samplot/src/samplot.py hg002.cram

## speed
 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
zcat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | cut -f 1-10 | perl -pe 's/HG00096/HG002/' | bgzip -c > ALL.wgs.merged-svs.vcf.gz

time duphold  -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t 3 -b hg002.cram -o all.vcf -v ALL.wgs.merged-svs.vcf.gz


for n in 100 1000 5000 10000 20000 40000 60000; do
time duphold  -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t 3 -b hg002.cram -o all.vcf -v ALL.wgs.merged-svs.vcf.gz.$n.vcf.gz
done
