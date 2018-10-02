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

samtools index hg002.cram

smoove call --genotype -o lf-dev/ -d -x -p 5 -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -n HG002 hg002.cram
DONE

mkdir -p eval/
rm -rf eval/esvr
rm -rf eval/lumpy-filter

python ~/Projects/src/truvari/truvari.py -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c lf-dev/HG002-smoove.genotyped.DEL.vcf.gz -o eval/unfiltered/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.95


bcftools view -i '(FMT/DHD[0] < 0)' lf-dev/HG002-smoove.genotyped.DEL.vcf.gz -O z -o lf-dev/HG002-smoove.genotyped.DEL.duphold-DHD.vcf.gz
tabix lf-dev/HG002-smoove.genotyped.DEL.duphold-DHD.vcf.gz
python ~/Projects/src/truvari/truvari.py -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c lf-dev/HG002-smoove.genotyped.DEL.duphold-DHD.vcf.gz -o eval/dhd/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.95

bcftools view -i '(FMT/DHBFC[0] < 0.75)' lf-dev/HG002-smoove.genotyped.DEL.vcf.gz -O z -o lf-dev/HG002-smoove.genotyped.DEL.duphold-DHBFC.vcf.gz
tabix lf-dev/HG002-smoove.genotyped.DEL.duphold-DHBFC.vcf.gz
python ~/Projects/src/truvari/truvari.py -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c lf-dev/HG002-smoove.genotyped.DEL.duphold-DHBFC.vcf.gz -o eval/dhbfc/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.95

bcftools view -i '(FMT/DHBFC[0] < 0.75) && (FMT/DHD[0] < 0)' lf-dev/HG002-smoove.genotyped.DEL.vcf.gz -O z -o lf-dev/HG002-smoove.genotyped.DEL.duphold-both.vcf.gz
tabix lf-dev/HG002-smoove.genotyped.DEL.duphold-both.vcf.gz
python ~/Projects/src/truvari/truvari.py -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c lf-dev/HG002-smoove.genotyped.DEL.duphold-both.vcf.gz -o eval/both/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.95


rm -rf eval/stringent
bcftools view -i '(FMT/DHBFC[0] < 0.65) && (FMT/DHD[0] == -2)' lf-dev/HG002-smoove.genotyped.DEL.vcf.gz -O z -o lf-dev/HG002-smoove.genotyped.DEL.duphold-stringent.vcf.gz
tabix lf-dev/HG002-smoove.genotyped.DEL.duphold-stringent.vcf.gz
python ~/Projects/src/truvari/truvari.py -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c lf-dev/HG002-smoove.genotyped.DEL.duphold-stringent.vcf.gz -o eval/stringent/ --passonly --pctsim=0  -r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.95

python figure1.py eval/{unfiltered,dhd,dhbfc,both,stringent}/summary.txt


~/Projects/src/samplot/src/samplot_vcf.sh -o samplot-fp/ -v eval/stringent/fp.vcf -r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -S ~/Projects/src/samplot/src/samplot.py hg002.cram
# used: http://home.chpc.utah.edu/~u6000771/samplot-fp/DEL_4_36469923-36470237.png


bcftools view -i '(FMT/DHBFC[0] > 0.85) && (FMT/DHD[0] == 0)' lf-dev/HG002-smoove.genotyped.DEL.vcf.gz -O z -o lf-dev/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz
tabix lf-dev/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz
rm -rf eval-no-support
python ~/Projects/src/truvari/truvari.py -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c \
	lf-dev/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz -o eval-no-support --passonly --pctsim=0 \
   	-r 20 --giabreport -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa --no-ref \
	--includebed HG002_SVs_Tier1_v0.6.bed -O 0.95

rm -rf samplot-tp-missed
mkdir samplot-tp-missed
~/Projects/src/samplot/src/samplot_vcf.sh -o samplot-tp-missed/ -v eval-no-support/tp-call.vcf \
	-r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa \
	-S ~/Projects/src/samplot/src/samplot.py hg002.cram

## speed
 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
zcat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | cut -f 1-10 | perl -pe 's/HG00096/HG002/' | bgzip -c > ALL.wgs.merged-svs.vcf.gz

time duphold  -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t 3 -b hg002.cram -o all.vcf -v ALL.wgs.merged-svs.vcf.gz


for n in 100 1000 5000 10000 20000 40000 60000; do
time duphold  -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t 3 -b hg002.cram -o all.vcf -v ALL.wgs.merged-svs.vcf.gz.$n.vcf.gz
done
