set -euo pipefail

#ncftp -R ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/ .

#for f in $(find . -type f -name "*_R1_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 14 -c > hg002_R1.fastq.gz &
#for f in $(find . -type f -name "*_R2_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 14 -c > hg002_R2.fastq.gz &

#wait

ref=/data/human/g1k_v37_decoy.fa
CRAM=/data/human/hg002.cram
truth=/data/human/HG002_SVs_Tier1_v0.6.vcf.gz
truth_del=/data/human/HG002_SVs_Tier1_v0.6.DEL.vcf.gz

bcftools view -i 'SVTYPE="DEL"' -O z -o $truth_del --apply-filters "PASS,." $truth
tabix $truth_del
<<DONE

pip install progressbar2 python-levenshtein "intervaltree<2.1.0" pyvcf pyfaidx pysam python-dateutil

~/Projects/src/bwa/bwa mem -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tPU:HG002\tLB:HG002' -t 32 $ref hg002_R1.fastq.gz hg002_R2.fastq.gz \
	| samblaster \
	| samtools sort -@ 4 -m 15G --output-fmt CRAM --reference $ref -o hg002.cram

bcftools view -i 'REPTYPE="DUP"' -O z -o HG002_SVs_Tier1_v0.6.DUP.vcf.gz HG002_SVs_Tier1_v0.6.vcf.gz
tabix HG002_SVs_Tier1_v0.6.DUP.vcf.gz

DONE

set -euo pipefail

ODIR=lf-dev-docker
#docker run -v "$(pwd):/work" -v $(pwd)/src:/src/ -v $(dirname $CRAM):$(dirname $CRAM) -v $(dirname $ref):$(dirname $ref) brentp/smoove:v0.2.3 smoove call -F -d --genotype -o $ODIR/ -x -p 1 -f $ref -n HG002 $CRAM


ODIR=smoove-dev
rm -rf $ODIR
mkdir -p $ODIR
smoove call -d --genotype -o $ODIR/ -x -p 1 -f $ref -n HG002 $CRAM

bcftools view -i 'SVTYPE="DEL"' $ODIR/HG002-smoove.genotyped.vcf.gz -O z -o $ODIR/HG002-smoove.genotyped.DEL.vcf.gz
tabix $ODIR/HG002-smoove.genotyped.DEL.vcf.gz

ev=DEL
filt="< 0.7"

eval=eval-$ev
mkdir -p $eval/
rm -rf $eval/*

sizemax=15000000
sizemin=300
bed=/data/human/HG002_SVs_Tier1_v0.6.bed


python ~/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b $truth_del -c $ODIR/HG002-smoove.genotyped.$ev.vcf.gz -o $eval/unfiltered/ --passonly --pctsim=0  -r 20 --giabreport -f $ref --no-ref --includebed $bed -O 0.6

bcftools view -i "(FMT/DHBFC[0] $filt)" $ODIR/HG002-smoove.genotyped.$ev.vcf.gz -O z -o $ODIR/HG002-smoove.genotyped.$ev.duphold-DHBFC.vcf.gz
tabix $ODIR/HG002-smoove.genotyped.$ev.duphold-DHBFC.vcf.gz
python ~/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b $truth_del -c $ODIR/HG002-smoove.genotyped.$ev.duphold-DHBFC.vcf.gz -o $eval/dhbfc/ --passonly --pctsim=0  -r 20 --giabreport -f $ref --no-ref --includebed $bed -O 0.6

bcftools view -i "(FMT/DHFFC[0] $filt)" $ODIR/HG002-smoove.genotyped.$ev.vcf.gz -O z -o $ODIR/HG002-smoove.genotyped.$ev.duphold-DHFFC.vcf.gz
tabix $ODIR/HG002-smoove.genotyped.$ev.duphold-DHFFC.vcf.gz
python ~/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b $truth_del -c $ODIR/HG002-smoove.genotyped.$ev.duphold-DHFFC.vcf.gz -o $eval/dhffc/ --passonly --pctsim=0  -r 20 --giabreport -f $ref --no-ref --includebed $bed -O 0.6

bcftools view -i "(FMT/DHFC[0] $filt)" $ODIR/HG002-smoove.genotyped.$ev.vcf.gz -O z -o $ODIR/HG002-smoove.genotyped.$ev.duphold-DHFC.vcf.gz
tabix $ODIR/HG002-smoove.genotyped.$ev.duphold-DHFC.vcf.gz
python ~/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b $truth_del -c $ODIR/HG002-smoove.genotyped.$ev.duphold-DHFC.vcf.gz -o $eval/dhfc/ --passonly --pctsim=0  -r 20 --giabreport -f $ref --no-ref --includebed $bed -O 0.6

python paper/table1.py eval-DEL/{unfiltered,dhfc,dhbfc,dhffc}/summary.txt
exit

#python figure1.py eval-DUP/{unfiltered,dhfc,dhbfc,dhffc}/summary.txt
export SAMPLOT_COVERAGE_ONLY=FALSE
rm -rf ~/public_html/samplot-fps/
~/src/samplot/src/samplot_vcf.sh -O png -o ~/public_html/samplot-fps/ -v $eval/dhffc/fp.vcf -r $ref -S ~/src/samplot/src/samplot.py $CRAM
export SAMPLOT_COVERAGE_ONLY=FALSE
rm -rf ~/public_html/samplot-tp/
~/src/samplot/src/samplot_vcf.sh -O pdf -o ~/public_html/samplot-tp/ -v $eval/dhffc/tp-call.vcf -r $ref -S ~/src/samplot/src/samplot.py $CRAM

eval=missed
bcftools view -i '(FMT/DHFFC[0] > 0.85)' $ODIR/HG002-smoove.genotyped.DEL.vcf.gz -O z -o $ODIR/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz
tabix $ODIR/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz
rm -rf $eval-no-support
python ~/src/truvari/truvari.py --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b $truth_del -c \
	$ODIR/HG002-smoove.genotyped.DEL.duphold-no-support.vcf.gz -o $eval-no-support --passonly --pctsim=0 \
   	-r 20 --giabreport -f $ref --no-ref \
	--includebed $bed -O 0.95

rm -rf ~/public_html/samplot-tp-missed
~/src/samplot/src/samplot_vcf.sh -o ~/public_html/samplot-tp-missed/ -v $eval-no-support/tp-call.vcf \
	-r $ref \
	-S ~/src/samplot/src/samplot.py $CRAM

### DUPS

duphold -v dups.vcf.gz -o duphold.annotated.dups.vcf.gz -f $ref -t 6 -b $CRAM


## speed
 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
zcat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | cut -f 1-10 | perl -pe 's/HG00096/HG002/' | bgzip -c > ALL.wgs.merged-svs.vcf.gz

time duphold  -f $ref -t 3 -b $CRAM -o all.vcf -v ALL.wgs.merged-svs.vcf.gz

for n in 100 1000 5000 10000 20000 40000 60000; do
time duphold  -f $ref -t 3 -b $CRAM -o all.vcf -v ALL.wgs.merged-svs.vcf.gz.$n.vcf.gz
done


######## distribution of del, dup calls for truth:

struth=paper-results/giab-with-fake.vcf.gz
nim c -d:release paper/giab_ins_to_dup.nim 
# find what appear to be DUP calls
./paper/giab_ins_to_dup /data/human/HG002_SVs_Tier1_v0.6.vcf.gz /data/human/g1k_v37_decoy.fa | bgzip -c > $struth
# insert fake, 0/0 calls.
nim c -r paper/insert_regions.nim $struth /data/human/HG002_SVs_Tier1_v0.6.bed t.vcf.gz /data/human/g1k_v37_decoy.fa
mv t.vcf.gz $struth

mkdir -p paper-results/
for f in 500 5000 50000; do
  DUPHOLD_FLANK=$f ./src/duphold -f /data/human/g1k_v37_decoy.fa -v $struth -t 3 -o paper-results/giab.duphold.flank-$f.vcf.gz -b /data/human/hg002.cram &
done
wait

for step in 100 250 5000; do
  echo "DUPHOLD_GC_STEP=$step ./src/duphold -f /data/human/g1k_v37_decoy.fa -v $struth -t 3 -o paper-results/giab.duphold.flank-step-$step.vcf.gz -b /data/human/hg002.cram"
done | gargs -p 3 "{}"
