import time
import toolshed as ts

svt = "zcat {vcf} | svtyper -o t -B hg002.cram -T ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa"
dph = "duphold -f ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa -t {threads} -b hg002.cram -o all.vcf -v {vcf}"

for n in [100, 1000, 5000, 10000, 20000, 35000, 50000]:

    times = [float('inf'), float('inf'), float('inf')]

    for i in range(5):
        t0 = time.time()

        list(ts.nopen("|" + svt.format(vcf="ALL.wgs.merged-svs.vcf.gz." + str(n) + ".vcf.gz")))

        st = time.time() - t0
        times[0] = min(times[0], st)



        t0 = time.time()

        list(ts.nopen("|" + dph.format(threads=0, vcf="ALL.wgs.merged-svs.vcf.gz." + str(n) + ".vcf.gz")))
        dh0 = time.time() - t0
        times[1] = min(times[1], dh0)

        t0 = time.time()
        list(ts.nopen("|" + dph.format(threads=3, vcf="ALL.wgs.merged-svs.vcf.gz." + str(n) + ".vcf.gz")))
        dh3 = time.time() - t0
        times[2] = min(times[2], dh3)

    print "%d\t%.2f\t%.2f\t%.2f" % (n, times[0], times[1], times[2])

