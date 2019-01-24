import cyvcf2
import sys
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
sns.set_style('white')

colors = sns.color_palette()

metrics = {"DHFFC": {"DUP": [[], [], []],
         "DEL": [[], [], []]},
         "DHBFC": {"DUP": [[], [], []],
         "DEL": [[], [], []]}}

metric = "DHFFC"


for v in cyvcf2.VCF(sys.argv[1], gts012=True):
    if v.FILTER is not None: continue
    svtype = v.INFO.get("SVTYPE")
    if not svtype in ("DEL", "DUP"): continue
    gt = v.gt_types[0]

    val = v.format(metric)[0][0]
    if np.isnan(val):
        val = -1.0
        if gt == 0: continue
    metrics[metric][svtype][gt].append(min(3, val))

fix, ax = plt.subplots(2, 2)
ax[0, 0].hist(metrics[metric]["DUP"][0], 40, label="0/0", alpha=0.7, color=colors[0])
ax[0, 0].hist(metrics[metric]["DUP"][1], 40, label="0/1", alpha=0.7, color=colors[1])
ax[0, 0].hist(metrics[metric]["DUP"][2], 40, label="1/1", alpha=0.7, color=colors[2])
ax[0, 0].set_ylabel("Count")
ax[0, 0].set_xlim(0, 3)
ax[0, 0].text(2, 10, "Duplications")
ax[0, 0].legend()

ax[1, 0].hist(metrics[metric]["DEL"][0], 40, label="0/0", alpha=0.7, color=colors[0])
ax[1, 0].hist(metrics[metric]["DEL"][1], 40, label="0/1", alpha=0.7, color=colors[1])
ax[1, 0].hist(metrics[metric]["DEL"][2], 40, label="1/1", alpha=0.7, color=colors[2])
ax[1, 0].text(2, 10, "Deletions")
ax[1, 0].set_xlim(0, 3)
ax[1, 0].set_xlabel(metric)
ax[1, 0].set_ylabel("Count")

from sklearn.metrics import auc, roc_curve

L = ["xx", "0/1", "1/1"]

for i, ev in enumerate(("DUP", "DEL")):

    for alts in (1, 2):
        truth = [0] * len(metrics[metric][ev][0]) + [1] * len(metrics[metric][ev][alts])
        scores = metrics[metric][ev][0] + metrics[metric][ev][alts]
        if ev == "DEL":
            scores = [-s for s in scores]
        fpr, tpr, rscores = roc_curve(truth, scores)
        if ev == "DEL":
            rscores = -rscores
            idx = np.searchsorted(rscores, 0.7)
        else:
            idx = np.searchsorted(-rscores, -1.25)
        ax[i, 1].plot(fpr, tpr, color=colors[alts])
        ax[i, 1].plot([0, 1], [0, 1], '--', color='gray')
        ax[i, 1].text(0.4, 0.04 + (2 - alts) * 0.1, "0/0 vs %s AUC: %.2f" % (L[alts], auc(fpr, tpr)),
                color=colors[alts])
        ax[i, 1].plot([fpr[idx]], [tpr[idx]], marker='o', color=colors[alts])

ax[1, 1].set_xlabel("1 - Sensitivity")
ax[1, 1].set_ylabel("Specificity")
ax[0, 1].set_ylabel("specificity")

plt.tight_layout()
plt.show()
