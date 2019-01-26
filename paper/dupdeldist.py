import cyvcf2
import sys
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
sns.set_style('white')

colors = sns.color_palette()

metrics = {"DUP": [[], [], []],
           "DEL": [[], [], []]}

metric_dup = "DHBFC"
metric_del = "DHFFC"

gcs = {"DUP":[], "DEL": []}

for v in cyvcf2.VCF(sys.argv[1], gts012=True):
    if v.FILTER is not None: continue
    svtype = v.INFO.get("SVTYPE")
    if not svtype in ("DEL", "DUP"): continue
    gt = v.gt_types[0]

    gcs[svtype].append(v.INFO.get("GCF"))
    #if v.INFO.get("GCF") == 0.0: continue

    metric = metric_dup if svtype == "DUP" else metric_del

    val = v.format(metric)[0][0]
    if np.isnan(val):
        val = -1.0
        if gt == 0: continue
    metrics[svtype][gt].append(min(3, val))

fix, ax = plt.subplots(2, 2)
ax[0, 0].hist(metrics["DUP"][0], 40, label="0/0", alpha=0.7, color=colors[0])
ax[0, 0].hist(metrics["DUP"][1], 40, label="0/1", alpha=0.7, color=colors[1])
ax[0, 0].hist(metrics["DUP"][2], 40, label="1/1", alpha=0.7, color=colors[2])
ax[0, 0].set_ylabel("Count")
ax[0, 0].set_xlim(0, 2.5)
ax[0, 0].set_xlabel(metric_dup)
ax[0, 0].text(0.6, 0.24, "Duplications", transform=ax[0, 0].transAxes)
ax[0, 0].legend()

ax[1, 0].hist(metrics["DEL"][0], 40, label="0/0", alpha=0.7, color=colors[0])
ax[1, 0].hist(metrics["DEL"][1], 40, label="0/1", alpha=0.7, color=colors[1])
ax[1, 0].hist(metrics["DEL"][2], 40, label="1/1", alpha=0.7, color=colors[2])
ax[1, 0].text(0.6, 0.24, "Deletions", transform=ax[1, 0].transAxes)
ax[1, 0].set_xlim(0, 2.5)
ax[1, 0].set_xlabel(metric_del)
ax[1, 0].set_ylabel("Count")

from sklearn.metrics import auc, roc_curve

L = ["xx", "0/1", "1/1"]

for i, ev in enumerate(("DUP", "DEL")):

    for alts in (1, 2):
        truth = [0] * len(metrics[ev][0]) + [1] * len(metrics[ev][alts])
        scores = metrics[ev][0] + metrics[ev][alts]
        if ev == "DEL":
            scores = [-s for s in scores]
        fpr, tpr, rscores = roc_curve(truth, scores)
        if ev == "DEL":
            rscores = -rscores
            idx = np.searchsorted(rscores, 0.7)
        else:
            idx = np.searchsorted(-rscores, -1.3)
        ax[i, 1].plot(fpr, tpr, color=colors[alts])
        ax[i, 1].plot([0, 1], [0, 1], '--', color='gray')
        ax[i, 1].text(0.4, 0.04 + (2 - alts) * 0.1, "0/0 vs %s AUC: %.2f"
                % (L[alts], auc(fpr, tpr)),
                color=colors[alts])
        ax[i, 1].plot([fpr[idx]], [tpr[idx]], marker='o', color=colors[alts])
    ax[i, 1].text(0.4, 0.04 + 2 * 0.1, "Duplications" if ev == "DUP" else "Deletions")

ax[0, 1].set_xlabel("1 - Sensitivity")
ax[1, 1].set_xlabel("1 - Sensitivity")
ax[1, 1].set_ylabel("Specificity")
ax[0, 1].set_ylabel("specificity")

plt.tight_layout()
plt.show()
plt.close()
