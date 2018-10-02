import numpy as np

data = """\
100     12.13   367.19  213.85
1000    39.48   381.16  214.17
5000    158.51  381.67  215.11
10000   313.93  373.39  222.55
20000   614.61  373.57  222.99
35000   1065.35 377.65  220.66
50000   1530.08 376.60  225.07\
""".split("\n")

data = [[int(x[0])] + map(float, x[1:]) for x in (l.split() for l in data)]

xs = np.array([x[0] for x in data])
data = np.array([x[1:] for x in data])

from matplotlib import pyplot as plt
import seaborn as sns

plt.plot(xs, data[:, 0], ls='--', marker='o', label="svtyper")
plt.plot(xs, data[:, 1], ls='--', marker='o', label="duphold (1 thread)")
plt.plot(xs, data[:, 2], ls='--', marker='o', label="duphold (3 threads)")

plt.xlabel("Number of variants")
plt.ylabel("Time (seconds)")
plt.legend()
plt.show()

