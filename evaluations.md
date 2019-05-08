evaluations are preformed on hg002 from genome in a bottle truthset. see [readme](https://github.com/brentp/duphold#accuracy)

for events > 300bp using output from [smoove](https://github.com/brentp/smoove) v0.2.3 with the --noextrafilters flag:

| method      |   FDR |   FN |   FP |   TP-call |   precision |   recall |   recall-% |    FP-% |
|:------------|------:|-----:|-----:|----------:|------------:|---------:|-----------:|--------:|
| unfiltered  | 0.054 |  276 |   86 |      1496 |       0.946 |    0.844 |    100.000 | 100.000 |
| DHBFC < 0.7 | 0.018 |  298 |   27 |      1474 |       0.982 |    0.832 |     98.529 |  31.395 |
| DHFFC < 0.7 | 0.021 |  289 |   32 |      1483 |       0.979 |    0.837 |     99.131 |  37.209 |


When run with the default filters, the DHFFC loses 1 TP (1482 vs 1483) and reduces the FPs from 32 to 29):

| method     |   FDR |   FN |   FP |   TP-call |   precision |   recall |   recall-% |    FP-% |
|:-----------|------:|-----:|-----:|----------:|------------:|---------:|-----------:|--------:|
| unfiltered | 0.044 |  282 |   69 |      1490 |       0.956 |    0.841 |    100.000 | 100.000 |
| dhfc       | 0.017 |  295 |   26 |      1477 |       0.983 |    0.834 |     99.128 |  37.681 |
| dhbfc      | 0.015 |  304 |   23 |      1468 |       0.985 |    0.828 |     98.523 |  33.333 |
| dhffc      | 0.019 |  290 |   29 |      1482 |       0.981 |    0.836 |     99.463 |  42.029 |


