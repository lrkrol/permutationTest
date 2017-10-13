# permutationTest
A permutation test (or randomisation test) for MATLAB, supporting one- and two-tailed tests and capable of visualising the outcome using a histogram. Provides a p-value, the observed difference, and the effect size.

For example, the following tests two samples (n=5000 each) from normal distributions against each other, one shifted .1 to the right (two-tailed, 10000 permutatinos). It produces the subsequent plot.

```
sample1 = randn(1,5000);
sample2 = randn(1,5000) + .1;
permutationTest(sample1, sample2, 10000, 'plot', 1);
```

![Output plot](./permutationTest.png)
