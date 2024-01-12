# Block-Majorization-Minimization-with-Extrapolation
Block Majorization Minimization with Extrapolation and Application to $\beta$-NMF

This MATLAB software reproduces the results from the following paper:

```
@misc{LeThiKhanhHien_BMMex2023,
      title={Block Majorization Minimization with Extrapolation and Application to $\beta$-NMF}, 
      author={Le Thi Khanh Hien and Valentin Leplat and Nicolas Gillis},
      year={2023},
      eprint={},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```

## Acknowledgements

The baseline algorithms used in the manuscript are courtesy of their respective authors.


## Content
 
 - /Libraries : contains helpful libaries; in particular Libraries/group_robust_NMF/ contains the code of group robust NMF proposed in Ref[1] and Algorithm 2 from Ref[2] aimed at solving min-vol KL NMF.
 
 - /Datasets : contains test data sets.

 - /Utils : contains helpful files and MatLab routines to run the demos.
   
 - /Methods: contains the MatLab implementations of the Algorithms developped in the paper. 

 - test files detailed in next Section
   
## Test files
 
 Test files are available. To proceed, open and start one of the following files:
 
- test_Example_in_Intro_CBCL.m : run demo for generating Figure 1 from Section 1.3 of the paper.
- test_Hyperspectral_Images_betaNMF32.m : run demo for $\beta=3/2$-NMF for hyperspectral imaging, see Section 5.1 of the paper. 
- test_Images_DocClassification.m : run benchmark for KL-NMF for Images data sets and document classification (data sets specified within the file), see section 5.1 of th paper.
- test_Audio_minvolKLNMF.m : run demo for min-vol KL NMF for blind audio SS tasks, section 5.2 of the paper.
 
## Tuning the parameters
 
 There are several parameters that you can choose:
 - $r$: the factorization rank
 - $\tilde{\lambda}$: the relative weight parameter for min-vol penalty, see Section XX for more details
 - the maximum computation time for benchmarking
 
For benchmarked approaches, the parameters have been tuned according to the original works:
- for audio SS tasks; we followed Ref [3].
- for cbcl and document classification; we followed Ref [4].
 
## References:
[1]: C. Févotte and N. Dobigeon, "Nonlinear Hyperspectral Unmixing With Robust Nonnegative Matrix Factorization," in IEEE Transactions on Image Processing, vol. 24, no. 12, pp. 4810-4819, Dec. 2015, doi: 10.1109/TIP.2015.2468177.

[2]: V. Leplat,  N. Gillis and J. Idier, "Multiplicative Updates for NMF with β-Divergences under Disjoint Equality Constraints", SIAM J. on Matrix Analysis and Applications 42 (2), 730-752, 2021. 

[3]: V. Leplat, N. Gillis and A.M.S. Ang, "Blind Audio Source Separation with Minimum-Volume Beta-Divergence NMF", IEEE Transactions on Signal Processing 68, pp. 3400-3410, 2020.  

[4]: L.T.K. Hien and N. Gillis. "Algorithms for nonnegative matrix factorization with the Kullback-Leibler divergence." Journal of Scientific Computing, (87):93, 2021.
