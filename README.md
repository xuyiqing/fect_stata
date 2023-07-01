# fect for STATA

## A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data

---

**Authors:** Licheng Liu (MIT); Ye Wang (NYU); Yiqing Xu (Stanford); Ziyi Liu (Uchicago); Shijian Liu (Duke)

**Maintainer:** Yiqing Xu [<yiqingxu@stanford.edu>]  

Installation
=======

As pre-requisites, the `reghdfe` package and `_gwtmean` package need to be installed. 

```Stata
. ssc install reghdfe, replace
. ssc install _gwtmean, replace
```

To install the `fect` package with Stata v14 or greater

```Stata
. cap ado uninstall fect //in-case already installed
. net install fect, from(https://raw.githubusercontent.com/xuyiqing/fect_stata/master/) replace
```

If the above link is unaccessible for your device, use the "Download Zip" button above, unzip it to a local directory, and then replace the above `net install` with

```Stata
. net install fect, from(full_local_path_to_files) replace
```

How to use
=======

 [Examples](https://yiqingxu.org/packages/fect/stata/fect_md.html)
 
 ---

**Reference:** "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data" Available at SSRN: https://papers.ssrn.com/abstract=3555463

---

