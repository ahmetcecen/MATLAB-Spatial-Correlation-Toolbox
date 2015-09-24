MATLAB Spatial Correlation Toolbox
==========================
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.28859.svg)](http://dx.doi.org/10.5281/zenodo.28859)

Table of contents
=================

  * [Introduction](#introduction)
  * [Correlation and Convolution](#correlation-and-convolution)
  * [Use for 2-Point Statistics](#use-for-2-point-statistics)
    * [Setting Up a Problem](#setting-up-a-problem)
    * [Converting the Problem to a Computation](#converting-the-problem-to-a-computation)
    * [TwoPoint](#twopoint)

# Introduction

A toolbox designed specifically for computing spatial correlations of gigantic datasets, with support for regular sized datasets as well. The toolbox takes advantage of the memory mapping functionality in MATLAB to operate on a chunk of the data at a time. The overal strategy is ineffective for parrallelization as it involves tremendous overhead, but it is ideal for "sequentialization", when the algorithm needs to be able to run on a simple everyday machine, and it is okay for it to take a bit longer than the optimal calculation. I refer to it as the Patched Correlation Method since it uses patches of data at a time, although you are free to not call it that.

# Correlation and Convolution

`GG = CorrMaster('full','auto',cutoff,H1)` and `GG = CorrMaster('full','cross',cutoff,H1,H2)` will compute the plain old circular ('full' option) auto or cross correlation using FFTs. It will return the correlation results trimmed by a cutoff, with the final output being size 2*cutoff-1 in all dimensions, and the 0 index shifted to the center.

`GG = CorrMaster('patched','auto',cutoff,DataFile,winmulti)` and `GG = CorrMaster('patched','cross',cutoff,DataFile,winmulti,DataFile2)` will compute the patched auto or cross correlation. The "datafile" is the path to a matfile that contains the data you need to process, **already zero padded by exactly cutoff on BOTH sides** (you can use padarray for this or the memory mapped version below) and the name of the variable has to be H1. An example H1 for 100x100 image is with a cutoff of 20 is a variable H1 inside a data.m file that is 140x140 with 20 pixels zero padded on each side. 

There is a heuristic in place to calculate an optimal patch size for minimum memory requirement to divide up the data, and winmulti is the multiplier for that size. So a winmulti of 2 will use twice the optimal size. This is useful if you can afford the extra memory, since the higher the winmulti, the faster the computation time (as long as you don't go out of bounds into a pagefile).

There are two files to pad and element-wise multiply 2D or 3D data for use with this code called ewpmfile and padmfile, which should be straight forward to use since they share the same input terminology, with the exception of being able to handle any variable inside a matfile chosen by you, not only H1.

# Use for 2-Point Statistics

The following article describes in detail the theory necessary to properly utilize this toolbox to calculate 2-pt statistics: UNDER REVIEW  

Below will be a less detailed guide.

## Setting Up a Problem

The user is responsible with finding a suitable statistical question they want answered, with a clearly defined event and a clearly defined population. Depending on the question, the quantity they need to calculate and its pre-processing stages can change significantly. Only some of many considerations are exampled below:

If your problem requires **non-periodic** boundaries you need to pad your image on all dimensions with **0** values in the size of the longest vector of interest. Let us assume you have a microstructure volume represented in a **100x100x100** voxel grid, and you are interested in spatial statistics of vectors of size that fit inside a **20x20x20** voxel grid; then you would pad **20** voxels of **0** values on all 3 dimensions of your data, resulting in a **120x120x120** voxel microstructure.

If your problem requires **masking** you need to create a second image with the exact size of your actual input, with values representing mask values. This feature can be used when there is an uncertainty associated with each pixel, or when the image represents non-rectangular data to strike out pixels that don't belong in the data etc.. 

If you need to perform element-wise multiplication or padding of very large image files, use the memory mapped versions mentioned above.

## Converting the Problem to a Computation

Here is an example procedure that calculates a non-periodic memory mapped auto correlation without a mask:

I want to find the probability of finding every possible vector that is shorter than cutoff in this data, while discarding any vector that crosses the data boundaries.

```matlab
padmfile('yourdir\yourdata.m',myvar,cutoff,'mydata.m') <- Pad the image by cutoff. mydata will be saved with a variable named H1 for convenience.
GG = CorrMaster('patched','auto',cutoff,'mydata',2.23) <- Calculate auto-correlation.

padmfile('normalizer.m',myvar,cutoff,'mynorm.m') <- Since there is no mask, the normalization is a matrix of ones everywhere, then padded by cutoff.
BB = CorrMaster('patched','auto',cutoff,'mynorm',2.23) <- Calculate normalization.
 
AutoCorrelation=GG./BB <- Because I asked the question this way, this is the corresponding 2-pt statistics.
```

## TwoPoint

`TwoPoint` is a wrapper for a handful of common two point statistic calculation schemes. It calculates the two point statistics using the straightforward convolution method without patching. Use it if you have small datasets or an abundance of memory. 

`GG = TwoPoint('auto',cutoff,periodicity,H1)` will calculate the autocorrelation of H1.

`GG = TwoPoint('auto',cutoff,periodicity,H1,M1)` will calculate the autocorrelation of H1 with a mask M1.

`GG = TwoPoint('cross',cutoff,periodicity,H1,H2)` will calculate the crosscorrelation of H1 with H2.

`GG = TwoPoint('cross',cutoff,periodicity,H1,H2,M1)` will calculate the crosscorrelation of H1 with H2, with a uniform mask M1.

`GG = TwoPoint('cross',cutoff,periodicity,H1,H2,M1,M2)` will calculate the crosscorrelation of H1 with mask M1 and H2 with mask M2.


