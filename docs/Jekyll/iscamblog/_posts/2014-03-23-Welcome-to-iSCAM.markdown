---
layout: post
title:  "Welcome to iSCAM!"
date:   2014-03-23 06:50:11
categories: iscam overview
---

# integrated Statistical Catch Age Model (iSCAM)
Initially developed by Steve Martell at the University of British Columbia and has now evolved to a more substantial project that is used and maintained by a broader community.

The source code for the project is maintained at a [Github](https://github.com/smartell/iSCAM) repository.

# Brief History
This project originated in 2010 as a contract with DFO to develop a new statistical catch-age model for the assessment of Pacific herring stocks in British Columbia.  The approach was to develop a generic model that could be adapted to accommodate different data and assumptions from each of five separate herring stocks.

The original model consisted of a single stock with one sex, could be fit to multiple abundance indexes and age-composition data.  Features that were unique to iSCAM that were not found in other modeling platforms include:

- Estimation of total variance and variance partitioning.
- Self-weighting of age-composition data.
- Parametric and non-parametric selectivity functions.
- Time-varying selectivity.
- 2-deminsional bicubic splines.
- Time varying natural mortality using cubic spline interpolation.
- Explicit modeling of roe-based fisheries.
- MSY-based reference points for multiple fleets with different selectivities.
- Command line options for retrospective analysis & simulation testing.

Over time a library of statistical functions (stats.cxx) was developed that eventually ended up being directly incorporated into [AD Model Builder](http://admb-project.org) and the stats.cxx library was deprecated from iSCAM.



# iSCAM moves to the IPHC
In September 2012, I took a position with the International Pacific Halibut Commission to begin working on harvest policy issues for Pacific halibut.  To this end, iSCAM has evolved considerably since this time to incorporate sex-specific parameters, spatially explicitly polygons representing management areas, and additional flexibility for investigating sub-stock structure, cumulative effects of size-selective fishing, or even multi-species assessments.


