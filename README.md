# Connectivity Change Point Detection

Connectivity Change Point Detection (CCPD) is a fully automated two-stage approach which pools information across multiple subjects to estimate change points in functional connectivity and subsequently estimates the brain networks within each
state phase lying between consecutive change points. The number and positioning of the change points are unknown and learned from the data in the first stage, by modeling a time-dependent connectivity metric under a fused lasso approach. In the second stage, the brain functional network for each state phase is inferred via sparse inverse covariance matrices. The step could be done through the Matlab function ‘CCPD’.

We have also extended the proposed CCPD approach to infer subject level change points without assuming  similar change points locations across subjects. The change points for each subject and corresponding graphs can be estimated using Matlab function ‘CCPD_single’.

## Usage

```
[change_points,graph] = CCPD(Y,count)

%Input:
%    Y: 	Data set, format as T*V*Nsub
%     where
%     T is number of time points
%     V is number of ROIs
%     Nsub is number of subjects
%    count: number of subsampling (default value 100)
%
%Output:
%     change_points: change points for all subjects
%     graph: graph of each bin, format is P*P*Nbin.
%        The max value of Nbin is 20
```

```
[change_points,graph] = CCPD_single(Y)
%input:
%     Y: Data set, format is T*V*Nsub
%      T is number of time points
%      V is number of ROIs
%      Nsub is number of subjects
%
%Output:
%     change_points: change points for each subject. %        format is Nsub*20, where 20 is the max %        number of possible change points for
%        each subject.
%     graph: graph of each bin for each subject, %        format is P*P*Nbin*Nsub.
%        The max value of Nbin is 20
```


## Requirements

QUIC: http://www.cs.utexas.edu/~sustik/QUIC/

## Reference

Kundu, S., Ming, J., Pierce, J., McDowell, J. and Guo, Y., Estimating Dynamic Brain Functional Networks Using Multi-subject fMRI Data.

## Contact

Please send comments and bug reports to: suprateek.kundu@emory.edu /jming2@emory.edu
