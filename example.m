%Example of CCPD toolbox
%Provided by Suprateek Kundu and Jin Ming @ Emory University

%QUIC is needed for our CCPD, you could install QUIC by following the
%README file in QUIC folder
%mex -lmwlapack QUIC.cpp QUIC-mex.cpp -output QUIC.mexmaci64

%Multi-subject case
load('sample_data.mat')
%There are 60 subjects of 200 time points with 10 ROIs, there are 3 chagne points [40, 60, 175]
[change_points,graph]=BG_final(Y,count);

%Single subject case
[change_points,graph]=CCPD_single(Y);