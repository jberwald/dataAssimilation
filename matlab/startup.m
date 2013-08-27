addpath(genpath('matlab'))

set(0,'DefaultLineLinewidth',2)
%set(0,'DefaultTextInterpreter','latex')
set(0,'defaultaxesfontname','helvetica')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)

%% Don't fuck me with the random numbers, motherfucker
randn('state',sum(100*clock))
rand('state',sum(100*clock))

%% Fuck you, MATLAB default colours!
set(0,'DefaultAxesColorOrder',cbrewer('qual','Set1',8))
set(0,'DefaultFigureColormap',cbrewer('seq','YlOrRd',64)); 
%%%%% Note: Permuted the first 3 colors of Set1 using: %%%%%
% colorbrewer.qual.Set1{8}(1:3,:) = ...
%  circshift(colorbrewer.qual.Set1{8}(1:3,:),-1);

