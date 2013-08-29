% defaultParams.m
% the default parameters for the Eisenman sea ice model.

dF=0;   % imposed surface heat flux (annual warming rate)
Fbot=2;
ai=0.68;
ao=0.2;
ki=2;
Hml=50;
Tlin=0; % linearity of ice surface temperature: 0 for sea ice model, 1 for as mixed layer (linear)
v0=0.1; % ice export

keys = {'dF','Fbot','ai','ao','ki','Hml','Tlin','v0'};
values = {10,2,0.68,0.2,2,50,0,0.1};


defaultParams = containers.Map(keys,values);
