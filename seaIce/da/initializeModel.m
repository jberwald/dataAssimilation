

fprintf('Initializing model\n')

D = 1;                                      % number of model variables

% z0 = rand(D,1);                           % non spun-up init. condition
z0 = 0;                                     % non spun-up init. condition
% tSpinup = 1;                                % spinup time length (years)
% zSpinup = run_sea_ice(z0,tSpinup,-tSpinup); % spin that mofo up!
% z0 = zSpinup(end,2);