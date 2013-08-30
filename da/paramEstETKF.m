% function [] = twinExperimentETKF
% twinExperimentETKF.m : Assimilate some real data
% LM,JBB,TB 8/26/13 ASU
% Input: No input
% Output: A whole DA twin experiment


%% Initialisation
clear all
randn('state',100*sum(clock));  % so we don't get fucked by MATLAB's rng.
addpath(genpath('../'));



% sim_no = str2num(getenv('PASS')); % for farmed-out parallel jobs
sim_no = 1;

fprintf('Simulation: %d\n',sim_no)

nState = 1;         % dimension of model state vector
nParam = 1;         % number of parameters to estimate
D = nState+nParam;  % total variables (size of model + parameter space)
z0 = 0;             % default initial condition

%% ETKF ensemble size and forecast length

N_ens = 10;   % number of ensemble members
dt_f = 1/12;  % forecast length (years)


%% Inflation and Localization

%%%%% multiplicative and additive inflation factors %%%%%
% dInfMult = str2num(getenv('dInf'))/100  % for farmed-out shit
dInfMult = 1.0;
dInfAdd = 0.0;

%%%%% no localization necessary...
rho = ones(D);

% %% load SIE data
% 
% SIEData = load('../seaIce/data/sea_ice_extent_1979_2009.txt');
% 
% t = load('../seaIce/data/SIE_dates.txt');
% 
% N_obs = D;            % number of observations
% % H = eye(N_obs);       % simple observation operator
% H = -1;               % linearized observation operator
% R = 0.05*eye(N_obs);  % specified ob error covariance matrix
% 
% y = SIEData(:,4)';

%% Make truth and observations (toy twin experiment)

loadDefaultParams;  % load default parameters into struct 'defaultParams'

%%%%% set the true parameters we wish to estimate %%%%%
% defaultParams('dF') = 50; % imposed surface heat flux
defaultParams('Fbot') = 10; % imposed bottom heat flux

%%%%% spin up for 1 year %%%%%
zInit = run_sea_ice_params(z0,1,-1,defaultParams);
z0 = zInit(end,2);


Tf = 10;  % length of simulation (years)
t = 0:dt_f:Tf;

fprintf('Creating truth\n')
zt(:,1) = z0;
for n = 2:length(t)
%   thisForecast = run_sea_ice(zt(:,n-1),dt_f,t(n-1));
%     defaultParams('hAlpha') = n/length(t);
  thisForecast = run_sea_ice_params(zt(:,n-1),dt_f,t(n-1), defaultParams );
  zt(n) = thisForecast(end,2);
end



N_obs = nState;            % number of observations
% H = eye(N_obs);          % simple observation operator
H = [-1 zeros(1,nParam)];  % linearized observation operator
R = 0.1*eye(N_obs);        % specified ob error covariance matrix

% y = H*zt + sqrtm(R)*randn(size(zt(1:N_obs,:)));
y = enthalpyToSIE(zt) + sqrtm(R)*randn(size(zt(1:N_obs,:)));




%% Make initial ensemble with guess at F


% %%%%% real observations version %%%%%
% % zb = SIEToEnthalpy(y(1) + R*randn(D,N_ens));  % perturb observation
% zb = SIEToEnthalpy(y(1)) + 3*randn(D,N_ens);  % perturb observation
% zb0 = zb;                                     % store initial zb


%%%%% toy experiment version %%%%%
zb(1:nState,:) = ...
  zt(:,1)*ones(1,N_ens)+randn(nState,N_ens); % perturb true inital state
% zb(nState+1:D,:) = randn(nParam,N_ens);      % guess at parameter
zb(nState+1:D,:) = rand(nParam,N_ens);      % guess at parameter (uniform)
zb0 = zb;   

%% ETKF

zf_m = zeros(D,length(t));
za_m = zeros(D,length(t));

za_mean = mean(zb,2);
za_m(:,1) = za_mean;

fprintf('Doin'' the ETKF\n')
tic


modelValues = {0,0.5,2,0.68,0.2,2,50,0,0.1};

modelParams = containers.Map(keys,modelValues);


% for n_a = 2:length(t)
for n_a = 2:length(t)

  if(~mod(round(n_a/length(t)*100),10))
    tElapsed = toc;
    fprintf('%d%% complete, %ds elapsed\n',...
      round(n_a/length(t)*100),round(tElapsed)); 
  end

  %%%%% Forecast Step %%%%%
  zf = zeros(D,N_ens);
  for k = 1:N_ens
%     modelParams('dF') = zb(nState+1,k);
%     modelParams('hAlpha') = zb(nState+1,k);
    modelParams('Fbot') = zb(nState+1,k);
    this_zf = run_sea_ice_params(zb(1:nState,k),dt_f,t(n_a-1),modelParams);
    zf(1:nState,k) = this_zf(end,2);
    zf(nState+1:D,k) = zb(nState+1:D,k);
  end
  
  Pf(:,:,n_a) = cov(zf');
  zf_m(:,n_a) = mean(zf,2);
  
  %%%%% adaptive inflation (Wang & Bishop 2003) %%%%%
  dOMB = (y(:,n_a) - enthalpyToSIE(zf_m(:,n_a)));
  dInfMult = max((dOMB'*dOMB - trace(R))/trace(H*full(Pf(:,:,n_a))*H'),1);
  adInf(n_a) = dInfMult;
  
%   %%%%% give obs same weight as forecast %%%%%
%   dInfMult = trace(R)/trace(H*Pf(n_a)*H');
%   adInf(n_a) = max(dInfMult,1);

  
  [za,Pa(:,:,n_a)] = ...
    ETKF_onestep_obsOperator(zf,y(:,n_a),...
    @enthalpyToSIE,H,R,dInfMult,dInfAdd,rho);
  
  za_m(:,n_a) = mean(za,2);
  
  zb = za;
  
end
toc



% startPointRMSE = round(size(za_m,2)*0.2);
% fprintf('RMSE: %.2f\n',RMSE(za_m(:,startPointRMSE:end) - zt_o(1:I,startPointRMSE:end)));

%   save(['long/N',num2str(N_ens),'_l',num2str(locLength),'_d',num2str(dInfMult),'_',num2str(sim_no),'.mat'])


%% Pretty pictures
figure(1)
clf(1)
% plot(t,zt(1,:),t(2:end),za_m(1:end-1),t,y,'x')
% plot(t,zt(1,:),t,za_m,t,y,'x')
plot(t,zt(1,:),t,za_m(1,:),t,SIEToEnthalpy(y),'x',t(2:end),zf_m(1,2:end))
xlabel('Time (years)')
ylabel('Energy')

figure(2)
clf(2)
% plot(t,zt(1,:),t(2:end),za_m(1:end-1),t,y,'x')
% plot(t,zt(1,:),t,za_m,t,y,'x')
% plot(t,defaultParams('dF')*ones(1,length(t)),'k--',t,za_m(2,:))
% plot(t,defaultParams('hAlpha')*ones(1,length(t)),'k--',t,za_m(2,:))
plot(t,defaultParams('Fbot')*ones(1,length(t)),'k--',t,za_m(2,:))
% plot(t,defaultParams('dF')*ones(1,length(t)),'k--',...
%   t,defaultParams('Fbot')*ones(1,length(t)),'k-.',...
%   t,za_m(2:end,:))
% plot(t,defaultParams('dF')*ones(1,length(t)),'k--',...
%   t,defaultParams('Fbot')*ones(1,length(t)),'k-.')
% hold on
% errorbar(t,za_m(2,:),squeeze(Pa(2,2,:)))
% errorbar(t,za_m(3,:),squeeze(Pa(3,3,:)),'r')
% hold off
xlabel('Time (years)')
ylabel('Parameter')



