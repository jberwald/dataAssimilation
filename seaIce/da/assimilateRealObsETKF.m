% function [] = twinExperimentETKF
% twinExperimentETKF.m : Do a made-up sea ice experiment
% LM,JBB,TB 8/26/13 ASU
% Input: No input
% Output: A whole DA twin experiment


%% Initialisation
clear all
randn('state',100*sum(clock));  % so we don't get fucked by MATLAB's rng.
% addpath('../JJBCode');



% sim_no = str2num(getenv('PASS')); % for farmed-out parallel jobs
sim_no = 1;

fprintf('Simulation: %d\n',sim_no)


initializeModel;  % spinup, eliminate transient behaviour


%% ETKF ensemble size and forecast length

N_ens = 10;   % number of ensemble members
% dt_f = 1/12;  % forecast length (years)


%% Inflation and Localization

%%%%% multiplicative and additive inflation factors %%%%%
% dInfMult = str2num(getenv('dInf'))/100  % for farmed-out shit
dInfMult = 1.5;
dInfAdd = 0.0;

%%%%% no localization necessary...
rho = ones(D);

%% load SIE data

SIEData = load('sea_ice_extent_1979_2009.txt');

t = load('SIE_dates.txt');

N_obs = D;            % number of observations
% H = eye(N_obs);       % simple observation operator
H = -1;               % linearized observation operator
R = 0.05*eye(N_obs);  % specified ob error covariance matrix

y = SIEData(:,4)';



%% Make initial ensemble with guess at F

% zb = SIEToEnthalpy(y(1) + R*randn(D,N_ens));  % perturb observation
zb = SIEToEnthalpy(y(1)) + 3*randn(D,N_ens);  % perturb observation
zb0 = zb;                                     % store initial zb

%% ETKF

zf_m = zeros(D,length(t));
za_m = zeros(D,length(t));

za_mean = mean(zb,2);
za_m(:,1) = za_mean;

fprintf('Doin'' the ETKF\n')
tic

nSkip = 30; % number of obs to skip each time
% for n_a = 2:length(t)
for n_a = nSkip+1:nSkip:length(t)

  if(mod(n_a,100)==1)
    fprintf('%d\n',n_a);
    toc
    tic
  end

  %%%%% Forecast Step %%%%%
  zf = zeros(D,N_ens);
  dt_f = (t(n_a)-t(n_a-nSkip))/365; % time between obs (in years)
  currentTime = (t(n_a-nSkip)-t(1))/365;
  ensFcst = [];
  for k = 1:N_ens
    this_zf = run_sea_ice(zb(:,k),dt_f,currentTime);
    zf(:,k) = this_zf(end,2);
    ensFcst{k} = this_zf(:,2);
  end
  
  Pf(n_a) = cov(zf');
  zf_m(:,n_a) = mean(zf,2);
  
%   %%%%% adaptive inflation (Wang & Bishop 2003) %%%%%
%   dOMB = (y(:,n_a) - enthalpyToSIE(zf_m(:,n_a)));
%   dInfMult = (dOMB'*dOMB - trace(R))/trace(H*Pf(n_a)*H');
%   adInf(n_a) = max(dInfMult,1);
  
  %%%%% give obs same weight as forecast %%%%%
  dInfMult = trace(R)/trace(H*Pf(n_a)*H');
  adInf(n_a) = max(dInfMult,1);

  
  [za,Pa(n_a)] = ...
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
% plot(t(2:end),zf_m(2:end),t(2:end),za_m(2:end),t,SIEToEnthalpy(y),'x')
plot(t(nSkip+1:nSkip:end),zf_m(nSkip+1:nSkip:end),...
  t(nSkip+1:nSkip:end),za_m(nSkip+1:nSkip:end),...
  t([1 nSkip+1:nSkip:length(t)]),...
  SIEToEnthalpy(y([1 nSkip+1:nSkip:length(y)])),'x')
xlabel('Time (years)')
ylabel('Energy')



figure(2)
clf(2)
% plot(t(2:200),log10(Pf(2:200)),t(2:200),log10(Pa(2:200)))
plot(t(nSkip+1:nSkip:length(Pf)),log10(Pf(nSkip+1:nSkip:end)),...
  t(nSkip+1:nSkip:length(Pa)),log10(Pa(nSkip+1:nSkip:end)),...
  t(nSkip+1:nSkip:length(Pa)),...
  log10(Pa(nSkip+1:nSkip:end).*adInf(nSkip+1:nSkip:end)))
% plot(t(nSkip+1:nSkip:length(Pf)),log10(Pf(nSkip+1:nSkip:end)),...
%   t(nSkip+1:nSkip:length(Pa)),log10(Pa(nSkip+1:nSkip:end)))
  
