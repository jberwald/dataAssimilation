

%%%%% ETKF_onestep.m : Do one step of ETKF (analysis step only!)

%%%%% LM 1/15/13 UVM

function [za,Pa] = ETKF_onestep(zf,y,H,Hlin,R,dInfMult,dInfAdd,rhoLoc)
% ETKF_onestep - perform one step of ETKF
%
% Input arguments:
%   zf:       forecast ensemble (DxN_ens)
%   y:        observation for this time
%   H:        observation operator (function)
%   Hlin:     linearised observation operator
%   R:        observation error covariance matrix
%   dInfMult: multiplicative inflation factor (>1, scalar)
%   dInfAdd:  additive inflation factor (scalar or vector)
%   rhoLoc:   Gaspari-Cohn localization




N_ens = size(zf,2);


zf_mean = mean(zf,2);
zf_pert = zf - zf_mean*ones(1,N_ens);

%%%%% additive inflation on ensemble %%%%%
zf(end,:) = zf(end,:) + dInfAdd*(rand(1,N_ens)-0.5);

Pf = cov(zf');

%%%%% additive inflation on Pf %%%%%
Pf = dInfMult*Pf;  % multiplicative inflation


%%%%% Localization %%%%%
Pf = rhoLoc.*Pf;


%%%%% Analysis Step %%%%%
K = Pf*Hlin'/(Hlin*Pf*Hlin' + R);

za_mean = zf_mean + K*(y - H(zf_mean));


Pa = (eye(length(za_mean)) - K*Hlin)*Pf;


%%%%% ET Step (Bishop 2001) %%%%%
%       U = zf_pert'*(H'*inv(R)*H + H_F'*inv(R_F)*H_F)*zf_pert/(N_ens - 1);
U = zf_pert'*(Hlin'*inv(R)*Hlin)*zf_pert/(N_ens - 1);

if(~isnan(U))
  
  [~,Gamma,C] = svd(U);
  S = C*inv(sqrt(eye(N_ens) + Gamma))*C';
  
  %           [C,Gamma] = eig(U);
  %             S = C/(sqrtm(eye(N_ens) + Gamma));
  %           [~,zero_eig] = min(diag(Gamma));
  %           C(:,zero_eig) = [];
  %           Gamma(:,zero_eig) = [];
  %           Gamma(zero_eig,:) = [];
  %           S = (C/(sqrtm(eye(N_ens-1) + Gamma)))*C';
  
  za_pert = zf_pert*S;
  
  %%%%% Update zb = za step %%%%%
  za = za_mean*ones(1,N_ens) + za_pert;
  
  
end

