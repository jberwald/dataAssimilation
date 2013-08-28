function [ z_forecast ] = sea_ice_EW09_params( z_background, delta_t, t0, params )
% [final_year, full_time_series]=sea_ice_model_EW09_params( arguments )
%
% Model of vertical sea ice thermodynamics with no snow and no shortwave
% penetration, with quasi-stationary approximation for internal ice
% temperature (linear temperature profile). When there is no ice
% the ocean mixed layer thermally evolves.
% Prognostic variable is E=-Li*h+c*H*Tml (surface enthalpy).
%   Output: final_year or full_time_series = [t E Tsrf -E/Li E/(cml*Hml) Ftop F0 FT Fsw]
%
% Arguments:
%    
%    z_background : background state
%
%    delta_t : length of time to evolve model (in years)
%
%    t0 : start time, for restarting after assimilation step
%
%    params : hash table containing all parameters that the parameter
%    estimation technique is modifying. These will be updated below.
%    Eg., params.keys = {'kD', 'F0',...}
%         params.values = { 5, 6, ... }	  
%
% Reference:
%   Ian Eisenman and J.S. Wettlaufer, 2009. Nonlinear threshold behavior
%   during the loss of Arctic sea ice. Proceedings of the National Academy
%   of Sciences USA 106, 28-32.
%
% Ian Eisenman (ian@gps.caltech.edu), 2009
%
% Modifications: Jesse Berwald, Tom Bellsky, Lewis Mitchell, August
% 2013 (contact: jberwald@gmail.com)
% 
% 1. Removed plotting
% 2. modified to accept background state 'z_background' and variable length
% forecast length 'delta_t'

global Fbot ai ao ki Li Hml cml F0input Fswinput FTinput tinput v0 tanha Tlin Tsmelt intmeth

% === Parameters (these are defaults!) === 
% Measuring time in years for d/dt, using W/m^2 for fluxes and then
% multiplying d/dt by number of seconds in year.
yr=3.16e7; % seconds in year
Fbot=2;
ai=0.68;
ao=0.2;
ki=2;
Li=3*10^8/yr; % kJ-->W : W = 1000kJ/s
Hml=50;
cml=4*10^6/yr;
Tsmelt=0; % sfc temperature for onset of melt
Tlin=0; % linearity of ice surface temperature: 0 for sea ice model, 1 for as mixed layer (linear)
v0=0.1; % ice export
dF=0;   % imposed surface heat flux (annual)


% Reset parameters from 'params' here
for k=params.keys
    thisKey = char( k ); % matlab stores keys as cell objects, get at pure string
    eval( [ thisKey, '=', int2str( params( thisKey ) ) ';'] );
end

tanha=0.5*Li; % if not equal to zero, gives width of albedo dependence


% = For computing Ftop =
% atmmod: 0 for surface fluxes specified, 1 for active atmosphere with KLW
% from cloud fraction, 2 for specified F0(tinput) and FT(tinput).
atmmod=1;
tinput=(0.5:11.5)/12; % times that forcing is input at
% Stefan-Boltzmann linearization
sigma0=316;
sigmaT=3.9;
% surface fluxes from Maykut & Untersteiner (1971)
Fswinput=15.9*[0 0 1.9 9.9 17.7 19.2 13.6 9.0 3.7 0.4 0 0];
SH=15.9*[1.18 0.76 0.72 0.29 -0.45 -0.39 -0.30 -0.40 -0.17 0.10 0.56 0.79];
LH=15.9*[0 -0.02 -0.03 -0.09 -0.46 -0.70 -0.64 -0.66 -0.39 -0.19 -0.01 -0.01];
LWd=15.9*[10.4 10.3 10.3 11.6 15.1 18.0 19.1 18.7 16.5 13.9 11.2 10.9];
% atmospheric model parameters
% Cloud fraction from Maykut & Church (1973)
Cl=[4.9 4.7 4.2 5.9 8.1 8.0 8.3 8.8 8.9 8.4 7.1 5.2]/10;
% 0-70N mean temp based on NCEP-NCAR reanalysis
Tsouth=18*(1-1/3*cos((tinput-20/365)*2*pi));
kD=2.9;   % 2.7 in the article (LM)
tau0=0.5; % optical thickness independent of cloud fraction
tauc=3.6; % optical thickness dependence on cloud fraction

%JJB: max_dur == delta_t
max_dur = delta_t; %50; % max duration of integration 

% JJB: each integration typically begins at new bakground state
if z_background == 0
   E0 =  -3.1*Li; % E=-Li*hi+cml*Hml*Tml
else
  E0 = z_background; % initial value of E=-Li*hi+cml*Hml*Tml
end

% JJB: Is this an intermediate time step?


intmeth='linear'; % method for interpolating input data in time
silent=0; % if 1, don't display error or final value
RelTol=1e-7; % for ODE solver
AbsTol=1e-6; % for ODE solver

odesolver='ode45'; % which ode solver to use

% === EDIT PARAMETERS BELOW HERE ===
% parameter value changes as input to function (batch mode)
%if nargin>0, for j=1:length(varargin), eval([varargin{j} ';']), end, end
%if nargout==0
    % = Enter changes to default parameter values here for interactive run =
    %dF=22;
%end

% = atmospheric fluxes =
if atmmod==0 % specified surface fluxes
    F0input=-LH-SH-LWd+sigma0;
    FTinput=sigmaT+zeros(1,12);
elseif atmmod==1 % active equilibrium atmospheric model with KLW from cloud fraction, v0>0
    KLW=1./(tau0+tauc*Cl);
    F0input=KLW*sigma0-dF-kD*Tsouth/2;
    FTinput=KLW*sigmaT+kD/2;
elseif atmmod==2 % F0,FT specified as input
    F0input=F0;
    FTinput=FT;
end

% ======

% === integration ===
% run integration one yr at time until spun_up condition or max
% duration
% JJB: integrate from t0 to t0 + delta_t

%intyr=1; 
E=[]; t=[];
options=odeset('RelTol',RelTol,'AbsTol',AbsTol);

%while intyr<=max_dur

if strcmp(odesolver,'ode45'), [tt,EE]=ode45(@dEdt, t0 + [0 max_dur ], E0,options); end
if strcmp(odesolver,'ode23t'), [tt,EE]=ode23t(@dEdt, t0 + [0 max_dur ], E0,options); end
if strcmp(odesolver,'ode15s'), [tt,EE]=ode15s(@dEdt, t0 + [0 max_dur ],E0,options); end
E=[E; EE]; t=[t; tt];

%intyr2=intyr;
E0=EE(end,:); %intyr=intyr+1;
%if ( abs(EE(end,:)-EE(1,:)) < spun_up ) && spin_up, intyr=max_dur+1; end

%end
%if silent==0 % display final value/error in command window
 %   if intyr2<max_dur
 %       disp(['Integration ended at year ' num2str(intyr2) ' with final value E0=' mat2str(EE(end),2)])
 %   else
 %       disp(['Integration ended at year ' num2str(intyr2) ' with final value E0=' mat2str(EE(end),2)])
 %       disp(['ERROR: Integral did not converge during ' num2str(intyr2) ' years.'])
 %   end
%end

% === output ===
% find equilibrium surface temperatures, etc from t,E data
[F Tsrf Ftop F0 FT Fsw] = dEdt(tt,EE);
Y=[tt EE Tsrf -EE/Li EE/(cml*Hml) Ftop F0 FT Fsw];
final_year=Y;

% E is the concatenation of all EE's during the max_dur integration
[F Tsrf Ftop F0 FT Fsw]=dEdt(t,E);

% update time vector to append to 
%t = t0 + t;
Y=[t E Tsrf -E/Li E/(cml*Hml) Ftop F0 FT Fsw];
z_forecast = Y;

% ===================


% === model equations ===
function [F Tsrf Ftop F0 FT Fsw]=dEdt(t,E)
global Fbot cml Hml ki Li F0input Fswinput FTinput tinput ai ao intmeth v0 tanha Tlin Tsmelt intmeth
F0=interp1([tinput(end)-1 tinput tinput(1)+1],F0input([end 1:end 1]),mod(t,1),intmeth);
Fsw=interp1([tinput(end)-1 tinput tinput(1)+1],Fswinput([end 1:end 1]),mod(t,1),intmeth);
FT=interp1([tinput(end)-1 tinput tinput(1)+1],FTinput([end 1:end 1]),mod(t,1),intmeth);

Tsrf0=-(E.*F0-E*(1-ai).*Fsw)./(E.*FT-ki*Li); % ice surface temperature when no surface melt
Tsrfi=Tsrf0.*(Tsrf0<Tsmelt);
Tsrf=Tsrfi.*(E<0)+E/(cml*Hml).*(E>=0);

if Tlin==1 % for sea ice as ocean mixed layer with only albedo changing when Ts<0
    Tsrf=E/(cml*Hml).*(E<0)+E/(cml*Hml).*(E>=0);
end

a = tanh_albedo( ai, ao, E );
% a = linear_albedo( ai, ao, E );
% original_tanh( ai, ao, E, tanha );

F=-F0+(1-a).*Fsw-FT.*Tsrf+Fbot-v0*E.*(E<0);
Ftop=F0+a.*Fsw+FT.*Tsrf; % output for diagnostic plots
