

% dInfMult = str2num(getenv('dInf'))/100
dInfMult = 1
dInfAdd = 0;

N_ens = 50;

dt_f = 0.05; % forecast length




%% observation operator
N_obs = I/3;

H = zeros(N_obs,I);
for kk = 1:size(H,1)
  H(kk,(kk-1)*floor(I/N_obs)+1) = 1;  % obs every N_obs, equally spaced
end

obLocs = 1:round(I/N_obs):I;


R = 0.2*eye(N_obs);

%% Gaspari-Cohn localisation

% windowWidth = str2num(getenv('rLoc'))
windowWidth = I;



%%%%% no localization
rho = ones(I);

for yy = 1:I
 for xx = 1:I
   dd = min([abs(xx-yy) abs(xx-yy-I) abs(xx-yy+I)]);
%     rr = dd/(sqrt(10/3)*1.8);
   rr = dd/windowWidth;
   if(rr<=1)
     rho(yy,xx) = 1-0.25*rr^5+0.5*rr^4+(5/8)*rr^3-(5/3)*rr^2;
   elseif(rr>1 && rr <= 2)
     rho(yy,xx) = (1/12)*rr^5-0.5*rr^4+(5/8)*rr^3+(5/3)*rr^2 ...
       -5*rr+4-2/3/rr;
   else
     rho(yy,xx) = 0;
   end
 end
end

