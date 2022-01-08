clear all; close all;
%% Vary CR

%% Constants
n = 1e3;
vec = ones(1,n);
vec2 = ones(1,2*n);

% Temperatures
THin = vec2*(240 + 273);     % [K] Given
TCout = vec2*(110 + 273);    % [K] Given

% Heat and mass flux
Qdot = 50*10^6;      % [J/s] Given
mdotC = 150;         % [kg/s] Given
cpH = 4450;          % [J/K] Not actually constant
cpC = 4250;          % [J/K] Constant

Cc = vec2*mdotC*cpC;

%% Questions
TCin = TCout - Qdot/(mdotC*cpC); % [K] Calculate cold-side inlet temperature

CR = linspace(0,1,n); % [-] Choose a value for heat capacity ratio CR NIET HOGER DAN 1
CR = [CR,flip(CR)];
hi = zeros(1,2*n); % When is Ch higher than Cc
hi(n+1:end) = 1;

Ch = zeros(1,2*n);
Ch(hi==1) = Cc(hi==1)./CR(hi==1); % When Ch > Cc
Ch(hi==0) = Cc(hi==0).*CR(hi==0); % When Ch < Cc

Cmin = min(Cc,Ch);
Cmax = max(Cc,Ch);

mdotH = Ch/cpH;
THout = THin - Cc./Ch.*(TCout-TCin);
THout(THout < TCin) = nan;

% For counterflow:
dT1 = THout - TCin;
dT2 = THin - TCout;
dT_lm = (dT2 - dT1)./log(dT2./dT1);
    
%%
UA = Qdot./dT_lm;
NTU = UA./Cmin;
eps = (1-exp(-NTU.*(1-CR)))./(1-CR.*exp(-NTU.*(1-CR)));


eps2 = zeros(1,2*n);
eps2(hi==0) = (THin(hi==0)-THout(hi==0))./(THin(hi==0)-TCin(hi==0));
eps2(hi==1) = (TCout(hi==1)-TCin(hi==1))./(THin(hi==1)-TCin(hi==1));
eps2(eps2>1) = nan;

%% Plot
figure()
plot(CR(hi==0), eps(hi==0)); hold on
plot(CR(hi==1), eps(hi==1))
title('eps')
legend('Ch < Cc', 'Ch > Cc')
grid on
xlabel('CR')

figure()
plot(CR(hi==0), eps2(hi==0)); hold on
plot(CR(hi==1), eps2(hi==1))
title('eps')
legend('Ch < Cc', 'Ch > Cc')
grid on
xlabel('CR')

figure()
plot(CR(hi==0), NTU(hi==0)); hold on
plot(CR(hi==1), NTU(hi==1))
title('NTU')
legend('Ch < Cc', 'Ch > Cc')
grid on
xlabel('CR')


% figure(3)
% plot(CR(hi==0), UA(hi==0)); hold on
% plot(CR(hi==1), UA(hi==1))
% title('UA')
% legend('Ch < Cc', 'Ch > Cc')
% grid on
% xlabel('CR')
% 
% figure(4)
% plot(CR(hi==0), mdotH(hi==0)); hold on
% plot(CR(hi==1), mdotH(hi==1))
% title('mdotH')
% legend('Ch < Cc', 'Ch > Cc')
% grid on
% xlabel('CR')
% 
% figure(5)
% plot(CR(hi==0), Ch(hi==0)); hold on
% plot(CR(hi==1), Ch(hi==1))
% title('Ch')
% legend('Ch < Cc', 'Ch > Cc')
% grid on
% xlabel('CR')

figure()
plot(NTU(hi==0),eps2(hi==0)); hold on
plot(NTU(hi==1),eps2(hi==1))
title('eps-NTU')
legend('Ch < Cc', 'Ch > Cc')
grid on
xlabel('NTU')
ylabel('eps2')
text(NTU(n-1),eps2(n-1),"Hier is CR = 1")
text(NTU(n+2),eps2(n+2),"Hier is CR = 1")
