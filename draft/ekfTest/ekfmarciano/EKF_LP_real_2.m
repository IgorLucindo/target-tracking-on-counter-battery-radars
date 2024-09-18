
close all
clear 
clc
%============ Morteiro 120mm - TIROs 2 e 3 28/01/22 ======================================
 
%medidas81T2_280122=readmatrix('R.CampanhaMorteiro81_Tiro 2_280122.xlsx'); % DADOS ADOUR Tiro 2  
%medidas81T2_280122= readmatrix('R.Campanha_Morteiro120_Tiro2_221222.xlsx');
medidas81T2_280122=readmatrix('R.CampanhaMorteiro81_Tiro 2_280122.xlsx'); % DADOS ADOUR  Tiro 3
medidas81T2_SISROT= medidas81T2_280122(1:120,:); % DADOS SISROT
medidas81T2_ADOUR=  medidas81T2_280122(121:end,:); % DADOS ADOUR 



%medidas81T3_280122=readmatrix('R.Campanha_Morteiro120_Tiro2_221222.xlsx'); % DADOS ADOUR  Tiro 3
% Tiro 2
% ========================≠≠≠=============================================================
%top1 = Top: 14:47:36.543342(disparo)  to1 = 14:47:41:353109 (início rastreio)
latenciaS = medidas81T2_SISROT(:,4);
latenciaA = medidas81T2_ADOUR(:,4);
VelDist = medidas81T2_ADOUR(:,13);
T_top2 = datenum( '16:56:54.303807', 'HH:MM:SS.FFF' ) .* (24*60*60) - ...
datenum( '00:00:00.000', 'HH:MM:SS.FFF' ) .* (24*60*60); % disparo
T_2 = datenum( '16:57:00.904133', 'HH:MM:SS.FFF' ) .* (24*60*60) - ...
datenum( '00:00:00.000', 'HH:MM:SS.FFF' ) .* (24*60*60); % início do rastreio 
%
T_2S = datenum( '16:56:54.908709', 'HH:MM:SS.FFF' ) .* (24*60*60) - ...
datenum( '00:00:00.000', 'HH:MM:SS.FFF' ) .* (24*60*60); % início do rastreio

% ADOUR
delta_to2 = T_2-T_top2;
tdados2=medidas81T2_ADOUR(:,3) -latenciaA; % instantes com estatos de "válidos"   = modo automático
t2=tdados2-tdados2(1)+delta_to2;% mudança de escala

% SISROT
delta_to2S = T_2S-T_top2;
tdados2S=medidas81T2_SISROT(:,3) - latenciaS;       % instantes com estatos de "válidos" = modo automático
t2S=tdados2S-tdados2S(1)+delta_to2S;       % mudança de escala
%==========================================================================
%Posição dos sensores [xEast,yNorth,zUp] = geodetic2enu(lat,lon,h,lat0,lon0,h0,spheroid)
wgs84 = wgs84Ellipsoid('meters');
[xlanc,ylanc,zlanc] = geodetic2enu(-23.049181,-43.633323,2.5,-23.049181,-43.633323,0,wgs84); % Lançadora
[xrad,yrad,zrad] = geodetic2enu(-23.049118,-43.632974,3.9,-23.049118,-43.632974,0,wgs84); % Radar
%==========================================================================
ele2A=medidas81T2_ADOUR(:,5);  % ELEVAÇÃO 
azi2A=medidas81T2_ADOUR(:,6);  % AZIMUTE
dist2A=medidas81T2_ADOUR(:,7); % DISTÂNCIA
SNR81 =medidas81T2_ADOUR(:,20);
v_rad =medidas81T2_ADOUR(:,13);
%Reamostragem Uniforme====================================================
%[XOUT,YOUT] = msresample(X,Y,N)
T = 50e-3;
Fs = 1/T;
[ele2Ar, t1r] = resample(ele2A,  t2,Fs); % medidas reamostradas
[azi2Ar, t2r] = resample(azi2A,  t2,Fs);
[dist2Ar, t3r] = resample(dist2A,t2,Fs);
[SNR81, t4r] = resample(SNR81,t2,Fs);

ultima_amostra = 511;% length(azi2A);%612;%length(azi2A);%;612
tr2A = t1r(1:ultima_amostra);
SNR81 = SNR81(1:ultima_amostra);
%tr2A = t2;
%==========================================================================
% Conversão de cordenadas Sistema Esférico Local AER para ENU
[xE2A,yN2A,zU2A] = aer2enu(azi2Ar*180/pi,ele2Ar*180/pi,dist2Ar); % Transformação de AER (esféricas) para ENU
%   xE2A = abs(xE2A);
%   yN2A = abs(yN2A);
xE2A = xE2A(1:ultima_amostra);
yN2A = yN2A(1:ultima_amostra);
zU2A = zU2A(1:ultima_amostra);

%==========================================================================
% Simulações de Monte Carlo
ENS=1;
delta_Tracking_total = tr2A(end) - tr2A(1);
%delta_Tracking = 0.2*delta_Tracking_total; % 4seg. de tracking\\TOTAL 23.85 Atenção: o tempo total de rastreio é 23,8s //10%= 2.3850   //20%= 4.7700 //25%= 5.9625// 30%=7.1550  // 35%=8.3475 //40%=9.5400 //45%= 10.7325 // 50%= 11.9250
delta_Tracking =1*delta_Tracking_total;%0.5*delta_Tracking_total;
amostra_parada = 20*delta_Tracking;
amostra_parada = round( amostra_parada);
tic
auxENS = zeros(ENS,1);
%Variáveis auxiliares para ponto de impcto previsto
xaux1 = zeros(ENS,1);
yaux1 = zeros(ENS,1);
zaux1 = zeros(ENS,1);
%Variáveis auxiliares para ponto de impcto estimado
xaux2 = zeros(ENS,1);
yaux2 = zeros(ENS,1);
zaux2 = zeros(ENS,1);
erro=zeros(1000,1);
erro_rms = zeros(ENS,1);
erro_rse = 0;
tempo_processamento = 0;


% SNR81_linear = 10.^(SNR81/10);
% potencia_media_x = sum(xE2A.^2) / length(xE2A);
% potencia_media_y = sum(yN2A.^2) / length(yN2A);
% potencia_media_z = sum(zU2A.^2) / length(zU2A);
% 
% N_x = potencia_media_x./SNR81_linear;
% N_y = potencia_media_y./SNR81_linear;
% N_z = potencia_media_z./SNR81_linear;
% 
% sigma_nx = sqrt(N_x);
% sigma_ny = sqrt(N_y);
% sigma_nz = sqrt(N_z);
% sigma_nx = sigma_nx(1:amostra_parada);
% sigma_ny = sigma_ny(1:amostra_parada);
% sigma_nz = sigma_nz(1:amostra_parada);

for ens=1:ENS
% 
% wx = sigma_nx.*randn(1,amostra_parada)';
% wx = mean(wx) + sigma_nx.*randn(1,amostra_parada)';
%  
% wy = sigma_ny.*randn(1,amostra_parada)';
% wy = mean(wy) + sigma_ny.*randn(1,amostra_parada)';
% 
% wz = sigma_nz.*randn(1,amostra_parada)';
% wz = mean(wz) + sigma_nz.*randn(1,amostra_parada)';
% 
wx = 0*wgn(length(xE2A), 1, 1, 'real');
wy = 0*wgn(length(xE2A), 1, 1, 'real');
wz = 0*wgn(length(xE2A), 1, 1, 'real');

% wx = 0*wgn(length(xE2A), 1, 0.001, 'real');
% wy = 0*wgn(length(xE2A), 1, 0.001, 'real');
% wz = 0*wgn(length(xE2A), 1, 0.001, 'real');

xE2A = xE2A + wx;
yN2A = yN2A + wy;
zU2A = zU2A + wz;

% Teste Matriz de Conversão
% xE2A2 = dist2Ar.*sin(ele2Ar);
% yN2A2 = dist2Ar.*sin(azi2Ar);
% zU2A2 = dist2Ar.*sqrt((1 - (sin(ele2Ar).^2) - (sin(azi2Ar).^2) ));


% 2) Filtro Média Móvel

tamJanela = 50 ;

bm = (1/tamJanela)*ones(1,tamJanela);
x_filtMM = filtfilt(bm,1,xE2A(1:end));
y_filtMM = filtfilt(bm,1,yN2A(1:end));
z_filtMM = filtfilt(bm,1,zU2A(1:end));
% 
% x_filtMM = filter(bm,1,xE2A);
% y_filtMM = filter(bm,1,yN2A);
% z_filtMM = filter(bm,1,zU2A);

% 
%  xMM2A = abs(x_filtMM);
%  yMM2A = abs(y_filtMM);
%  zMM2A = abs(z_filtMM);

% 
 xMM2A = x_filtMM(1:end);
 yMM2A = y_filtMM(1:end);
 zMM2A = z_filtMM(1:end);
 tr3A =tr2A(1:end);
% % 
% xMM2A = xE2A;
% yMM2A = yN2A;
% zMM2A = zU2A;


% Cálculo das Velocidades =================================================
% Método Comando diff
dx_dt = diff(xMM2A)./diff(tr3A);
dy_dt = diff(yMM2A)./diff(tr3A);
dz_dt = diff(zMM2A)./diff(tr3A);
% dx_dt = diff(xE2A)./diff(tr2A);
% dy_dt = diff(yN2A)./diff(tr2A);
% dz_dt = diff(zU2A)./diff(tr2A);
vxA = dx_dt;
vyA = dy_dt;
vzA = dz_dt; 
v = sqrt((vxA.^2)+(vyA.^2)+(vzA.^2));

v_ = [vxA';vyA';vzA'];


% Cálculo das Acelerações =================================================
dvx_dt = diff(vxA)./diff(tr3A(1:end-1));
dvy_dt = diff(vyA)./diff(tr3A(1:end-1));
dvz_dt = diff(vzA)./diff(tr3A(1:end-1));

ax = dvx_dt; 
ay = dvy_dt; 
az = dvz_dt; 
a = diff(v)./diff(tr3A(1:end-1)); % amodC = diff(vmodC)./diff(tr(1:end-1)');
%ax +ay +az;
a_ = [ax';ay';az'];


% Filtro de Kalman EKF=====================================================

T = 0.05;                % Período de Amostragem Radar (T =50ms)
t =tr2A(1:end);            %0: T: total_time; % vetor de tempo
n_estados = 7;               % Quantidade de estados / 7 estados e 3 medidas
po =  [xE2A(1); yN2A(1); zU2A(1)]; %  posição inicial do projétil
v_ = v_(:,1:end-1);          % velocidade inicial do projétil
tam=length(xE2A)*2;

% Inicialização das variáveis==============================================
alpha = zeros(1,tam);        % parâmetro de arrasto - alpha = S*Cd/m (m^2/kg)
alpha(1) = 10^-4;%0.7e-04;%3.9000e-04;% 10^-6;           % inicialização do parâmetro de arrasto
%lamb = zeros(1,tam);
 p0 = 1.203411;                % densidade local
%p0 = 0.002378;
h = 1.0361*10^-4;    %3.158*10^-5 (1/feet) = 1.0361*10^-4 (1/m) ---- ft = m x 3.28084
lamb(1) =alpha(1)*p0;  %9.8761e-05
m = 4.5; %m = 13.685;
d = 81.03e-3; %d = 119.530e-3;%
S = pi*(d^2)/4;
xk_k = zeros(n_estados,tam);   % vetor de estados x
f_xk = zeros(n_estados,tam);   % vetor de estados x' (derivada do vetor x)
Cd_est = zeros(tam,1);
% px = zeros(tam,1);
% py = zeros(tam,1);
% pz = zeros(tam,1);
% Vx = zeros(tam,1);
% Vy = zeros(tam,1);
% Vz = zeros(tam,1);
%Vlamb = zeros(tam,1);

% Matriz de covariâncias do ruído de medidas ============================
% Conforme as medidas radar, temos os seguintes valores de desvio padrão.
% =========================================================================
erroEle = medidas81T2_ADOUR(1:end,8)*180/pi;  %*180/pi;
erroAzi = medidas81T2_ADOUR(1:end,9)*180/pi;  %*180/pi;
erroDist = medidas81T2_ADOUR(1:end,10);
sigma_Ele = std(erroEle); sigma_Azi = std(erroAzi); sigma_Dist = std(erroDist);
% t2 = t2(1:length(t2)-110);
% 
% intervalo_tempo= [t2(1) t2(end)];

% % % % 
% figura13=figure(13);
% subplot(3,1,1)
% plot(t2,erroEle,'k','linewidth',1); title('Erro das Medidas AER do Radar','Interpreter','latex','FontSize',14);
% hold on; plot(t2,sigma_Ele*ones(1, length(t2))); stem(intervalo_tempo, [sigma_Ele(1) sigma_Ele(end)]); grid
% legEle = legend('$e_{\phi}(t)$','$\sigma_{\phi}=0,19^{\circ}$','Location','best');
% set(legEle,'Interpreter','latex');
% set(legEle,'FontSize',12);xlabel('tempo(s)','Interpreter','latex','FontSize',14); ylabel('Eleva\c c\~ao ($^\circ$)','Interpreter','latex','FontSize',16)  %legend('$\hat{\psi}$','Interpreter','latex')
% 
% subplot(3,1,2)
% plot(t2,erroAzi,'k','linewidth',1); %title('Erros medidas Azimute'); 
% hold on; plot(t2,sigma_Azi*ones(1, length(t2))); stem(intervalo_tempo,[sigma_Azi(1) sigma_Azi(end)]); grid
% legAzi = legend('$e_{\theta}(t)$','$\sigma_{\theta}=0,33^{\circ}$','Location','best');
% set(legAzi,'Interpreter','latex');
% set(legAzi,'FontSize',12);xlabel('tempo(s)','Interpreter','latex','FontSize',14); ylabel('Azimute($^\circ$)','Interpreter','latex','FontSize',16)
% 
% subplot(3,1,3)
% plot(t2,erroDist,'k','linewidth',1); %title('Erros medidas Distância');
% hold on; plot(t2,sigma_Dist*ones(1, length(t2))); stem(intervalo_tempo,[sigma_Dist(1) sigma_Dist(end)]); grid
% legDis = legend('$e_{R}(t)$','$\sigma_{R}=9.52 m$','Location','southwest');
% set(legDis,'Interpreter','latex');
% set(legDis,'FontSize',12);xlabel('tempo(s)','Interpreter','latex','FontSize',14); ylabel('Dist\^ancia (m)','Interpreter','latex','FontSize',16)
% exportgraphics(figura13, 'fig.medErro_Radar.pdf', 'ContentType', 'vector');
% figura2=figure(2);
% 
% plot(tr2A,VelDist,'k','linewidth',1); title('Medidas de Velocidade em Distancia do Radar','Interpreter','latex','FontSize',14);
% hold on;  grid %plot(tr2A,var(VelDist)*ones(1, length(tr2A)));
% % legEle = legend('$e_{Vel_Dist}(t)$');
% % set(legEle,'Interpreter','latex'); set(legEle,'FontSize',12);
% xlabel('tempo(s)','Interpreter','latex','FontSize',14); ylabel('Velocidade em Dist.(m/s)','Interpreter','latex','FontSize',16)  %legend('$\hat{\psi}$','Interpreter','latex')



%=========================================
% Matriz de covariância do Ruído da TRANSFORMAÇÃO das medidas -
% Representa a conversão mais o ruído
sig2E = var(erroEle);  sig2A = var(erroAzi);  sig2R = var(erroDist);  % variância dos erros AER

phi = ele2Ar; theta = azi2Ar; pho = dist2Ar;                            % medidas AER

sig2x = (cos(theta).^2) .* (sig2R.*cos(phi).^2 + (pho.^2).*sig2E.*sin(phi).^2) + (pho.^2).*sig2A.*(sin(theta).^2).*cos(phi).^2; % variância de x

sig2y = (sin(theta).^2) .* (sig2R.*cos(phi).^2 + (pho.^2).*sig2E.*sin(phi).^2) + (pho.^2).*sig2A.*(cos(theta).^2).*cos(phi).^2; % variância de y

sig2z = (sig2R.*sin(phi).^2 + (pho.^2).*sig2E.*cos(phi).^2); % variância de z

sigxy = 0.5*sin(2*theta).*((sig2R-(pho.^2).*sig2A).*cos(phi).^2 + (pho.^2)*sig2E.*sin(phi).^2); % covariância de xy

sigxz = (0.5*cos(theta).*sin(2*phi).*(sig2R-(pho.^2).*sig2E)); % covariância de xz

sigyz = (0.5*sin(theta).*sin(2*theta).*(sig2R-(pho.^2).*sig2E)); % covariância de yz


% 
% sig2x = abs(sig2x);
% 
% sig2y = abs(sig2y);
% 
% sig2z = abs(sig2z);
% 
% sigxy = abs(sigxy);
% 
% sigxz = abs(sigxz);
% 
% sigyz = abs(sigyz);

% 
% 
 R = zeros(3,3); % Matriz de Covariância do erro das medidas (inclui os erros da transformação de coordenadas e ruídos/incertezas dos sensores de medição AER)
sig2x = 10*ones(1,tam);
sig2y = sig2x;
sig2z = sig2x;
sigxy = 1*rand(1,tam);
sigxz = sigxy;
sigyz = sigxy;

R_o = [sig2x(1) sigxy(1) sigxz(1); 
       sigxy(1) sig2y(1) sigyz(1);
       sigxz(1) sigyz(1) sig2z(1)]; 
% =================== EKF =================================================
% Inicialização do vetor de estados, 7 estados: x(0|0) = [x(0) y(0) z(0) vx(0) vy(0) vz(0) CD(0)]
%v_(:,1) = [100;100;100];

xk_k(:,1)  = [po; v_(:,1) ; lamb(1)];    % Definir estimativa inicial do estado
sig2v= (150)^2;       %       sig2v = var(VelDist)
   
%======== Modelagem Coeficiente de Arrasto como Processo de Wiener =========
% Nw = tam;
% Tw = T;
% dt = Tw/Nw;
% dW = zeros(1,Nw);               % preallocate arrays ...
% W = zeros(1,Nw);                % for efficiency
% dW(1) = sqrt(dt)*randn;         % first approximation outside the loop ...
% W(1) = dW(1);                   % since W(0) = 0 
% for jw = 2:Nw
% dW(jw) = sqrt(dt)*randn;        % general increment
% W(jw) = W(jw-1) + dW(jw);
% end
%alpha = 10^-4*abs(W);
%alpha =(10^-4)*exp(alpha)
%beta = alpha*diff(abs(W))/T;
%sig2alpha =  var(alpha);
%sig2alpha =  (10^-4)^2;% 0.001^2;
sig2lamb =  10^-11;%0.00001^2; %(6.0907e-06)^2;7.4030e-04;

% Inicialização da Matriz de covariância de estimativa do estado Re_o(0|0) ou P (0|0) - Definir covariância de erro inicial 
P_k_k = [R_o zeros(3,3) zeros(3,1); %(7x7)
      zeros(3,3) sig2v*eye(3,3) zeros(3,1);
      zeros(1,3) zeros(1,3) sig2lamb];


% Vetor das medidas reais==================================================
%z_med = [xE2A(1:tam) yN2A(1:tam) zU2A(1:tam)]'; % Medidas
z_med = [xE2A yN2A zU2A]'; % Medidas
H=[eye(3) zeros(3,4)];  % -  Dim:(3x7) Matriz de saída H
% Matriz de Ganho T  -  Dim: (7 x 4)
T7 = [((T^2)/2)*eye(3)    zeros(3,1);
       T*eye(3)         zeros(3,1); 
       zeros(1,3)         1       ];


% Matriz F de transição de estados=========================================
%==========================================================================
vv = (10^-7)*wgn(4, tam, 1); % Ruído do processo v  0.00001 ,'real'
%w  = wgn(3, tam, 1); % Ruído das medidas w
%w =  0*wgn(3, tam, 10); %w = 10*wgn(3, tam, 10 ,'dBm');
w  = [wx wy wz]';
%w=zeros(3,tam);
%w  = 30*randn(3, tam, 40);
% =========================================================================
% Algorítmo EKF ===========================================================

% Sejam as seguintes considerações para os sete estados:
% xk_k(1,k) = x(k)  (posição na coordenada x) 1º elemento do vetor de estado x(k|k)
% xk_k(2,k) = y(k)  (posição na coordenada y) 2º elemento do vetor de estado x(k|k)
% xk_k(3,k) = z(k)  (posição na coordenada z) 3º elemento do vetor de estado x(k|k)
% xk_k(4,k) = vx(k) (velocidade na coordenada x) 4º elemento do vetor de estado x(k|k)
% xk_k(5,k) = vy(k) (velocidade na coordenada y) 5º elemento do vetor de estado x(k|k)
% xk_k(6,k) = vz(k) (velocidade na coordenada z) 6º elemento do vetor de estado x(k|k)
% xk_k(7,k) = alpha(k) (parâmetro de arrasto) 7º elemento do vetor de estado x(k|k)

VV = zeros(1,tam);
paux =zeros(tam,1);
for k = 1:tam
   
% Predição do Estado
V_ = [xk_k(4,k) xk_k(5,k) xk_k(6,k)]; % posições das velocidades vx vy e vz
V = norm(V_);
% For Densidade do ar nível do mar Army Standard Metro:
%                             po = 0.0751265 lbs./cubic foot = 1.203411 kg/m3
%                                  1 kg/m3 =   0.06242796 lb/ft3
%                             ///  1 lb/ft3 = 16.02 kg/m3
VV(k)=V;
%                             h = 3.158*10^-5 (1/feet) = 1.0361*10^-4 (1/m)
g = 9.80665;                  % aceleração da gravidade m/s^2
p0 = 1.203411;                % densidade local
%p0 = 0.002378;
h = 1.0361*10^-4; %h = 2.926*10^-5 + (10^-10*xk_k(3,k-1))*3.28084;
p = p0*exp(-h*xk_k(3,k));      % modelo exponencial da densidade do ar
paux(k) = p;
Pp = -(1/2)*p*V;            % variável P para auxilar nos cálculos

% R = [sig2x(k) sigxy(k) sigxz(k);
%       sigxy(k) sig2y(k) sigyz(k);
%       sigxz(k) sigyz(k) sig2z(k)];

% Matriz de Covariância do Processo
% Q_k = [sqrt(sig2x(k-1))        0          0        0;
%       0           sqrt(sig2y(k-1))        0        0;
%       0              0         sqrt(sig2z(k-1))    0;
%       0              0            0     sqrt(sig2alpha)];



% Q_k = [var(xk_k(1,:))        0          0        0; %10^-8
%       0           var(xk_k(2,:))         0        0;
%       0              0         var(xk_k(3,:))     0;
%       0              0            0     var(alpha)]; %5*10^-8



q_x = 10^-3;  %0.9*10^3;q_x = 10^-3; 3.8^2;1.0000e-06
q_y =q_x;
q_z =q_x;
q_lamb =1e-12;%*exp(-0.00005*xk_k(3,k));%10^-9;
%q_lamb =  (7.25*10^-8)*exp(-0.00005*xk_k(3,k));%10^-6; %  1.5600e-04;   1*10^-4

Q_k = [q_x        0          0        0;  % 10^-8
      0           q_y        0        0;
      0           0          q_z      0;
      0           0          0        q_lamb ]; %5*10^-8 %(7.25*10^-13)*exp(-0.00005*xk_k(3,k-1))

Q = T7*Q_k*T7';
% Q =[     0         0         0         0         0         0         0
%          0         0         0         0         0         0         0
%          0         0         0         0         0         0         0
%          0         0         0         q_x       0         0         0
%          0         0         0         0         q_y       0         0
%          0         0         0         0         0         q_z       0
%          0         0         0         0         0         0         q_lamb ];
% Cálculos dos elemetos da Matriz Jacobiana F em x_est(k|k)============================================================
 
vx = xk_k(4,k);      % Os elementos de F são dependentes das componentes da velocidade
vy = xk_k(5,k);
vz = xk_k(6,k);
%alpha(k) = xk_k(7,k);
lamb = xk_k(7,k);
%alpha(k-1) = xk_k(7,k-1);
ee = 0;
%===== Elementos da Matriz Jacobiana===================================================================================

% f43 = vx*V*h*lamb/2;   % xk_k(7,k) = alpha(k)
% f44 = -(lamb/2)*(V+(vx^2)/(V+ee));
% f45 = -lamb*vx*vy/(2*V+ee);
% f46 = -lamb*vx*vz/(2*V+ee);
% f47 = -V*vx/2;
% 
% f53 = vy*V*h*lamb/2;   % xk_k(7,k) = alpha(k)
% f54 = -lamb*vy*vx/(2*V+ee);
% f55 = -(lamb/2)*(V+(vy^2)/(V+ee));
% f56 = -lamb*vy*vz/(2*V+ee);
% f57 = -V*vy/2;
% 
% f63 = vz*V*h*lamb/2;   % xk_k(7,k) = alpha(k)
% f64 = -lamb*vx*vz/(2*V+ee);
% f65 = -lamb*vy*vz/(2*V+ee);
% f66 = -lamb*(vz^2)/(2*V+ee);
% f67 = -V*vz/2;
% 
% f73 = (h^2)*vz*lamb;   % xk_k(7,k) = alpha(k)
% f74 = 0;
% f75 = 0;
% f76 = -lamb*h;
% f77 = -h*vz;






f43 = vx*V*h*xk_k(7,k)/2;   % xk_k(7,k) = alpha(k)
%f43 = 0;
f44 = -(xk_k(7,k)/2)*(V+((vx^2)/(V+ee)));
f45 = -xk_k(7,k)*vx*vy/(2*V+ee);
f46 = -xk_k(7,k)*vx*vz/(2*V+ee);
f47 = -V*vx/2;

f53 = vy*V*h*xk_k(7,k)/2;   % xk_k(7,k) = alpha(k)
%f53 = 0;
f54 = -xk_k(7,k)*vy*vx/(2*V+ee);
f55 = -(xk_k(7,k)/2)*(V+((vy^2)/(V+ee)));
f56 = -xk_k(7,k)*vy*vz/(2*V+ee);
f57 = -V*vy/2;

f63 = vz*V*h*xk_k(7,k)/2;   % xk_k(7,k) = alpha(k)
%f63 = 0;
f64 = -xk_k(7,k)*vx*vz/(2*V+ee);
f65 = -xk_k(7,k)*vy*vz/(2*V+ee);
f66 = -xk_k(7,k)*((V^2)+vz^2)/(2*V+ee);
f67 = -V*vz/2;

f73 = (h^2)*vz*xk_k(7,k);   % xk_k(7,k) = alpha(k)
%f73 = 0;
f74 = 0;
f75 = 0;
f76 = -xk_k(7,k)*h;
f77 = -h*vz;


F_=[f43 f44 f45 f46 f47;
    f53 f54 f55 f56 f57;
    f63 f64 f65 f66 f67;
    f73 f74 f75 f76 f77];

F = [zeros(3) eye(3) zeros(3,1);
     zeros(4,2) F_];

% ESPAÇO de Estados x' (corresponde a derivada em relação tempo do vetor de estado x)
x_1 = xk_k(4,k);
x_2 = xk_k(5,k);
x_3 = xk_k(6,k);
x_4 = -(1/2)*xk_k(7,k)*V*x_1;
x_5 = -(1/2)*xk_k(7,k)*V*x_2;
x_6 = -(1/2)*xk_k(7,k)*V*x_3 - g;
x_7 = -xk_k(7,k)*h*x_3;

f_xk = [x_1; x_2; x_3; x_4; x_5; x_6; x_7];
%xk_k(:,k) = xk_k(:,k-1) + f_xk*T +  F*f_xk*(T^2)/2 + T7*vv(:,k-1);
%xkm1_k = xk_k(:,k-1) + f_xk*T +  F*f_xk*(T^2)/2 + T7*vv(:,k-1);

   fi_k = eye(7) + F*T;
   P_km1_k = fi_k*P_k_k*fi_k'+Q; 
%  

xkm1_k = xk_k(:,k) + f_xk*T + F*f_xk*(T^2)/2+T7*vv(:,k);
%xkm1_k = xk_k(:,k) + f_xk*T;
 %P_km1_k = F*P_k_k*F'+Q;
 %P_km1_k = F*(P_k_k+Q)*F'; 
% Observações
%z_km1 = H*xk_k(:,k) + w(:,k-1); % z_km1(:,k) = H*xkm1_k;
 z_km1 = H*xkm1_k;% + w(:,k); % z_km1(:,k) = H*xkm1_k;
% Atualização medidas
  if k>=amostra_parada
     K=zeros(7,3); 
     
     xkm1_km1 = xkm1_k; 
  else   
      R = [sig2x(k) sigxy(k) sigxz(k);
      sigxy(k) sig2y(k) sigyz(k);
      sigxz(k) sigyz(k) sig2z(k)];
   K = P_km1_k*H'/(H*P_km1_k*H' + R);    % K = P_km1_k*H'*inv(H*P_km1_k*H' + R);               % Cálculo do ganho de Kalman
   xkm1_km1 = xkm1_k + K*(z_med(:,k) - z_km1); 
   
   
 end


%K = P_km1_k*H'/(H*P_km1_k*H' + R);     % K = P_km1_k*H'*inv(H*P_km1_k*H' + R);               % Cálculo do ganho de Kalman

%xk_k(:,k) = xk_k(:,k) + K*(z_med(:,k) - z_km1);
%xkm1_km1 = xkm1_k + K*(z_med(:,k) - z_km1);    % Atualiza os estados estimados
%P_km1_k = (eye(n_estados)-K*H)*P_km1_k; 
P_km1_km1 = (eye(n_estados)-K*H)*P_km1_k;           % Atualiza a matriz de Covariancia do erro
%P_km1_km1 = (eye(n_estados)-K*H)*P_km1_k*(eye(n_estados)-K*H)'+K*R*K'; 

xk_k(:,k+1) = xkm1_km1;
P_k_k = P_km1_km1;

% px(k)  = xk_k(1,k);
% py(k)  = xk_k(2,k);
% pz(k)  = xk_k(3,k);
% Vx(k)  = xk_k(4,k);
% Vy(k)  = xk_k(5,k);
% Vz(k)  = xk_k(6,k);
% Vlamb(k)  = xk_k(7,k);
%Cd_est(k) = Vlamb(k)*m/(S.*paux(k));

px(k)  = xkm1_km1(1);
py(k)  = xkm1_km1(2);
pz(k)  = xkm1_km1(3);
Vx(k)  = xkm1_km1(4);
Vy(k)  = xkm1_km1(5);
Vz(k)  = xkm1_km1(6);
Vlamb(k)  = xkm1_km1(7);
% Valpha(k)  = xkm1_km1(7);
Cd_est(k) = Vlamb(k)*m/(S.*paux(k)); %paux(k)
tk(k) = t(1) +k*T;

% px(k)  = xkm1_k(1);
% py(k)  = xkm1_k(2);
% pz(k)  = xkm1_k(3);
% Vx(k)  = xkm1_k(4);
% Vy(k)  = xkm1_k(5);
% Vz(k)  = xkm1_k(6);
% Valpha(k)  = xkm1_k(7);
%alpha(k) = xk_k(7);

if pz(k)<0
        iendek = k-1;
break
        % pz(k) = pz(k-1);
        % py(k) = py(k-1);
        % px(k) = px(k-1);
          
end
end
px = px(1:iendek);    py = py(1:iendek);  pz = pz(1:iendek); Cd_est = Cd_est(1:iendek); tk = tk(1:iendek);
% pz = nonzeros(pz);
% py = nonzeros(py);
% px = nonzeros(px);
%fig2=figure(2); % Gráfico dos Estados

% subplot(511);
% plot(t,xE2A(1:tam),'k','Linewidth',1); grid
% hold on
% plot(t(2:end),px(2:end),'k--','Linewidth',1);
% hold on
% plot(t,yN2A(1:tam),'r','Linewidth',1);
% hold on
% plot(t(2:end),py(2:end),'r--','Linewidth',1);
% hold on
% plot(t,zU2A(1:tam),'b','Linewidth',1);
% hold on
% plot(t(2:end),pz(2:end),'b--','Linewidth',1);
% title('Posições');
% xlabel('tempo(s)');
% legend('x_{Medido}','x_{Real}');
% 
% subplot(512);
% plot(t,Vx,'k','Linewidth',1); grid
% hold on;
% plot(t,Vy,'r','Linewidth',1);
% hold on;
% plot(t,Vz,'b','Linewidth',1);
% title('Velocidades Estimadas');
% xlabel('tempo(s)');
% legend('v_x','v_y','v_z' );
% fig2=figure(2); % Gráfico dos Estados
% subplot(311);
% %plot(t,alpha,'k','Linewidth',1);
% %hold on;
% plot(t(1:end),Vlamb(1:end),'g','Linewidth',1); grid
% title('Parâmetro de Arrasto');
% xlabel('tempo(s)');
% legend('$\hat{\alpha}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %exportgraphics(fig2, 'fig_kalman0.pdf', 'ContentType', 'vector');
% 
% subplot(312);
% %plot(t,alpha,'k','Linewidth',1);
% %hold on;
% plot(t(1:end),1./(Vlamb(1:end)),'b','Linewidth',1); grid
% title('Parâmetro Balístico');
% xlabel('tempo(s)');
% legend('$\hat{\beta}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %exportgraphics(fig2, 'fig_kalman0.pdf', 'ContentType', 'vector');
% 
% subplot(313);
% 
% plot(t(1:end),Vlamb(1:end)./0.0013,'k-','Linewidth',1); grid; hold on; 
% plot(t(1:end),Cd_est,'r-','Linewidth',1);
% title('Coeficiente de Arrasto');
% xlabel('tempo(s)');
% legend('$\hat{C}_{D}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %legend('\hat{C}_{D}');


% fig3=figure(3);
% 
% plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% plot3(xE2A(1:tam),yN2A(1:tam),zU2A(1:tam),'k','Linewidth',1); hold on;
% %plot3(xE2A2(1:tam),yN2A2(1:tam),zU2A2(1:tam),'k--','Linewidth',1); 
% %hold on;
% plot3(px(2:end),py(2:end),pz(2:end),'b-','Linewidth',1);
% % grid
% % xlabel('x (m)')
% % ylabel('y (m)')
% % zlabel('z (m)')
% % AZIMUTE=-74.8675; ELEVACAO=15.9026;    
% %  view(AZIMUTE, ELEVACAO);
% title('Trajetória 3D');
% legend('Medida','Estimada');
% legend('Lan\c cadora','Medidas $P_1$','Traj. Estimada $P_1$','Location','northeast','Interpreter','latex','FontSize',12);grid
% %legend('Modelo $P_1$','Medidas','','Lan\c cadora','Ponto de impacto','Interpreter','latex','FontSize',12);grid 
% xlabel('x (m)','Interpreter','latex','FontSize',12)
% ylabel('y (m)','Interpreter','latex','FontSize',12)
% zlabel('z (m)','Interpreter','latex','FontSize',12)
% title('Trajet\''orias 3D Estimadas','Interpreter','latex','FontSize',13)
% AZIMUTE=5.5678; ELEVACAO=-10.3897; 
% exportgraphics(fig3, 'fig_EKF3.pdf', 'ContentType', 'vector');



% fig4=figure(4); 
% subplot(2,1,1); % subplot(2,2,[3,4]);
% plot(t(1:end),Vlamb(1:end).*(m/(S*paux)),'k-','Linewidth',1); grid
% title('Coeficiente de Arrasto ','Interpreter','latex','FontSize',13);
% xlabel('tempo(s)','Interpreter','latex','FontSize',12) 
% legend('$\hat{C}_{D}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% 
% 
% subplot(2,1,2);
% plot(t(1:end),Vlamb(1:end).*(m/(S*paux)),'k-','Linewidth',1); grid
% title('Coeficiente de Arrasto ','Interpreter','latex','FontSize',13);
% xlabel('tempo(s)','Interpreter','latex','FontSize',12); %axis([t(1) t(end) 0 0.2])
% legend('$\hat{C}_{D}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% exportgraphics(fig4, 'fig_EKF4.pdf', 'ContentType', 'vector');

% %==========================================================================
% %MÉTODO DE RUNGE KUTTA - PARA EXTRAPOLAR TRAJETÓRIA
% h1e=0.05;    % h = t(i+1)-t(i)                                         % tamanho do passo
% tame=tam;
% Ne = tame*h1e;
% %tam = length(v);
% %tam = 1000;
% t1e = 0:h1e:Ne-1;                  % tamanho do intervalo de t
% t = t1e + tr2A(1);
% %v =v(1:end-1);
% vx1e = zeros(1,tame);
% vy1e = zeros(1,tame);
% vz1e = zeros(1,tame);
% vMe= zeros(1,tame);
% 
% % vMe(1)= sqrt((Vx(1)^2)+(Vy(1)^2)+(Vz(1)^2));
% % deltaX = (px(2)-px(1));
% % deltaY = (py(2)-py(1));
% % deltaZ = (pz(2)-pz(1));
% 
% %Módulo da Velocidade Inicial Estimada
% vMe(1)= sqrt((vxA(1)^2)+(vyA(1)^2)+(vzA(1)^2));
% deltaX = (xMM2A(2)-xMM2A(1));
% deltaY = (yMM2A(2)-yMM2A(1));
% deltaZ = (zMM2A(2)-zMM2A(1));
% %Angulos de elevação e Azimute
% ANG1e = deltaZ/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2)));
% theta_0e = asin(ANG1e); 
% ANG2 = deltaY/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2))*cos(theta_0e));
% phi_0 = -acos(ANG2); 
% %Componentes da Velocidade 
% % (Teste)
% % theta_0e = 252;
% % phi_0 = 34.72;
% theta_0e = 4.3810;
% phi_0 = 0.6060;
% %
% % vx1e(1) = vMe(1)*cos(theta_0e)*sin(phi_0); % condições iniciais
% % vy1e(1) = vMe(1)*cos(theta_0e)*cos(phi_0); %vy(1) = -0.268;
% % vz1e(1) = vMe(1)*sin(theta_0e);             % condições iniciais
% 
% % vx1e(1) = Vx(1); % condições iniciais da velocidade
% % vy1e(1) = Vy(1); %vy(1) = -0.268;
% % vz1e(1) = Vz(1);             % condições iniciais
% vx1e(1) = vxA(1); % condições iniciais da velocidade
% vy1e(1) = vyA(1); %vy(1) = -0.268;
% vz1e(1) = vzA(1);             % condições iniciais
% x1e = zeros(1,tame);
% y1e = zeros(1,tame);
% z1e = zeros(1,tame);
% % x1e(1) = px(1) ;     % condições iniciais  -1.4953   -0.2465    1.2056
% % y1e(1) = py(1);      % condições iniciais
% % z1e(1) = pz(1);      % condições iniciais
% % 
% x1e(1) = xMM2A(1);     % condições iniciais  -1.4953   -0.2465    1.2056
% y1e(1) = yMM2A(1);      % condições iniciais
% z1e(1) = zMM2A(1);      % condições iniciais
% 
% % Cd_e = Cd_est(150)*exp(W);
% % C_Dxe = p*S*Cd_e/(2*m);
% C_Dxe = paux'.*S.*Cd_est.*ones(1,length(paux))/2*m;%.*ones(1,length(Cd_est))/(2*m);
% %C_Dxe = mean(Vlamb(75:272)).*ones(1,length(Cd_est))/2;
% C_Dxe = Vlamb/2;%Vlamb.*ones(1,length(Cd_est))/2;
% C_Dye = C_Dxe;
% C_Dze = C_Dye;
% 
% 
% F_vvx = @(C_Dxe,vMe,vx1e) -C_Dxe.*vMe.*vx1e';       % ax
% F_vvy = @(C_Dye,vMe,vy1e) -C_Dye.*vMe.*vy1e' ;      % ay
% F_vvz = @(C_Dze,vMe,vz1e) -C_Dze.*vMe.*vz1e' - g ;  % az           % função EDO 
% 
% F_vx = @(vx1e) vx1e; % vx
% F_vy = @(vy1e) vy1e; % vy
% F_vz = @(vz1e) vz1e; % vz
% 
% 
% for i=1:(length(C_Dxe))                              % acelerações - velocidades comp X
%     k1 = h1e*F_vvx(C_Dxe(i),vMe(i),vx1e(i));
%     k2 = h1e*F_vvx(C_Dxe(i),vMe(i),vx1e(i)+0.5*k1);
%     k3 = h1e*F_vvx(C_Dxe(i),vMe(i),vx1e(i)+0.25*(k1+k2));
%     k4 = h1e*F_vvx(C_Dxe(i),vMe(i),vx1e(i)-(k2+2*k3));
%     vx1e(i+1) = vx1e(i) + (1/6)*(k1+4*k3+k4);        
% 
%     k1_ = h1e*F_vx(vx1e(i));                             % velocidades - posições comp X
%     k2_ = h1e*F_vx(vx1e(i)+0.5*k1_);
%     k3_ = h1e*F_vx(vx1e(i)+0.25*(k1_+k2_));
%     k4_ = h1e*F_vx(vx1e(i)-(k2_+2*k3_));
% 
%     x1e(i+1) = x1e(i) + (1/6)*(k1_+4*k3_+k4_);  
% 
% %==========================================================================
% 
%     Q1 = h1e*F_vvy(C_Dye(i),vMe(i),vy1e(i));            % acelerações - velocidades comp Y
%     Q2 = h1e*F_vvy(C_Dye(i),vMe(i),vy1e(i)+0.5*Q1);
%     Q3 = h1e*F_vvy(C_Dye(i),vMe(i),vy1e(i)+0.25*(Q1+Q2));
%     Q4 = h1e*F_vvy(C_Dye(i),vMe(i),vy1e(i)-(Q2+2*Q3));
%     vy1e(i+1) = (vy1e(i) + (1/6)*(Q1+4*Q3+Q4));      
% 
%     Q1_ = h1e*F_vy(vy1e(i));                            % velocidades - posições comp Y
%     Q2_ = h1e*F_vy(vy1e(i)+0.5*Q1_);
%     Q3_ = h1e*F_vy(vy1e(i)+0.25*(Q1_+Q2_));
%     Q4_ = h1e*F_vy(vy1e(i)-(Q2_+2*Q3_));
% 
%     y1e(i+1) = y1e(i) + (1/6)*(Q1_+4*Q3_+Q4_);
% 
% %==========================================================================
% 
%     J1 = h1e*F_vvz(C_Dze(i),vMe(i),vz1e(i));                % acelerações - velocidades comp Z
%     J2 = h1e*F_vvz(C_Dze(i),vMe(i),vz1e(i)+0.5*J1);
%     J3 = h1e*F_vvz(C_Dze(i),vMe(i),vz1e(i)+0.25*(J1+J2));
%     J4 = h1e*F_vvz(C_Dze(i),vMe(i),vz1e(i)-(J2+2*J3));
%     vz1e(i+1) = vz1e(i) + (1/6)*(J1+4*J3+J4);    
% 
%     J1_ = h1e*F_vz(vz1e(i));                                % velocidades - posições comp Z
%     J2_ = h1e*F_vz(vz1e(i)+0.5*J1_);
%     J3_ = h1e*F_vz(vz1e(i)+0.25*(J1_+J2_));
%     J4_ = h1e*F_vz(vz1e(i)-(J2_+2*J3_));
% 
%     z1e(i+1) = z1e(i) + (1/6)*(J1_+4*J3_+J4_);
% 
% vMe(i+1) = sqrt((vx1e(i+1)^2)+(vy1e(i+1)^2)+(vz1e(i+1)^2)); 
% 
%     if z1e(i)<0
%         iende = i-1;
% 
%         z1e(i) = z1e(i-1);
%         y1e(i) = y1e(i-1);
%         x1e(i) = x1e(i-1);
%         x1e = x1e(1:iende);    y1e = y1e(1:iende);  z1e = z1e(1:iende);  vMe= vMe(1:iende); 
% break
% 
%     end
% 
% end

x1e = px;
y1e = py;
z1e = pz;

xaux1(ens) = x1e(end);
yaux1(ens) = y1e(end);
zaux1(ens) = z1e(end);

xaux2(ens) = px(end);
yaux2(ens) = py(end);
zaux2(ens) = pz(end);


%erro = abs((C_d(1:length(Cd_est))'-Cd_est')./C_d(1:length(Cd_est))')*100 ; %+erro(1:length(Cd_est))
%erro_rse=sqrt((C_d(1:length(Cd_est))' - Cd_est').^2) + erro_rse; %sqrt(mean((observed - predicted).^2));
%erro_aux = erro+erro_aux
%erroX = abs((xMM2A(1:length(x1))-x1')./xMM2A(1:length(x1)))*100 + erroX;
tempo_processamento = toc + tempo_processamento;
end


% Pontos de Impacto EKF completo
matriz_aux2 = [xaux2 yaux2];
[idx,cord_ponto_central] = kmeans(matriz_aux2,1);
% 
% figPIP_EKF_completo = figure('Position', [100, 100, 900, 700]);
% plot(xaux2,yaux2, 'r.','linewidth',5); hold on; grid
% plot(cord_ponto_central(1,1),cord_ponto_central(1,2),'b+','linewidth',7)
% xlabel('x (m)','Interpreter','latex','FontSize',15)
% ylabel('y (m)','Interpreter','latex','FontSize',15)
% legend('Pontos Extrapolados \textit{open-loop} EKF','Centroide (algoritmo \textit{K-Means})','Interpreter','latex','FontSize',16);
% title('Trajet\''orias 3D Estimadas','Interpreter','latex','FontSize',16)
% %Fim do Loop Monte Carlo
% %erro = erro/ENS;
% %erro_rms=sum(erro_rms)/ENS;
 tempo_processamento = tempo_processamento/ENS

% fig3=figure('Position', [100, 100, 700, 500]);
% 
% plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k.','Linewidth',0.5); hold on;
% %plot3(xE2A2(1:tam),yN2A2(1:tam),zU2A2(1:tam),'k--','Linewidth',1); 
% %hold on;
% plot3(px(2:end),py(2:end),pz(2:end),'b-','Linewidth',1);
% % grid
% % xlabel('x (m)')
% % ylabel('y (m)')
% % zlabel('z (m)')
% % AZIMUTE=-74.8675; ELEVACAO=15.9026;    
% %  view(AZIMUTE, ELEVACAO);
% title('Trajetória 3D');
% legend('Medida','Estimada');
% legend('Lan\c cadora','Medidas $P_1$','Traj. Estimada $P_1$','Location','northeast','Interpreter','latex','FontSize',12);grid
% %legend('Modelo $P_1$','Medidas','','Lan\c cadora','Ponto de impacto','Interpreter','latex','FontSize',12);grid 
% xlabel('x (m)','Interpreter','latex','FontSize',12)
% ylabel('y (m)','Interpreter','latex','FontSize',12)
% zlabel('z (m)','Interpreter','latex','FontSize',12)
% title('Trajet\''orias 3D Estimadas','Interpreter','latex','FontSize',13)
% AZIMUTE=5.5678; ELEVACAO=-10.3897; 
% exportgraphics(fig3, 'fig_EKFREAL.pdf', 'ContentType', 'vector');

% 
% fig5 = figure('Position', [100, 100, 900, 700]);
% plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k.','Linewidth',0.1); hold on;
% %plot3(px(1:end),py(1:end),pz(1:end),'b-','Linewidth',1); hold on;
% plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r.','Linewidth',1); hold on;
% plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
% title('Trajetória 3D Extrapolada');
% legend('Medida','Estimada');
% legend('Lan\c cadora','Medidas $P_2$','Traj. Extrapolada $P_2$','PIP estimado $P_2$','Location','northeast','Interpreter','latex','FontSize',14);grid
% %legend('Modelo $P_1$','Medidas','','Lan\c cadora','Ponto de impacto','Interpreter','latex','FontSize',12);grid 
% xlabel('x (m)','Interpreter','latex','FontSize',14)
% ylabel('y (m)','Interpreter','latex','FontSize',14)
% zlabel('z (m)','Interpreter','latex','FontSize',14)
% title('Trajet\''orias 3D Estimadas','Interpreter','latex','FontSize',15)
% AZIMUTE= 11.0613; ELEVACAO=-80.7600; 
% view(ELEVACAO, AZIMUTE);
% exportgraphics(fig5, 'fig_PIP_EKF_REAL_120.pdf', 'ContentType', 'vector');
% 
% 
% 



%======================== Plotando TLE  ===================================
% x_ipV = -5.56782e+03; %cord_ponto_central(1,1); 
% y_ipV = -1.78891e+03; %cord_ponto_central(1,2);

%Tiro 1========================================
% x_ipV = -5.49604e+03; %cord_ponto_central(1,1); 
% y_ipV = -1.75012e+03; %cord_ponto_central(1,2);

%Tiro 2========================================
x_ipV = -5.57301e+03; %cord_ponto_central(1,1); 
y_ipV = -1.79074e+03; %cord_ponto_central(1,2);


% 
% x_ip = x_ipV; 
% y_ip = y_ipV;


x_ip = cord_ponto_central(1,1); 
y_ip = cord_ponto_central(1,2);


x_pd = xaux1;
y_pd = yaux1;

% x_pd = rmmissing(xaux1);
% y_pd = rmmissing(yaux1);
% 
% x_ip = mean(x_pd);
% y_ip = mean(y_pd);

%mu_x ==x_ip;
%mu_y ==y_ip;
% var_x_pd = var(x_pd);
% cov_x_pd_y_pd = sum((x_pd - mean(x_pd)) .* (y_pd - mean( y_pd))) / (length(x_pd) - 1);
% var_y_pd =  var(y_pd);


vetor_ip = [x_ip, y_ip];  % Vetor de médias %meanVector = [mu_x,mu_y];
%covarianceMatrix = [var_x_pd, cov_x_pd_y_pd; cov_x_pd_y_pd, var_y_pd];  % Matriz de covariância
% erro_x =  abs(x_ip - x_pd);
% erro_y =  abs(y_ip - y_pd);

erro_x =  (x_ip - x_pd);
erro_y =  (y_ip - y_pd);

%covarianceMatrix = cov(x_pd,y_pd); 
covarianceMatrix = cov(erro_x,erro_y); 
 [V1, D1] = eig(covarianceMatrix);
 lambda_max = max(diag(D1));
% 
 covarianceMatrixC = [lambda_max, 0; 0, lambda_max];
% Nível de confiança (por exemplo, 0.95 para 95% de confiança)
confidenceLevel = 0.991;
confidenceLevelC = 0.99;


% % Calcula os autovalores e autovetores da matriz de covariância
% [V, D] = eig(covarianceMatrix);
% [Vc, Dc] = eig(covarianceMatrixC);

% Calcula os comprimentos dos semi-eixos
% semiAxisLengths = sqrt(diag(D) * chi2inv(confidenceLevel, 2));
% semiAxisLengthsC = sqrt(diag(Dc) * chi2inv(confidenceLevelC, 2));
% % Gera pontos da elipse
% theta = linspace(0, 2*pi, 10000);
% %points = [cos(theta'), sin(theta')] * (V * sqrtm(D) * diag(semiAxisLengths));
% points = [cos(theta'), sin(theta')] * (V * diag(semiAxisLengths));
% pointsC = [cos(theta'), sin(theta')] * (Vc * sqrtm(Dc) * diag(semiAxisLengthsC));
% % Translada e rotaciona a elipse para a posição correta
% ellipsePoints = bsxfun(@plus, points, vetor_ip);
% ellipsePointsC = bsxfun(@plus, pointsC, vetor_ip);

%======================== Cáculo do TLE  ==================================
r_TLE_s = sqrt((x_pd-x_ipV).^2+(y_pd-y_ipV).^2);
r_TLE_s_ord = sort(r_TLE_s);
TLE_porcent = 85;
TLE_amostra_porcent = (TLE_porcent*length(x_pd))/100;
TLE_amostra_porcent = round(TLE_amostra_porcent);
%TLE_amostra_porcent = (TLE_porcent*ens)/100;
r_TLE_final = r_TLE_s_ord(TLE_amostra_porcent );

theta_TLE = linspace(0,2*pi);
x_TLE = r_TLE_final*cos(theta_TLE) + x_ipV;
y_TLE = r_TLE_final*sin(theta_TLE) + y_ipV;

%======================== Cáculo 10% alcance máximo do projétil  ==================================
alc_max = sqrt((x_ipV)^2+(x_ipV)^2);

x_max = alc_max*cos(theta_TLE) + x_ipV;
y_max = alc_max*sin(theta_TLE) + y_ipV;




%============
% Calculo dos vetores de erro e medidas de dispersão

% x_pd = xaux1;
% y_pd = yaux1;
vet_erroX = x_ipV - x_pd;
vet_erroY = y_ipV - y_pd;
Mod_Erro = sqrt(vet_erroX.^2 + vet_erroY.^2); % módulo dos vetores de erro
% Calcula o desvio padrão dos erros
desv_erros = std(Mod_Erro);
% % Calcula o IQR dos erros
% Q1_erro = prctile(Mod_Erro, 25);  % Primeiro quartil (25%)
% Q3_erro = prctile(Mod_Erro, 75);  % Terceiro quartil (75%)
% iqr_erro = Q3 - Q1;

% Calcula o MSE dos erros
mse = mean(Mod_Erro.^2);
desv_mse_erros = std(mse);
rmse=sqrt(mse);
fig_errosRMSE_R=figure('Position', [100, 100, 700, 500]); %20%= 4.7700 //25%= 5.9625// 30%=7.1550  // 35%=8.3475 //40%=9.5400 //45%= 10.7325 // 50%= 11.9250
load '2sR1.mat'
load '4sR1.mat'
load '6sR1.mat'
load '8sR1.mat'
load '10sR1.mat'
load '12sR1.mat'
load '14sR1.mat'
%====Projétill P2 - Tiro 2
load '2sR2.mat'
load '4sR2.mat'
load '6sR2.mat'
load '8sR2.mat'
load '10sR2.mat'
load '12sR2.mat'
load '14sR2.mat'

%title('Devio Padr\~ao do m''{o}dulo do erro  $\tilde{x}$','Interpreter','latex','FontSize',13);
%title('Desvio Padrão do módulo do erro de ', 'Interpreter', 'latex', 'FontSize', 13);
% Exemplo de vetor com valores de SNR em dB
snr_db2_4 = SNR81(40:80); % Substitua com seus valores de SNR
snr_db4_6 = SNR81(80:120);
snr_db6_8 = SNR81(120:160);
snr_db8_10 = SNR81(160:200);
snr_db10_12 = SNR81(200:240);
snr_db12_14 = SNR81(240:280);
% Converter SNR de dB para linear
snr_linear2_4 = 10.^(snr_db2_4 / 10);
snr_linear4_6 = 10.^(snr_db4_6 / 10);
snr_linear6_8 = 10.^(snr_db6_8 / 10);
snr_linear8_10 = 10.^(snr_db8_10 / 10);
snr_linear10_12 = 10.^(snr_db10_12 / 10);
snr_linear12_14 = 10.^(snr_db12_14 / 10);

% Calcular a média dos valores lineares
media_snr_linear2_4 = mean(snr_linear2_4);
media_snr_linear4_6 = mean(snr_linear4_6);
media_snr_linear6_8 = mean(snr_linear6_8);
media_snr_linear8_10 = mean(snr_linear8_10);
media_snr_linear10_12 = mean(snr_linear10_12);
media_snr_linear12_14 = mean(snr_linear12_14);

% Converter a média linear de volta para dB
media_snr_db2_4_ = 10 * log10(media_snr_linear2_4);
media_snr_db4_6_ = 10 * log10(media_snr_linear4_6);
media_snr_db6_8_ = 10 * log10(media_snr_linear6_8);
media_snr_db8_10_ = 10 * log10(media_snr_linear8_10);
media_snr_db10_12_ = 10 * log10(media_snr_linear10_12);
media_snr_db12_14_ = 10 * log10(media_snr_linear12_14);

% figure
% load SNR81_new.mat
% load SNR120_new.mat
% 
% plot(tr2A_81_new,SNR81_new); hold on; 
% plot(tr2A_120_new,SNR120_new)  
% 
% xlabel('intervalo de rastreio (s)','Interpreter','latex','FontSize',16); grid
% ylabel('$\overline{SNR}_{dB}$','Interpreter','latex','FontSize',18);
% title('$\overline{SNR}_{dB}$ para cada intervalo de 2s','Interpreter','latex','FontSize',20); %title('Erro quadr\''{a}tico m\''{e}dio de $\tilde{x}$','Interpreter','latex','FontSize',13);
% legend('P$_1$ ','P$_2$ ','Location','North East','Interpreter','latex','FontSize',12);



% figure
% % Exibir o resultado
% %disp(['A SNR média em dB é: ', num2str(media_snr_db)]);
% subplot(2,1,1)
% plot([2 4 6 8 10 12 14],[ sqrt(mse2R1) sqrt(mse4R1) sqrt(mse6R1) sqrt(mse8R1) sqrt(mse10R1) sqrt(mse12R1) sqrt(mse14R1)/1  ],'k-.*','linewidth',1); hold on; grid %plot([4.7700 5.9625 7.1550 8.3475 9.5400 10.7325 11.9250],[mse20 mse25 mse30 mse35 mse40 mse45 mse50],'k-.o','linewidth',1); grid
% plot([2 4 6 8 10 12 14],[ sqrt(mse2R2) sqrt(mse4R2) sqrt(mse6R2) sqrt(mse8R2) sqrt(mse10R2) sqrt(mse12R2) sqrt(mse14R2)  ],'b-.*','linewidth',1); hold on; 
% legend('Tiro P$_1$ ','Tiro P$_2$','Location','North East','Interpreter','latex','FontSize',16);
% xlabel('intervalo de rastreio (s)','Interpreter','latex','FontSize',16)
% ylabel('erro (m)','Interpreter','latex','FontSize',18);
% title('Erro$_{(absoluto)}$','Interpreter','latex','FontSize',20); %title('Erro quadr\''{a}tico m\''{e}dio de $\tilde{x}$','Interpreter','latex','FontSize',13);
% 
%  subplot(2,1,2)
% %plot([2 4 4 6 6 8 8 10 10 12 12 14],[ media_snr_db2_4 media_snr_db2_4 media_snr_db4_6 media_snr_db4_6 media_snr_db6_8 media_snr_db6_8    media_snr_db8_10 media_snr_db8_10 media_snr_db10_12 media_snr_db10_12 media_snr_db12_14 media_snr_db12_14  ],'k--+','linewidth',1); hold on; 
% load SNRT1.mat
% plot([2 4 6 8 10 12 14],[ media_snr_db2_4    media_snr_db4_6 media_snr_db6_8 media_snr_db8_10 media_snr_db10_12 media_snr_db12_14 media_snr_db12_14  ],'k-.*','linewidth',1); hold on; 
% plot([2 4 6 8 10 12 14],[ media_snr_db2_4_  media_snr_db4_6_  media_snr_db6_8_   media_snr_db8_10_ media_snr_db10_12_ media_snr_db12_14_ media_snr_db12_14_  ],'b-','linewidth',1); hold on; 
% 
% xlabel('intervalo de rastreio (s)','Interpreter','latex','FontSize',16); grid
% ylabel('$\overline{SNR}_{dB}$','Interpreter','latex','FontSize',18);
% title('$\overline{SNR}_{dB}$ para cada intervalo de 2s','Interpreter','latex','FontSize',20); %title('Erro quadr\''{a}tico m\''{e}dio de $\tilde{x}$','Interpreter','latex','FontSize',13);
% legend('P$_2$ Tiro 1','P$_2$ Tiro 2','Location','North East','Interpreter','latex','FontSize',12);
% exportgraphics(fig_errosRMSE_R, 'fig_errosRMSE_R2.pdf', 'ContentType', 'vector');

% [x2,y2,z2,t2] = RK4(x1e(end)-x1e(1), y1e(end)-y1e(1), z1e(1),-Vx(iendek),-Vy(iendek), Vz(1), Vlamb/2,tam);
%[x2, y2, z2, vxr2,vyr2, vzr2, vMr2, t22] = RK4(xE2A(1),yN2A(1),zU2A(1), 1*Vx(1),1*Vy(1), 1*Vz(1), Vlamb/2);
[x2, y2, z2, vxr2,vyr2, vzr2, vMr2, t22] = RK4(xE2A(1),yN2A(1),zU2A(1), vxA(1),vyA(1), vzA(1), Vlamb/2);
% tamJanela = 20; 
% 
% bm = (1/tamJanela)*ones(1,tamJanela);
% x_filtMM = filtfilt(bm,1,xE2A);
% y_filtMM = filtfilt(bm,1,yN2A);
% z_filtMM = filtfilt(bm,1,zU2A);


% %taux = zeros(length(x1e));
tz0 = -z_filtMM(5)/(z_filtMM(10)-z_filtMM(5)); 
xLP = x_filtMM(5) +tz0*(x_filtMM(10)-x_filtMM(5));
yLP = y_filtMM(5) +tz0*(y_filtMM(10)-y_filtMM(5));
zLP = z_filtMM(5) +tz0*(z_filtMM(10)-z_filtMM(5));


% 
% xLP = 10;
% yLP = 10;
% zLP = 2.5;
%Estimativas de Ponto de Lançamento=============
% taux = zeros(length(x1e));
% tz0 = -z1e(5)/(z1e(10)-z1e(5)); 
% xLP = x1e(5) +tz0*(x1e(10)-x1e(5));
% yLP = y1e(5) +tz0*(y1e(10)-y1e(5));
% zLP = z1e(5) +tz0*(z1e(10)-z1e(5));
%======================================================================================================================================================================



% [z_max, index] = max(z1e); % Encontra o valor máximo em z e o índice correspondente
% t_max = t(index); % Usa o índice para encontrar o tempo correspondente
% vz_00 = g*t_max;
% vx_00 = x1e(end)/(2*t_max); % testar com t(end) 2*t_max
% vy_00 = y1e(end)/(2*t_max);
% x00 = x1e(end) - vx_00*2*t_max;
% y00 = y1e(end) - vy_00*2*t_max;
% %z00 = z_max - 0.5*g*(t(end)/2)^2;
% z00 = 0;



% figParametros=figure('Position', [100, 100, 700, 500]);
% subplot(4,1,1)
% plot(tk,Cd_est); hold on; plot(t,Cd_est(1:length(t)));grid
% subplot(4,1,2)
% plot(tk,pz); hold on; plot(tr2A,zU2A); plot(tk(1:length(z1e)),z1e(1:end));grid
% subplot(4,1,3)
% plot(tk,py); hold on; plot(tr2A,yN2A); plot(tk(1:length(z1e)),y1e(1:end));grid
% subplot(4,1,4)
% plot(tk,px); hold on; plot(tr2A,xE2A); plot(tk(1:length(z1e)),x1e(1:end));grid
% 
% fig4D=figure('Position', [100, 100, 700, 500]);
% scatter3(xE2A, yN2A, zU2A, 5, tr2A, 'filled'); hold on; grid
% scatter3(px, py, pz, 5, tk, 'filled'); hold on
% scatter3(x1e, y1e, z1e, 5, tk(1:length(z1e)), 'filled'); hold on
% colormap(jet);  % Utiliza o mapa de cores 'jet' para representar o tempo
% 
% %colorbar;       % Adiciona uma barra de cores para indicar a relação com o tempo
% cb = colorbar;%('northoutside');  % Coloca a barra de cores na parte superior e na horizontal;  % Adiciona uma barra de cores
% cb.Label.Interpreter = 'latex';
% cb.Label.String = 'Tempo (s)';  % Legenda da barra de cores
% % Ajustando a espessura da barra de cores
% cb.Position = cb.Position .* [1.1 1 0.15 1];  % Reduz a espessura pela metade
% cb.Label.FontSize = 14;
% 
% xlabel('x (m)','Interpreter','latex','FontSize',14)
% ylabel('y (m)','Interpreter','latex','FontSize',14)
% zlabel('z (m)','Interpreter','latex','FontSize',14)
% title('Trajet\''orias: Medidas e Estimadas','Interpreter','latex','FontSize',15)
% grid on;
% 
% load Cd_81_tiro_2_real.mat
% load Cd_120_tiro_2_real.mat
% figCd_est_real=figure('Position', [100, 100, 700, 500]);
% subplot(2,1,1)
% plot(t_120,Cd_est_120(1:length(t_120)), 'k'); grid
% xlabel('t (s)','Interpreter','latex','FontSize',16)
% ylabel('$\hat{C}_\mathrm{d}$','Interpreter','latex','FontSize',16)
% title('Coeficiente de Arrasto Estimado para P$_1$','Interpreter','latex','FontSize',16)
% %legend('P$_1$', 'Interpreter','latex','Location','northwest','Interpreter','latex','FontSize',16);
% 
% subplot(2,1,2)
% plot(t_81,Cd_est_81(1:length(t_81)), 'k'); grid
% 
% xlabel('t (s)','Interpreter','latex','FontSize',16)
% ylabel('$\hat{C}_\mathrm{d}$','Interpreter','latex','FontSize',16)
% title('Coeficiente de Arrasto Estimado para P$_2$','Interpreter','latex','FontSize',16)
% %legend('P$_2$', 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',16);
% exportgraphics(figCd_est_real, 'figCd_est_real.pdf', 'ContentType', 'vector');


%Método Levenberg-Marquardt================================================================================================




xm = xE2A;
ym = yN2A;
zm = zU2A;

% xm = x1e;
% ym = y1e;
% zm = z1e;


C_Daux =  Vlamb(end)/2;%C_D(end)

sigma_x = wgn(length(xm), 1, 1); 
sigma_y = wgn(length(xm), 1, 1);
sigma_z = wgn(length(xm), 1, 1);
x0 = [xLP, yLP, zLP, Vx(4), Vy(4), Vz(4), C_Daux]; % Valores iniciais dos parâmetros
%x0 = [0, 0, 0, Vx(1), Vy(1), Vz(1), C_Daux];
% % Defina os limites inferiores e superiores para cada parâmetro
  % ub = [Inf, Inf, 10, Inf, Inf, Inf]; % Limites superiores, ajuste conforme necessário
% Limites dos parâmetros [x_0, y_0, z_0, vx_0, vy_0, vz_0]
% % Variações dos parametros do vetor X_0 (LIMITES)
% delta_x = 100;
% delta_y = 100;
% delta_z  =100;
% delta_vx = 100;
% delta_vy = 100;
% delta_vz = 100;
%  lb = [-delta_x,-delta_y, -delta_z, -Inf, -Vy(1), -Vz(1), -Inf]; % Limite inferior (exemplo: [-10, -10, -10, -5, -5, -5])
%  ub = [Inf, Inf, Inf, Inf, Inf, 100,Inf]; % Limite superior (exemplo: [10, 10, 10, 5, 5, 5])

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
%[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x0) funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z), x0, [], [], options);
[x] = lsqnonlin(@(x0) funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z, C_Daux, tr2A), x0, [], [], options);
  %X0_opt = [x1(1) y1(1) z1(1) vxr(1) vyr(1) vzr(1)];
  X0_LM = x;

%====== Estimativas PL===================================================================================================================================================================================
figPIP=figure('Position', [150, 150, 750, 550]);

subplot(2,1,1)
plot3(xlanc,ylanc,zlanc,'bx','linewidth',2); hold on; 
plot3(xE2A(1),yN2A(1),zU2A(1),'kx','Linewidth',2); hold on;
plot3(x2,y2,z2,'k','Linewidth',1); hold on;
[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(X0_LM(1), X0_LM(2),X0_LM(3), X0_LM(4), X0_LM(5), X0_LM(6),X0_LM(7));
plot3(xPL,yPL,zPL,'b--','linewidth',2); hold on;
%plot3(xE2A(1:amostra_parada),yN2A(1:amostra_parada),zU2A(1:amostra_parada),'g-','Linewidth',1); hold on;
tr2A_ = tr2A - tr2A(1);
scatter3(xE2A(1:amostra_parada), yN2A(1:amostra_parada), zU2A(1:amostra_parada), 3, tr2A_(1:amostra_parada), 'filled'); hold on;
%plot3([xLP x1e(1)],[yLP y1e(1)],[0 z1e(1)],'r','Linewidth',1); hold on;
%plot3([x00 x1e(1)],[y00 y1e(1)],[ z00 z1e(1)],'b','Linewidth',1); hold on;

% plot3(x2,y2,z2,'r','linewidth',2); hold on;
%scatter3(px, py, pz, 5, tk, 'filled'); hold on
%scatter3(x1e, y1e, z1e, 5, tk(1:length(z1e)), 'filled'); hold on
colormap(jet);  % Utiliza o mapa de cores 'jet' para representar o tempo

% Supondo que você queira representar faixas específicas do mapa de cores
hold on;
scatter3(xE2A(amostra_parada),yN2A(amostra_parada),zU2A(amostra_parada), 30, 'g', 'filled');
%scatter3(nan, nan, nan, 100, 'g', 'filled', 'DisplayName', 'Faixa de Cor 2');
%scatter3(nan, nan, nan, 100, 'b', 'filled', 'DisplayName', 'Faixa de Cor 3');
legend('show');

%colorbar;       % Adiciona uma barra de cores para indicar a relação com o tempo
cb = colorbar;%('northoutside');  % Coloca a barra de cores na parte superior e na horizontal;  % Adiciona uma barra de cores
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Tempo (s)';  % Legenda da barra de cores
% Ajustando a espessura da barra de cores
cb.Position = cb.Position .* [1.1 1 0.15 1];  % Reduz a espessura pela metade
cb.Label.FontSize = 14;
plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k.','Linewidth',0.1); hold on;
%plot3(px(1:end),py(1:end),pz(1:end),'b-','Linewidth',1); hold on;
plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r-.','Linewidth',1); hold on;
plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
% centro_x = x_ipV;  % Centro da circunferência em x
% centro_y = y_ipV;    % Centro da circunferência em y
% raio = 30;       % Raio da circunferência
% theta = linspace(0, 2*pi, 100);  % Ângulos para definir a circunferência
% circ_x = centro_x + raio * cos(theta);
% circ_y = centro_y + raio * sin(theta);
% fill(circ_x, circ_y, 'r', 'FaceAlpha', 0.3);
%fill(10*x_TLE,10*y_TLE, 'r', 'FaceAlpha', 0.3);

legend('Medida','Estimada');
legend('Lan\c cadora','Ponto de detec\c{c}\~ao','Medidas do radar ','Fim da detec\c{c}\~ao','Traj. Extrapolada','PI$_\mathrm{E}$ de P$_2$','Location','North East','Interpreter','latex','FontSize',12);grid %'Medidas radar ($\Delta$t = 12s)',
%legend('Modelo $P_1$','Medidas','','Lan\c cadora','Ponto de impacto','Interpreter','latex','FontSize',12);grid 
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14)
zlabel('z (m)','Interpreter','latex','FontSize',14)
title('Trajet\''oria 14s do rastreio ','Interpreter','latex','FontSize',15); 
grid on;
%view(3);
AZIMUTE= 11.0613; ELEVACAO=-80.7600; 
view(ELEVACAO, AZIMUTE);
%exportgraphics(fig5, 'fig_PIP_EKF_REAL_120.pdf', 'ContentType', 'vector');

subplot(2,1,2)
%plot_ellipse(xaux1,yaux1,x1(end),y1(end))    ; hold on
plot(x_pd,y_pd,'rx','linewidth',3);axis equal;
%legend('PIP Previsto com $\Delta t = 40$ \quad seg','Interpreter','latex','FontSize',14);
%grid
hold on; 
plot(x_ipV,y_ipV,'b+','linewidth',3);axis equal;hold on; 
%ellipse(vetor_ip, covarianceMatrixC, confidenceLevelC);axis equal;hold on; 
%ellipse(vetor_ip, covarianceMatrix, confidenceLevel); axis equal; hold on; 

quiver(x_ipV, y_ipV, x_pd-x_ipV, y_pd-y_ipV, 0, 'k','Linewidth',1); hold on;
fill(x_TLE,y_TLE, 'g', 'FaceAlpha', 0.1); axis equal; hold on;% O vetor é desenhado em vermelho
%fill(x_max,y_max, 'r', 'FaceAlpha', 0.1); axis equal; hold on;% O vetor é desenhado em vermelho
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14);
Legend_TLE = num2str(r_TLE_final,'Vetor erro(r = %.2f m)');
%legend('PI$_\mathrm{E}$ ($\Delta t = 12$ seg.)','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',12); %%'95\% ellipse (intervalo de confian\c{c}a)',
legend('PI$_\mathrm{E}$ ','PI$_\mathrm{V}$', Legend_TLE,'C\''{i}rculo do vetor', 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',12);
grid on;
exportgraphics(figPIP, 'figPIP_REAL_14sR2.pdf', 'ContentType', 'vector');

% Garantir que as arrays tenham o mesmo tamanho antes de calcular os resíduos
function res = funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z, C_Daux,tr2A2)
    % Calcula a trajetória com os parâmetros atuais   
    [xc, yc, zc, vcx, vcy, vcz, vcM, tc]= RK4(x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7));  %%unction[x1, y1, z1, vxr,vyr, vzr, vMr, t] = RK4(xo, yo, zo, vxro, vyro, vzro, C_D)

 
% Valor a ser encontrado
valorProcurado = tr2A2(1);

% Calcula a diferença absoluta
diferencas = abs(tc - valorProcurado);

% Encontra o índice do valor mais próximo
[~, indice] = min(diferencas);

% % Valor mais próximo em t
% valorMaisProximo = tc(indice);

xc = xc(indice:end);
yc = yc(indice:end);
zc = zc(indice:end);



    % Ajustar o tamanho das arrays
    tamanhoMaximo = max([length(xm), length(ym), length(zm), length(xc), length(yc), length(zc)  ]);%, length(sigma_x), length(sigma_y), length(sigma_z)]);
    xm = padArrayToLength(xm, tamanhoMaximo);
    ym = padArrayToLength(ym, tamanhoMaximo);
    zm = padArrayToLength(zm, tamanhoMaximo);
    xc = padArrayToLength(xc, tamanhoMaximo);
    yc = padArrayToLength(yc, tamanhoMaximo);
    zc = padArrayToLength(zc, tamanhoMaximo);
    % sigma_x = padArrayToLength(sigma_x, tamanhoMaximo);
    % sigma_y = padArrayToLength(sigma_y, tamanhoMaximo);
    % sigma_z = padArrayToLength(sigma_z, tamanhoMaximo);

    sigma_x = std(sigma_x);
    sigma_y = std(sigma_y);
    sigma_z = std(sigma_z);

     res_x  = (xm - xc) ./ (sqrt(2)*sigma_x);
     res_y  = (ym - yc) ./ (sqrt(2)*sigma_y);
     res_z  = (zm - zc) ./ (sqrt(2)*sigma_z);

     % res_x  = (xm - xc) ;
     % res_y  = (ym - yc) ;
     % res_z  = (zm - zc) ;


    % Combina os resíduos em um único vetor
    res = [res_x; res_y; res_z];

       end


function paddedArray = padArrayToLength(array, targetLength)
    if length(array) < targetLength
        % Verificar se o array é coluna ou linha e ajustar os zeros adequadamente
        if iscolumn(array)
            paddedArray = [array; zeros(targetLength - length(array), 1)];
        else
            paddedArray = [array, zeros(1, targetLength - length(array))];
        end
    else
        paddedArray = array;
    end
end