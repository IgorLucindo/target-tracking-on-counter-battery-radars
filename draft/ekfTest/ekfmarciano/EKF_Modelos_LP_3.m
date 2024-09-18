close all
clear 
clc


%==========================================================================
% Simulações de Monte Carlo
 % Inicia a cronometragem
 
ENS=50;
auxENS = zeros(ENS,1);

delta_Tracking = 2; %  28.5500//  35.95  //4seg. de tracking\\TOTAL 23.85 Atenção: o tempo total de rastreio é 23,8s //10%= 2.3850   //20%= 4.7700 //25%= 5.9625// 30%=7.1550  // 35%=8.3475 //40%=9.5400 //45%= 10.7325 // 50%= 11.9250
amostra_parada = 20*delta_Tracking;
amostra_parada = round(amostra_parada);
%Variáveis auxiliares para ponto de impcto previsto
xaux1 = zeros(ENS,1);
yaux1 = zeros(ENS,1);
zaux1 = zeros(ENS,1);
%Variáveis auxiliares para ponto de impcto estimado
xaux2 = zeros(ENS,1);
yaux2 = zeros(ENS,1);
zaux2 = zeros(ENS,1);
xaux3 = zeros(ENS,1);
yaux3 = zeros(ENS,1);
zaux3 = zeros(ENS,1);
xaux4 = zeros(ENS,1);
yaux4 = zeros(ENS,1);
zaux4 = zeros(ENS,1);
erro=zeros(1000,1);
erro_rms = zeros(ENS,1);
erro_rse = 0;
tempo_processamento = 0;
for ens=1:ENS
    
    try
disp(ens)
%========================Dados Projétil====================================
%======================ICAO Atmosphere=====================================
p = 1.203411;
%h1 = 2.926*10^-5 + (10^-10*zMM2A(1:end-2))*3.28084;
%p = po*exp(-h1.*zMM2A(1:end-2));
%======================================================================
d = 81*10^-3;
S = pi*(d^2)/4; % área
m = 4.5; % massa
g = 9.81;

%MÉTODO DE RUNGE KUTTA - 1º Passo: Criar Trajetória teória para Validar Modelo
hr=0.05;    % h = t(i+1)-t(i)  % tamanho do passo
tam = 1000;%850; % tamanho amostra
N = tam*hr; % deta t

%t = 0:hr:N;                  % tamanho do intervalo de t
% t = t + tr2A(1);
vxr = zeros(1,tam);
vyr = zeros(1,tam);
vzr = zeros(1,tam);
vMr= zeros(1,tam);

theta_0 = 45;%260.6385; 
phi_0 = 60; 
v_0 = 200;
range_0 = 500;
% deltaX = (xE2A(3)-xE2A(1));
% deltaY = (yN2A(3)-yN2A(1));
% deltaZ = (zU2A(3)-zU2A(1));
% ANG1 = deltaZ/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2)));
% theta_0 = asin(ANG1); 
% ANG2 = deltaY/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2))*cos(theta_0));
% phi_0 = -acos(ANG2); 
% v_0 = 200;%veloc(1);
vxr(1) = v_0*cosd(theta_0)*sind(phi_0); % condições iniciais
vyr(1) = v_0*cosd(theta_0)*cosd(phi_0); %vy(1) = -0.268;
vzr(1) = v_0*sind(theta_0);             % condições iniciais
vMr(1)= v_0;
x1 = zeros(1,tam);
y1 = zeros(1,tam);
z1 = zeros(1,tam);
% x1(1) = 1500;     % condições iniciais  -1.4953   -0.2465    1.2056
% y1(1) = 250;      % condições iniciais
% z1(1) = 1200;      % condições iniciais
% 
% x1(1) = range_0*cosd(theta_0)*sind(phi_0); % condições iniciais
% y1(1) = range_0*cosd(theta_0)*cosd(phi_0); %vy(1) = -0.268;
% z1(1) = range_0*sind(theta_0);             % condições iniciais
%[x1(1) y1(1) z1(1)]
x1(1) = 25;     % condições iniciais  -1.4953   -0.2465    1.2056
y1(1) = 10;      % condições iniciais
z1(1) = 5;      % condições iniciais





%===================== Processo de Wiener: Método para adicionar aleatóridade a trajetória ========================
Tw = 0.05; Nw = tam;
dt = Tw/Nw;
dW = zeros(1,Nw); % preallocate arrays ...
W = zeros(1,Nw); % for efficiency
dW(1) = (sqrt(dt)*randn); % first approximation outside the loop ...
W(1) = dW(1); % since W(0) = 0 is not allowed
for jw = 2:Nw
dW(jw) = sqrt(1*dt)*randn; % general increment
W(jw) = W(jw-1) + dW(jw);
end

% ===== Método II Desvio a partir de um valor nominal a priori conhecido Cd* 
C_d0 = 0.15; 
%delta = wgn(tam,1,0);

%
%C_d = C_d0*exp(W(1:tam));
C_d = C_d0*ones(1,Nw);
C_D = p*S*C_d/(2*m);
C_Dx = C_D;
C_Dy = C_D;
C_Dz = C_D;
[x1,y1,z1,vx1,vy1,vz1,vM11,t1] = RK4(x1(1), y1(1), z1(1), vxr(1), vyr(1), vzr(1), C_D); % [x1, y1, z1, vxr,vyr, vzr, vMr, t] = RK4(xo, yo, zo, vxro, vyro, vzro, C_D)


%==================== Equações Diferenciais================================

% F_vvx = @(C_Dx,vMr,vxr) -C_Dx.*vMr.*vxr;       % ax
% F_vvy = @(C_Dy,vMr,vyr) -C_Dy.*vMr.*vyr ;      % ay
% F_vvz = @(C_Dz,vMr,vzr) -C_Dz.*vMr.*vzr - g ;  % az           % função EDO 
% 
% F_vx = @(vxr) vxr; % vx
% F_vy = @(vyr) vyr; % vy
% F_vz = @(vzr) vzr; % vz
% 
% for i=1:(tam-1)                                 % acelerações - velocidades comp X
%     k1 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i));
%     k2 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i)+0.5*k1);
%     k3 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i)+0.25*(k1+k2));
%     k4 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i)-(k2+2*k3));
%     vxr(i+1) = vxr(i) + (1/6)*(k1+4*k3+k4);        
% 
%     k1_ = hr*F_vx(vxr(i));                      % velocidades - posições comp X
%     k2_ = hr*F_vx(vxr(i)+0.5*k1_);
%     k3_ = hr*F_vx(vxr(i)+0.25*(k1_+k2_));
%     k4_ = hr*F_vx(vxr(i)-(k2_+2*k3_));
% 
%     x1(i+1) = x1(i) + (1/6)*(k1_+4*k3_+k4_);        
% %==========================================================================
% 
%     Q1 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i));             % acelerações - velocidades comp Y
%     Q2 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i)+0.5*Q1);
%     Q3 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i)+0.25*(Q1+Q2));
%     Q4 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i)-(Q2+2*Q3));
%     vyr(i+1) = (vyr(i) + (1/6)*(Q1+4*Q3+Q4));      
% 
%     Q1_ = hr*F_vy(vyr(i));                            % velocidades - posições comp Y
%     Q2_ = hr*F_vy(vyr(i)+0.5*Q1_);
%     Q3_ = hr*F_vy(vyr(i)+0.25*(Q1_+Q2_));
%     Q4_ = hr*F_vy(vyr(i)-(Q2_+2*Q3_));
% 
%     y1(i+1) = y1(i) + (1/6)*(Q1_+4*Q3_+Q4_);     
% %==========================================================================
% 
%     J1 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i));               % acelerações - velocidades comp Z
%     J2 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i)+0.5*J1);
%     J3 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i)+0.25*(J1+J2));
%     J4 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i)-(J2+2*J3));
%     vzr(i+1) = vzr(i) + (1/6)*(J1+4*J3+J4);    
% 
%     J1_ = hr*F_vz(vzr(i));                              % velocidades - posições comp Z
%     J2_ = hr*F_vz(vzr(i)+0.5*J1_);
%     J3_ = hr*F_vz(vzr(i)+0.25*(J1_+J2_));
%     J4_ = hr*F_vz(vzr(i)-(J2_+2*J3_));
% 
%     z1(i+1) = z1(i) + (1/6)*(J1_+4*J3_+J4_);
% 
% 
% vMr(i+1) = sqrt((vxr(i+1)^2)+(vyr(i+1)^2)+(vzr(i+1)^2)); 
% 
%     if z1(i)<0
%         z1(i) = z1(i-1);
%         y1(i) = y1(i-1);
%         x1(i) = x1(i-1);
%         iend = i-1;
%         break
%     end
% 
% end
% % vm = sqrt((vx.^2)+(vy.^2)+(vz.^2));
% x1 = x1(1:iend);    y1 = y1(1:iend);  z1 = z1(1:iend);  t = t(1:iend);
%==========================================================================
delta_detec = 5; %  28.5500//  35.95  //4seg. de tracking\\TOTAL 23.85 Atenção: o tempo total de rastreio é 23,8s //10%= 2.3850   //20%= 4.7700 //25%= 5.9625// 30%=7.1550  // 35%=8.3475 //40%=9.5400 //45%= 10.7325 // 50%= 11.9250
amostra_inicio = 20*delta_detec;
xMM2A = x1(amostra_inicio:end)';
yMM2A = y1(amostra_inicio:end)';
zMM2A = z1(amostra_inicio:end)';
t=t1(amostra_inicio:end);
tr2A = t';


% xMM2A = x1';
% yMM2A = y1';
% zMM2A = z1';
% t=t(1:end);
% tr2A = t';
% Cálculo das Velocidades =================================================
% Método Comando diff
dx_dt = diff(xMM2A)./diff(tr2A);
dy_dt = diff(yMM2A)./diff(tr2A);
dz_dt = diff(zMM2A)./diff(tr2A);
% dx_dt = diff(xE2A)./diff(tr2A);
% dy_dt = diff(yN2A)./diff(tr2A);
% dz_dt = diff(zU2A)./diff(tr2A);
vxA = dx_dt;
vyA = dy_dt;
vzA = dz_dt; 
v = sqrt((vxA.^2)+(vyA.^2)+(vzA.^2));

v_ = [vxA';vyA';vzA'];


% Cálculo das Acelerações =================================================
dvx_dt = diff(vxA)./diff(tr2A(1:end-1));
dvy_dt = diff(vyA)./diff(tr2A(1:end-1));
dvz_dt = diff(vzA)./diff(tr2A(1:end-1));

ax = dvx_dt; 
ay = dvy_dt; 
az = dvz_dt; 
a = diff(v)./diff(tr2A(1:end-1)); % amodC = diff(vmodC)./diff(tr(1:end-1)');
%ax +ay +az;
a_ = [ax';ay';az'];
% Filtro de Kalman EKF=====================================================

wx = wgn(length(xMM2A), 1, 10);
wy = wgn(length(xMM2A), 1, 10);
wz = wgn(length(xMM2A), 1, 10);

wx_ = wgn(length(xMM2A), 1, 10);
wy_ = wgn(length(xMM2A), 1, 10);
wz_ = wgn(length(xMM2A), 1, 10);

% wx_ = 2*randn(length(xMM2A), 1);
% wy_ = 2*randn(length(xMM2A), 1);
% wz_ = 2*randn(length(xMM2A), 1);
% 


% wx = zeros(length(xMM2A),1);
% wy = zeros(length(xMM2A),1);
% wz = zeros(length(xMM2A),1);

xE2A = xMM2A + wx_;
yN2A = yMM2A + wy_;
zU2A = zMM2A + wz_;
% xE2A = xMM2A;
% yN2A = yMM2A;
% zU2A = zMM2A;

T = 0.05;                    % Período de Amostragem Radar (T =50ms)
%t =tr2A(1:end-2);            %0: T: total_time; % vetor de tempo
n_estados = 7;               % Quantidade de estados / 7 estados e 3 medidas
po =  [xE2A(1); yN2A(1); zU2A(1)]; %  posição inicial do projétil
v_ = v_(:,1);          % velocidade inicial do projétil
%v_ = v_(:,1:end-1);
%v_=[10;10;10];
tam=length(xE2A);

% Inicialização das variáveis==============================================
%alpha = zeros(1,tam);        % parâmetro de arrasto - alpha = S*Cd/m (m^2/kg)
alpha(1) = 10^-4;%S*C_d0/m;% 5.7256*10^-4;%3.9000e-04;% 10^-6;           % inicialização do parâmetro de arrasto
%lamb = zeros(1,tam);
p0 = 1.203411;                % densidade local
%p0 = 0.002378;
%h = 1.0361*10^-4;    %3.158*10^-5 (1/feet) = 1.0361*10^-4 (1/m) ---- ft = m x 3.28084 
%h = 2.926*10^-5 + (10^-10*xk_k(3,k))*3.28084;

%lamb(1) =  1*2.7483e-04;
%lamb(1) =0;
% m = 4.036
% d = 81.03e-3;
% S = pi*(d^2)/4;
xk_k = zeros(n_estados,tam);   % vetor de estados x
f_xk = zeros(n_estados,tam);   % vetor de estados x' (derivada do vetor x)

% px = zeros(tam,1);
% py = zeros(tam,1);
% pz = zeros(tam,1);
% Vx = zeros(tam,1);
% Vy = zeros(tam,1);
% Vz = zeros(tam,1);
% Valpha = zeros(tam,1);
% Cd_est = zeros(tam,1);
% 
% sig2x = wgn(1,tam,5);
% sig2y = wgn(1,tam,5);
% sig2z = wgn(1,tam,5);
% sigxy = wgn(1,tam,5);
% sigxz = wgn(1,tam,5);
% sigyz = wgn(1,tam,5);
% %
% Ruídos do Processo e das medidas
vv =  (10^-7)*wgn(4, tam, 1);  % Ruído do processo v 0.00001
w  = [wx wy wz]';%wgn(3, tam, 1); % Ruído das medidas w

sig2x = 100*ones(1,tam);
sig2y = sig2x;
sig2z = sig2x;
sigxy = 1*rand(1,tam);
sigxz = 1*rand(1,tam);
sigyz = 1*rand(1,tam);
 %load cov_medidas

%R = zeros(3,3); % Matriz de Covariância do erro das medidas (inclui os erros da transformação de coordenadas e ruídos/incertezas dos sensores de medição AER)

 R_o = [sig2x(1) sigxy(1) sigxz(1); 
        sigxy(1) sig2y(1) sigyz(1);
        sigxz(1) sigyz(1) sig2z(1)]; 

%R_o = eye(3)*sig2x(1);
% =================== EKF =================================================
% Inicialização do vetor de estados, 7 estados: x(0|0) = [x(0) y(0) z(0) vx(0) vy(0) vz(0) CD(0)]
%v_(:,1) = [100;100;100];
%h = 2.926*10^-5 + (10^-10*po(3,1))*3.28084;
h = 1.0361*10^-4;
lamb(1) =alpha(1)*p0*exp(-h*po(3,1));
xk_k(:,1)  = [po; v_(:,1) ; lamb(1)];    % Definir estimativa inicial do estado

sig2v= 50^2;       %
%sig2v = var(VelDist)
   
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
sig2lamb = 10^-4;

% Inicialização da Matriz de covariância de estimativa do estado Re_o(0|0) ou P (0|0) - Definir covariância de erro inicial 
P_k_k = [R_o zeros(3,3) zeros(3,1); %(7x7)
      zeros(3,3) sig2v*eye(3,3) zeros(3,1);
      zeros(1,3) zeros(1,3) sig2lamb];

 
%P_k_k = diag([10^6 10^6 10^6 10^6 10^6 10^6 12.86*(10^-14)*exp(-7.38*(10^-5)*xk_k(3,1))]);



% Vetor das medidas reais==================================================
z_med = [xE2A(1:tam) yN2A(1:tam) zU2A(1:tam)]'; % Medidas
H=[eye(3) zeros(3,4)];  % -  Dim:(3x7) Matriz de saída H
% Matriz de Ganho T  -  Dim: (7 x 4)
% T7 = [((T^2)/2)*eye(3)    zeros(3,1);
%        T*eye(3)         zeros(3,1); 
%        zeros(1,3)         1       ];

T7 = [((T^2)/2)*eye(3)    zeros(3,1);
       T*eye(3)         zeros(3,1); 
       zeros(1,3)         1       ];


% Matriz F de transição de estados=========================================
%==========================================================================

%w = zeros(3,tam);%wgn(3, tam, 1, 'dBm');
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
%paux =zeros(tam,1);

for k = 1:tam
  
% Predição do Estado
V_ = [xk_k(4,k) xk_k(5,k) xk_k(6,k)]; % posições das velocidades vx vy e vz
V = norm(V_);
%V = sqrt((xk_k(4,k-1)^2)+(xk_k(5,k-1)^2)+(xk_k(6,k-1)^2));
% For Densidade do ar nível do mar Army Standard Metro:
%                             po = 0.0751265 lbs./cubic foot = 1.203411 kg/m3
%                                  1 kg/m3 =   0.06242796 lb/ft3
%                             ///  1 lb/ft3 = 16.02 kg/m3
VV(k)=V;
%                             h = 3.158*10^-5 (1/feet) = 1.0361*10^-4 (1/m)
%g = 9.80665;                  % aceleração da gravidade m/s^2                % densidade local
%p0 = 0.002378;
p0 = 1.203411;                % densidade local
%p0 = 0.002378;
h = 1.0361*10^-4;    %3.158*10^-5 (1/feet) = 1.0361*10^-4 (1/m) ---- ft = m x 3.28084 
%h = 2.926*10^-5 + (10^-10*xk_k(3,k))*3.28084;
p = p0*exp(-h*xk_k(3,k));      % modelo exponencial da densidade do ar
paux(k) = p;

    
R = [sig2x(k) sigxy(k) sigxz(k);
     sigxy(k) sig2y(k) sigyz(k);
     sigxz(k) sigyz(k) sig2z(k)];

%R = cov(w');

% Matriz de Covariância do Processo
% Q_k = [sqrt(sig2x(k-1))        0          0        0;
%       0           sqrt(sig2y(k-1))        0        0;
%       0              0         sqrt(sig2z(k-1))    0;
%       0              0            0     sqrt(sig2alpha)];



% Q_k = [var(xk_k(1,:))        0          0        0; %10^-8
%       0           var(xk_k(2,:))         0        0;
%       0              0         var(xk_k(3,:))     0;
%       0              0            0     var(alpha)]; %5*10^-8

 

% q_x =1*10^1;  %0.9*10^3;
% q_y =1*10^1;
% q_z =2*10^1;
% q_lamb =10^-3; %(2.2*10^-3)*exp(-h*xk_k(3,k-1));

q_x =10^-3;
q_y =q_x;
q_z =q_x;
%q_lamb =(2.2*10^-3)*exp(-h*xk_k(3,k));% 10^-6
%q_lamb = 7.25*(10^-11)*exp(-0.00005*xk_k(3,k));

q_lamb = 10^-12;%  10^-3

%q_lamb =  (10^-3)*exp(-xk_k(3,k));%10^-6; %  1.5600e-04;   1*10^-4

Q_k = [q_x        0          0        0; %10^-8
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
ee =0;
%===== Elementos da Matriz Jacobiana===================================================================================

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

F_=[f43 f44 f45 f46 f47;
    f53 f54 f55 f56 f57;
    f63 f64 f65 f66 f67;
    f73 f74 f75 f76 f77];

%F_= zeros(4,5);

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

%xk_k(:,k) = xk_k(:,k) + f_xk*T +  F*f_xk*(T^2)/2 + T7*vv(:,k);
%xkm1_k = xk_k(:,k) + f_xk*T + F*f_xk*(T^2)/2;
% xkm1_k = xk_k(:,k-1) + f_xk*T + F*f_xk*(T^2)/2 ;
%  xkm1_k = xk_k(:,k-1) + F*T ;
  

  fi_k = eye(7) + F*T;

% if k==40
% F_open = F;
% fi_k_open = fi_k;
% f_xk_open = f_xk;
% xk_k_open = xk_k(:,k);
% end





P_km1_k = fi_k*P_k_k*fi_k'+Q; 
%P_km1_k = fi_k*(P_k_k+Q)*fi_k'; 
xkm1_k = xk_k(:,k) + f_xk*T + F*f_xk*(T^2)/2+ T7*vv(:,k);
%xkm1_k = xk_k(:,k) + f_xk*T;

%P_km1_k = F*P_k_k*F'+Q; 
%P_km1_k = F*(P_k_k+Q)*F';
% Observações
%z_km1 = H*xk_k(:,k) + w(:,k-1); % z_km1(:,k) = H*xkm1_k;
z_km1 = H*xkm1_k ; %w(:,k); % z_km1(:,k) = H*xkm1_k; 
% Atualização medidas

if k>=amostra_parada
    K=zeros(7,3); 
    xkm1_km1 = xkm1_k; 
    
 else   
   K = P_km1_k*H'/(H*P_km1_k*H' + R);    % K = P_km1_k*H'*inv(H*P_km1_k*H' + R);               % Cálculo do ganho de Kalman
   xkm1_km1 = xkm1_k + K*(z_med(:,k) - z_km1); 
 end

% if xkm1_km1(3)<0
    % k_ext=k;
    % break;
    %end

    





%  if k>=849
%  K=zeros(7,3);
%  %xkm1_km1 = xkm1_k;
% else
 % K = P_km1_k*H'/(H*P_km1_k*H' + R);    % K = P_km1_k*H'*inv(H*P_km1_k*H' + R);               % Cálculo do ganho de Kalman  
 %end
 %xkm1_km1 = xkm1_k + K*(z_med(:,k) - z_km1);    % Atualiza os estados estimados
%xk_k(:,k) = xk_k(:,k) + K*(z_med(:,k) - z_km1);

%P_km1_k = (eye(n_estados)-K*H)*P_km1_k; 
%P_km1_km1 = (eye(n_estados)-K*H)*P_km1_k;          % Atualiza a matriz de Covariancia do erro
P_km1_km1 = (eye(n_estados)-K*H)*P_km1_k*(eye(n_estados)-K*H)'+K*R*K'; 
xk_k(:,k+1) = xkm1_km1;
P_k_k = P_km1_km1;

% px(k)  = xk_k(1,k);
% py(k)  = xk_k(2,k);
% pz(k)  = xk_k(3,k);
% Vx(k)  = xk_k(4,k);
% Vy(k)  = xk_k(5,k);
% Vz(k)  = xk_k(6,k);
% Valpha(k)  = xk_k(7,k);
tk(k) = t(1) +k*T;
px(k)  = xkm1_km1(1);
py(k)  = xkm1_km1(2);
pz(k)  = xkm1_km1(3);
Vx(k)  = xkm1_km1(4);
Vy(k)  = xkm1_km1(5);
Vz(k)  = xkm1_km1(6);
Valpha(k)  = xkm1_km1(7);
Cd_est(k) = Valpha(k)*m/(S.*paux(k));

if pz(k)<0
    k_ext = k-1;
    px  = px(1:k_ext);
py  = py(1:k_ext) ;
pz  = pz(1:k_ext);
Vx  =Vx(1:k_ext);
Vy  = Vy(1:k_ext);
Vz  = Vz(1:k_ext);
Valpha  = Valpha(1:k_ext);
Cd_est = Cd_est(1:k_ext);
tk = tk(1:k_ext);
    break
end



% px(k)  = xkm1_k(1);
% py(k)  = xkm1_k(2);
% pz(k)  = xkm1_k(3);
% Vx(k)  = xkm1_k(4);
% Vy(k)  = xkm1_k(5);
% Vz(k)  = xkm1_k(6);
% Valpha(k)  = xkm1_k(7);
%alpha(k) = xk_k(7);
end

% x1e = nonzeros(x1e);
% y1e = nonzeros(y1e);
% z1e = nonzeros(z1e);
%Fim do Loop Rung Kutta
% 
x1e = px;
y1e = py;
z1e = pz;

xaux1(ens) = x1e(end);
yaux1(ens) = y1e(end);
zaux1(ens) = z1e(end);

xaux2(ens) = px(end);
yaux2(ens) = py(end);
zaux2(ens) = pz(end);

tic;
% Estimativas de Ponto de Lançamento=============
taux = zeros(length(x1e));
tz0 = -z1e(5)/(z1e(10)-z1e(5)); 
xLP = x1e(5) +tz0*(x1e(10)-x1e(5));
yLP = y1e(5) +tz0*(y1e(10)-y1e(5));
zLP = z1e(5) +tz0*(z1e(10)-z1e(5));
% 
% tz0 = -z1e(1)/(z1e(2)-z1e(1)); 
% xLP = x1e(1) +tz0*(x1e(2)-x1e(1));
% yLP = y1e(1) +tz0*(y1e(2)-y1e(1));
% zLP = z1e(1) +tz0*(z1e(2)-z1e(1));

% xLP = 30;     % condições iniciais  -1.4953   -0.2465    1.2056
% yLP = 5;      % condições iniciais
% zLP  = 2;      % condições iniciais


% xm = xE2A';
% ym = yN2A';
% zm = zU2A';
% 
% C_Daux =  Valpha(end)/2;%C_D(end)
% 
% sigma_x = wgn(length(xm), 1, 1); 
% sigma_y = wgn(length(xm), 1, 1);
% sigma_z = wgn(length(xm), 1, 1);
% %x0 = [xLP, yLP, zLP, Vx(1), Vy(1), Vz(1) C_Daux]; % Valores iniciais dos parâmetros
% x0 = [0, 0, 0, Vx(1), Vy(1), Vz(1), C_Daux];
% % % Defina os limites inferiores e superiores para cada parâmetro
%   % lb = [-Inf, -Inf, 0, -Inf, -Inf, -Inf]; % Limite inferior para z0 definido como 0
%   % ub = [Inf, Inf, 10, Inf, Inf, Inf]; % Limites superiores, ajuste conforme necessário
% % Limites dos parâmetros [x_0, y_0, z_0, vx_0, vy_0, vz_0]
% % % Variações dos parametros do vetor X_0 (LIMITES)
% 
% 

%Método Levenberg-Marquardt================================================================================================

xm = xE2A';
ym = yN2A';
zm = zU2A';

% xm = x1e;
% ym = y1e;
% zm = z1e;


C_Daux =  Valpha(end)/2;%C_D(end)

sigma_x = wgn(length(xm), 1, 1); 
sigma_y = wgn(length(xm), 1, 1);
sigma_z = wgn(length(xm), 1, 1);
x0 = [xLP, yLP, zLP, Vx(1), Vy(1), Vz(1), C_Daux]; % Valores iniciais dos parâmetros
%x0 = [0, 0, 0, Vx(1), Vy(1), Vz(1), C_Daux];
% % Defina os limites inferiores e superiores para cada parâmetro
   lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]; % Limite inferior para z0 definido como 0
  % ub = [Inf, Inf, 10, Inf, Inf, Inf]; % Limites superiores, ajuste conforme necessário
% Limites dos parâmetros [x_0, y_0, z_0, vx_0, vy_0, vz_0]
% % Variações dos parametros do vetor X_0 (LIMITES)
delta_x = 100;
delta_y = 100;
delta_z  =100;
delta_vx = 100;
delta_vy = 100;
delta_vz = 100;
 %lb = [-delta_x,-delta_y, -delta_z, iInf, -Vy(1), -Vz(1), -Inf]; % Limite inferior (exemplo: [-10, -10, -10, -5, -5, -5])
 ub = [Inf, Inf, Inf, Inf, Inf, 100,Inf]; % Limite superior (exemplo: [10, 10, 10, 5, 5, 5])

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
%[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x0) funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z), x0, [], [], options);
[x] = lsqnonlin(@(x0) funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z, C_Daux, tr2A), x0, [], [], options);
  X0_opt = [x1(1) y1(1) z1(1) vxr(1) vyr(1) vzr(1)];
  X0_LM = x;

% 


[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(X0_LM(1), X0_LM(2),X0_LM(3), X0_LM(4), X0_LM(5), X0_LM(6),X0_LM(7));

xaux3(ens) = xPL(end);
yaux3(ens) = yPL(end);
zaux3(ens) = zPL(end);

xaux4(ens) = X0_LM(1);
yaux4(ens) = X0_LM(2);
zaux4(ens) = X0_LM(3);



erro = abs((C_d(1:length(Cd_est))'-Cd_est')./C_d(1:length(Cd_est))')*100 ; %+erro(1:length(Cd_est))
%erro_rse=sqrt((C_d(1:length(Cd_est))' - Cd_est').^2) + erro_rse; %sqrt(mean((observed - predicted).^2));
%erro_aux = erro+erro_aux
%erroX = abs((xMM2A(1:length(x1))-x1')./xMM2A(1:length(x1)))*100 + erroX;

catch ME

   % Código para lidar com o erro
        disp(['Erro na simulação ' num2str(ens) ': ' ME.message]);
    end
tempo_processamento = toc + tempo_processamento;
end
%Fim do Loop Monte Carlo
%erro = erro/ENS;
erro_rms=sum(erro_rms)/ENS;
tempo_processamento = tempo_processamento/ENS
fprintf('Tempo de processamento: %f segundos\n', tempo_processamento);
%=====================Open Loop Kalman Filter===================================================================================
% xk_k_o = zeros(n_estados,tam-40);   % vetor de estados x
% xk_k_o(:,40) = xk_k_open;
% for kk=40:tam-40
% xkm1_kopen = xk_k_o(:,kk) + f_xk_open*T + F_open*f_xk_open*(T^2)/2 + T7*vv(:,k);
% xk_k_o(:,kk+1) = xkm1_kopen ;
% end
% Pontos de Impacto EKF completo
 matriz_aux2 = [xaux2 yaux2];
[idx,cord_ponto_central] = kmeans(matriz_aux2,1);

xaux3 = nonzeros(xaux3);
yaux3 =nonzeros(yaux3);
xaux4 = nonzeros(xaux4);
yaux4 = nonzeros(yaux4);


matriz_aux3 = [xaux3 yaux3];
[idx3,cord_ponto_central3] = kmeans(matriz_aux3,1);

matriz_aux4 = [xaux4 yaux4];
[idx4,cord_ponto_central4] = kmeans(matriz_aux4,1);

%======================== Plotando TLE  ===================================
% xaux1 = randn(100, 1);
% yaux1 = 2*xaux1 + 0.5*randn(100, 1);
% 
x_ipV = x1(end);
y_ipV = y1(end);
x_lpV = x1(1);
y_lpV = y1(1);
%================ PI EKF=======
x_ip = cord_ponto_central(1,1); 
y_ip = cord_ponto_central(1,2);
x_pd_k = (xaux1);
y_pd_k =(yaux1);
% var_x_pd = var(x_pd_k);
% cov_x_pd_y_pd = sum((x_pd_k - mean(x_pd_k)) .* (y_pd_k - mean( y_pd_k))) / (length(x_pd_k) - 1);
% var_y_pd =  var(y_pd_k);
vetor_ip = [x_ip, y_ip];  % Vetor de médias %meanVector = [mu_x,mu_y];
erro_x =  (x_ip - x_pd_k);
erro_y =  (y_ip - y_pd_k);
%covarianceMatrix = cov(x_pd,y_pd); 
covarianceMatrix = cov(erro_x,erro_y); 
%covarianceMatrix = [var_x_pd, cov_x_pd_y_pd; cov_x_pd_y_pd, var_y_pd];  % Matriz de covariância
%covarianceMatrix = cov(x_pd,y_pd); 
 [V1, D1] = eig(covarianceMatrix);
 lambda_max = max(diag(D1));% 
 covarianceMatrixC = [lambda_max, 0; 0, lambda_max];
% Nível de confiança (por exemplo, 0.95 para 95% de confiança)
confidenceLevel = 0.991;
confidenceLevelC = 0.99;

%===============PI ML ===============
x_ip3 = cord_ponto_central3(1,1); 
y_ip3 = cord_ponto_central3(1,2);
x_pd_LM = (xaux3);
y_pd_LM =(yaux3);
vetor_ip_LM = [x_ip3, y_ip3];  % Vetor de médias %meanVector = [mu_x,mu_y];
erro_x_3 =  (x_ip3 - x_pd_LM);
erro_y_3 =  (y_ip3 - y_pd_LM);
%covarianceMatrix = cov(x_pd,y_pd); 
covarianceMatrix_LM = cov(erro_x_3,erro_y_3); 
%covarianceMatrix = [var_x_pd, cov_x_pd_y_pd; cov_x_pd_y_pd, var_y_pd];  % Matriz de covariância
%covarianceMatrix = cov(x_pd,y_pd); 
 [V1_LM, D1_LM] = eig(covarianceMatrix_LM);
 lambda_max_LM = max(diag(D1_LM));% 
 covarianceMatrixC_LM = [lambda_max_LM, 0; 0, lambda_max_LM];
% Nível de confiança (por exemplo, 0.95 para 95% de confiança)


%===============PL ML ===============
x_ip4 = cord_ponto_central4(1,1); 
y_ip4 = cord_ponto_central4(1,2);
x_pd_LM_L = (xaux4);
y_pd_LM_L =(yaux4);
vetor_ip_LM_L = [x_ip4, y_ip4];  % Vetor de médias %meanVector = [mu_x,mu_y];
erro_x_4 =  (x_ip4 - x_pd_LM_L);
erro_y_4 =  (y_ip4 - y_pd_LM_L);
%covarianceMatrix = cov(x_pd,y_pd); 
covarianceMatrix_LM_L = cov(erro_x_4,erro_y_4); 
 [V1_LM_L, D1_LM_L] = eig(covarianceMatrix_LM_L);
 lambda_max_LM_L = max(diag(D1_LM_L));% 
 covarianceMatrixC_LM_L = [lambda_max_LM_L, 0; 0, lambda_max_LM_L];
% Nível de confiança (por exemplo, 0.95 para 95% de confiança)

% x_pd = rmmissing(xaux1);
% y_pd = rmmissing(yaux1);

% x_pd = nonzeros(xaux1);
% y_pd = nonzeros(yaux1);

% x_ip = mean(x_pd);
% y_ip = mean(y_pd);

%mu_x ==x_ip;
%mu_y ==y_ip;


% % Calcula os autovalores e autovetores da matriz de covariância
% [V, D] = eig(covarianceMatrix);
% [Vc, Dc] = eig(covarianceMatrixC);
% 
% % Calcula os comprimentos dos semi-eixos
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

%======================== Cáculo do TLE PI EKF ==================================
r_TLE_s = sqrt((x_pd_k-x_ipV).^2+(y_pd_k-y_ipV).^2);
r_TLE_s_ord = sort(r_TLE_s);
TLE_porcent = 85;
%TLE_amostra_porcent = (TLE_porcent*ens)/100;
TLE_amostra_porcent = (TLE_porcent*length(x_pd_k))/100;TLE_amostra_porcent=round(TLE_amostra_porcent);
r_TLE_final = r_TLE_s_ord(TLE_amostra_porcent );

theta_TLE = linspace(0,2*pi);
x_TLE = r_TLE_final*cos(theta_TLE) + x_ipV(end);
y_TLE = r_TLE_final*sin(theta_TLE) + y_ipV(end);


%======================== Cáculo do TLE PI ML ==================================
r_TLE_ML = sqrt((x_pd_LM-x_ipV).^2+(y_pd_LM-y_ipV).^2);
r_TLE_ML_ord = sort(r_TLE_ML);
TLE_porcent = 85;
%TLE_amostra_porcent = (TLE_porcent*ens)/100;
TLE_amostra_porcent_ML = (TLE_porcent*length(x_pd_LM))/100; TLE_amostra_porcent_ML=round(TLE_amostra_porcent_ML);
r_TLE_final_ML = r_TLE_ML_ord(TLE_amostra_porcent_ML);

theta_TLE = linspace(0,2*pi);
x_TLE_ML = r_TLE_final_ML*cos(theta_TLE) + x_ipV(end);
y_TLE_ML = r_TLE_final_ML*sin(theta_TLE) + y_ipV(end);

%======================== Cáculo do TLE PL ML ==================================
r_TLE_ML_L = sqrt((x_pd_LM_L-x_lpV).^2+(y_pd_LM_L-y_lpV).^2);
r_TLE_ML_L_ord = sort(r_TLE_ML_L);
TLE_porcent = 85;
%TLE_amostra_porcent = (TLE_porcent*ens)/100;
TLE_amostra_porcent_ML_L = (TLE_porcent*length(x_pd_LM_L))/100; TLE_amostra_porcent_ML_L=round(TLE_amostra_porcent_ML_L);
r_TLE_final_ML_L = r_TLE_ML_L_ord(TLE_amostra_porcent_ML_L);

theta_TLE = linspace(0,2*pi);
x_TLE_ML_L = r_TLE_final_ML_L*cos(theta_TLE) + x_lpV(end);
y_TLE_ML_L = r_TLE_final_ML_L*sin(theta_TLE) + y_lpV(end);


%============
% Calculo dos vetores de erro e medidas de dispersão

% x_pd = xaux1;
% y_pd = yaux1;
% vet_erroX = x_ipV - x_pd;
% vet_erroY = y_ipV - y_pd;
vet_erroX = x_ipV - x_pd_k;
vet_erroY = y_ipV - y_pd_k;
Mod_Erro = sqrt(vet_erroX.^2 + vet_erroY.^2); % módulo dos vetores de erro
% Calcula o desvio padrão dos erros
desv_erros = std(Mod_Erro);
% Calcula o IQR dos erros
% Q1_erro = prctile(Mod_Erro, 25);  % Primeiro quartil (25%)
% Q3_erro = prctile(Mod_Erro, 75);  % Terceiro quartil (75%)
% iqr_erro = Q3 - Q1;

% Calcula o MSE dos erros
mse = mean(Mod_Erro.^2);
desv_mse_erros = std(mse);
rmse=sqrt(mse);



% %Método Levenberg-Marquardt================================================================================================
% % xm = x1e;
% % ym = y1e;
% % zm = z1e;
% 
% xm = xE2A';
% ym = yN2A';
% zm = zU2A';
% 
% 
% C_Daux =  Valpha(end)/2;%C_D(end)
% 
% sigma_x = wgn(length(xm), 1, 1); 
% sigma_y = wgn(length(xm), 1, 1);
% sigma_z = wgn(length(xm), 1, 1);
% x0 = [xLP, yLP, zLP, Vx(1), Vy(1), Vz(1) C_Daux]; % Valores iniciais dos parâmetros
% %x0 = [0, 0, 0, Vx(1), Vy(1), Vz(1), C_Daux];
% % % Defina os limites inferiores e superiores para cada parâmetro
%   % lb = [-Inf, -Inf, 0, -Inf, -Inf, -Inf]; % Limite inferior para z0 definido como 0
%   % ub = [Inf, Inf, 10, Inf, Inf, Inf]; % Limites superiores, ajuste conforme necessário
% % Limites dos parâmetros [x_0, y_0, z_0, vx_0, vy_0, vz_0]
% % % Variações dos parametros do vetor X_0 (LIMITES)
% delta_x = 100;
% delta_y = 100;
% delta_z  =100;
% delta_vx = 100;
% delta_vy = 100;
% delta_vz = 100;
%  lb = [0, 0, 0, Vx(1), Vy(1), Vz(1), C_Daux-C_Daux*0.1]; % Limite inferior (exemplo: [-10, -10, -10, -5, -5, -5])
%  ub = [delta_x, delta_y, delta_z, Vx(1)+delta_vx, Vy(1)+delta_vy, Vz(1)+delta_vz,C_Daux+C_Daux*0.1]; % Limite superior (exemplo: [10, 10, 10, 5, 5, 5])
% 
% options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
% %[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x0) funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z), x0, [], [], options);
% [x] = lsqnonlin(@(x0) funcaoResiduos(x0, xm, ym, zm, sigma_x, sigma_y, sigma_z, C_Daux, tr2A), x0, [], [], options);
%    X0_opt = [x1(1) y1(1) z1(1) vxr(1) vyr(1) vzr(1)]
%   X0_LM = x
%   C_D_real = C_D(end)
%   CD_LM = X0_LM(7)
%   C_D_est = Valpha(end)/2
%   erroCD = C_D(end) - X0_LM(end); %output



fig_BOI = figure('Position', [450, 450, 1000, 800]);
subplot(2,2,1:2)
plot3(xE2A(1),yN2A(1),zU2A(1),'kx','Linewidth',5); hold on;
%plot3([xLP x1e(1)],[yLP y1e(1)],[zLP z1e(1)],'r','Linewidth',1); hold on;
[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(X0_LM(1), X0_LM(2),X0_LM(3), X0_LM(4), X0_LM(5), X0_LM(6),X0_LM(7)); 
scatter3(xLP,yLP,zLP, 50, 'm', 'filled');hold on;
%plot3([x00 x1e(1)],[y00 y1e(1)],[ z00 z1e(1)],'b','Linewidth',1); hold on;
%plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k','Linewidth',1); hold on;
tr2A_ = tr2A - tr2A(1);
scatter3(xE2A(1:amostra_parada), yN2A(1:amostra_parada), zU2A(1:amostra_parada), 10, tr2A_(1:amostra_parada), 'filled'); hold on;
colormap(jet);  % Utiliza o mapa de cores 'jet' para representar o tempo
% Supondo que você queira representar faixas específicas do mapa de cores
hold on;
scatter3(xE2A(amostra_parada),yN2A(amostra_parada),zU2A(amostra_parada), 50, 'green', 'filled');hold on;


%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(xMin, yMin, zMin, vxMin, vyMin, vzMin,Valpha/2); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(25,  15,    2,   93.45,   34.19,  173.14,C_D); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(30,  20,    0,   97.45,   33.19,  183.14,C_D); 
plot3(x1,y1,z1,'g','Linewidth',1); hold on;

%colorbar;       % Adiciona uma barra de cores para indicar a relação com o tempo
cb = colorbar;%('northoutside');  % Coloca a barra de cores na parte superior e na horizontal;  % Adiciona uma barra de cores
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Tempo (s)';  % Legenda da barra de cores
% Ajustando a espessura da barra de cores
cb.Position = cb.Position .* [1.1 1 0.25 1];  % Reduz a espessura pela metade
cb.Label.FontSize = 14;
%plot3(xE2A,yN2A,zU2A,'k','Linewidth',1); hold on;

plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r-.','Linewidth',1); hold on;
plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
plot3(xPL,yPL,zPL,'b--','linewidth',2); hold on;
plot3(xPL(end),yPL(end),zPL(end),'bO','Linewidth',2); hold on;

xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14)
zlabel('z (m)','Interpreter','latex','FontSize',14)
title('Trajet\''oria 3D','Interpreter','latex','FontSize',15); 
legend('Ponto de detec\c{c}\~ao','Estimativa inicial de PL','Medidas do radar','Fim da detec\c{c}\~ao','Trajet\''oria real','Trajet\''oria extrapolada EKF','PI EKF','Trajet\''oria extrapolada ML','PI ML','Location','North East','Interpreter','latex','FontSize',15); grid
AZIMUTE=105.9135; ELEVACAO=21.5202;
 view(AZIMUTE, ELEVACAO);

 subplot(2,2,3)

 plot3(xE2A(1),yN2A(1),zU2A(1),'kx','Linewidth',5); hold on;
%plot3([xLP x1e(1)],[yLP y1e(1)],[zLP z1e(1)],'r','Linewidth',1); hold on;
scatter3(xLP,yLP,zLP, 50, 'm', 'filled');hold on;
%plot3([x00 x1e(1)],[y00 y1e(1)],[ z00 z1e(1)],'b','Linewidth',1); hold on;
%plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k','Linewidth',1); hold on;
tr2A_ = tr2A - tr2A(1);
scatter3(xE2A(1:amostra_parada), yN2A(1:amostra_parada), zU2A(1:amostra_parada), 10, tr2A_(1:amostra_parada), 'filled'); hold on;
colormap(jet);  % Utiliza o mapa de cores 'jet' para representar o tempo
% Supondo que você queira representar faixas específicas do mapa de cores
hold on;
scatter3(xE2A(amostra_parada),yN2A(amostra_parada),zU2A(amostra_parada), 50, 'green', 'filled');hold on;
[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(X0_LM(1), X0_LM(2),X0_LM(3), X0_LM(4), X0_LM(5), X0_LM(6),X0_LM(7)); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(xMin, yMin, zMin, vxMin, vyMin, vzMin,Valpha/2); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(25,  15,    2,   93.45,   34.19,  173.14,C_D); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(30,  20,    0,   97.45,   33.19,  183.14,C_D); 
plot3(x1,y1,z1,'g','Linewidth',1); hold on;

% %colorbar;       % Adiciona uma barra de cores para indicar a relação com o tempo
% cb = colorbar;%('northoutside');  % Coloca a barra de cores na parte superior e na horizontal;  % Adiciona uma barra de cores
% cb.Label.Interpreter = 'latex';
% cb.Label.String = 'Tempo (s)';  % Legenda da barra de cores
% % Ajustando a espessura da barra de cores
% cb.Position = cb.Position .* [1.1 1 0.25 1];  % Reduz a espessura pela metade
% cb.Label.FontSize = 14;
% %plot3(xE2A,yN2A,zU2A,'k','Linewidth',1); hold on;

plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r-.','Linewidth',1); hold on;
plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
plot3(xPL,yPL,zPL,'b--','linewidth',2); hold on;
plot3(xPL(end),yPL(end),zPL(end),'bO','Linewidth',2); hold on;

xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14)
zlabel('z (m)','Interpreter','latex','FontSize',14)
title('Trajet\''oria no plano (x,y) ','Interpreter','latex','FontSize',15); grid
%legend('Ponto de detec\c{c}\~ao','Estimativa inicial de PL','Medidas do radar','Fim da detec\c{c}\~ao','Trajet\''oria real','Trajet\''oria extrapolada EKF','PI EKF','Trajet\''oria extrapolada ML','PI ML','Location','North East','Interpreter','latex','FontSize',12); grid
AZIMUTE=0; ELEVACAO=90;
 view(AZIMUTE, ELEVACAO);

  subplot(2,2,4)

plot3(xE2A(1),yN2A(1),zU2A(1),'kx','Linewidth',5); hold on;
%plot3([xLP x1e(1)],[yLP y1e(1)],[zLP z1e(1)],'r','Linewidth',1); hold on;
scatter3(xLP,yLP,zLP, 50, 'm', 'filled');hold on;
%plot3([x00 x1e(1)],[y00 y1e(1)],[ z00 z1e(1)],'b','Linewidth',1); hold on;
%plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k','Linewidth',1); hold on;
tr2A_ = tr2A - tr2A(1);
scatter3(xE2A(1:amostra_parada), yN2A(1:amostra_parada), zU2A(1:amostra_parada), 10, tr2A_(1:amostra_parada), 'filled'); hold on;
colormap(jet);  % Utiliza o mapa de cores 'jet' para representar o tempo
% Supondo que você queira representar faixas específicas do mapa de cores
hold on;
scatter3(xE2A(amostra_parada),yN2A(amostra_parada),zU2A(amostra_parada), 50, 'green', 'filled');hold on;
[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(X0_LM(1), X0_LM(2),X0_LM(3), X0_LM(4), X0_LM(5), X0_LM(6),X0_LM(7)); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(xMin, yMin, zMin, vxMin, vyMin, vzMin,Valpha/2); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(25,  15,    2,   93.45,   34.19,  173.14,C_D); 
%[xPL,yPL,zPL,vxPL, vyPL, vzPL,vmPL,tPL] = RK4(30,  20,    0,   97.45,   33.19,  183.14,C_D); 
plot3(x1,y1,z1,'g','Linewidth',1); hold on;

%colorbar;       % Adiciona uma barra de cores para indicar a relação com o tempo
cb = colorbar;%('northoutside');  % Coloca a barra de cores na parte superior e na horizontal;  % Adiciona uma barra de cores
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Tempo (s)';  % Legenda da barra de cores
% Ajustando a espessura da barra de cores
cb.Position = cb.Position .* [1.1 1 0.25 1];  % Reduz a espessura pela metade
cb.Label.FontSize = 14;
%plot3(xE2A,yN2A,zU2A,'k','Linewidth',1); hold on;

plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r-.','Linewidth',1); hold on;
plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
plot3(xPL,yPL,zPL,'b--','linewidth',2); hold on;
plot3(xPL(end),yPL(end),zPL(end),'bO','Linewidth',2); hold on;

xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14)
zlabel('z (m)','Interpreter','latex','FontSize',14)
title('Trajet\''oria no plano (x,z) ','Interpreter','latex','FontSize',15); grid
%legend('Ponto de detec\c{c}\~ao','Estimativa inicial de PL','Medidas do radar','Fim da detec\c{c}\~ao','Trajet\''oria real','Trajet\''oria extrapolada EKF','PI EKF','Trajet\''oria extrapolada ML','PI ML','Location','North East','Interpreter','latex','FontSize',12); grid
AZIMUTE=0; ELEVACAO=0;
 view(AZIMUTE, ELEVACAO);
exportgraphics(fig_BOI, 'figPIP_TRAJ1_nominal_ML2.pdf', 'ContentType', 'vector');



  % % %===========================================================================
fig_PIP_TRAJ_ML = figure('Position', [450, 450, 1000, 800]);

subplot(3,1,1)
plot(x_pd_k,y_pd_k,'ro','linewidth',1);daspect([1 1 1]);hold on; 
plot(x_ipV,y_ipV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
ellipse(vetor_ip, covarianceMatrixC, confidenceLevelC,'r'); daspect([1 1 1]);hold on; 
ellipse(vetor_ip, covarianceMatrix, confidenceLevel,'r'); hold on;
fill(x_TLE,y_TLE, 'g', 'FaceAlpha', 0.1);
title('Pontos de impacto estimados pelo EKF','Interpreter','latex','FontSize',15)
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14);
Legend_TLE = num2str(r_TLE_final,'TLE 85 (r = %.2fm)');
%legend('PI$_\mathrm{E}$ ','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',12); %%'95\% ellipse (intervalo de confian\c{c}a)',
legend('PI$_\mathrm{E}$ ','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE ,'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); %%'95\% ellipse (intervalo de confian\c{c}a)',
%axis([2600 3400 990 1310])
%exportgraphics(fig_PIP_TRAJ, 'figPIP_TRAJ_nominal_1.0.pdf', 'ContentType', 'vector');
%set(gcf,'outerposition',get(0,'screensize'))
%export_fig ('test_full_screen_export_macos','-png','-r300', '-transparent',gcf);

subplot(3,1,2)
plot(x_pd_LM,y_pd_LM,'ko','linewidth',1);daspect([1 1 1]);hold on; 
plot(x_ipV,y_ipV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
ellipse(vetor_ip_LM, covarianceMatrixC_LM, confidenceLevelC,'k'); daspect([1 1 1]);hold on; 
ellipse(vetor_ip_LM, covarianceMatrix_LM, confidenceLevel,'k'); 
fill(x_TLE_ML,y_TLE_ML, 'g', 'FaceAlpha', 0.1);
title('Pontos de impacto estimados pelo crit\''{e}rio ML','Interpreter','latex','FontSize',15)
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14);
Legend_TLE = num2str(r_TLE_final_ML,'TLE 85 (r = %.2fm)');
legend('PI$_\mathrm{E}$ ','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); %%'95\% ellipse (intervalo de confian\c{c}a)',
%axis([2600 3400 990 1310])
%exportgraphics(fig_PIP_TRAJ, 'figPIP_TRAJ_nominal_1.0.pdf', 'ContentType', 'vector');
%set(gcf,'outerposition',get(0,'screensize'))
%export_fig ('test_full_screen_export_macos','-png','-r300', '-transparent',gcf);

subplot(3,1,3)
plot(x_pd_k,y_pd_k,'ro','linewidth',1);daspect([1 1 1]);hold on; 
plot(x_ipV,y_ipV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
ellipse(vetor_ip, covarianceMatrixC, confidenceLevelC,'r'); daspect([1 1 1]);hold on; 
ellipse(vetor_ip, covarianceMatrix, confidenceLevel,'r'); hold on;
%fill(x_TLE,y_TLE, 'g', 'FaceAlpha', 0.1);
plot(x_pd_LM,y_pd_LM,'ko','linewidth',1);daspect([1 1 1]);hold on; 
%plot(x_ipV,y_ipV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
ellipse(vetor_ip_LM, covarianceMatrixC_LM, confidenceLevelC,'k'); daspect([1 1 1]);hold on; 
ellipse(vetor_ip_LM, covarianceMatrix_LM, confidenceLevel,'k'); 
title('Comparativo das estimativas: Estimador ML e EKF','Interpreter','latex','FontSize',15)
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14);
legend('PI$_{\mathrm{E}}(\mathrm{EKF})$ ','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\%  $ \chi^2_2$ (EKF)','Elipse 99\% $\chi^2_2 $ (EKF)','PI$_{\mathrm{E}}$ (ML)','C\''{i}rculo de Alerta 99\%  $ \chi^2_2$ (ML)','Elipse 99\% $\chi^2_2 $ (ML)', 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); %%'95\% ellipse (intervalo de confian\c{c}a)',

exportgraphics(fig_PIP_TRAJ_ML, 'figPIP_TRAJ_nominal_ML2.pdf', 'ContentType', 'vector');



fig_PLP_TRAJ = figure%('Position', [450, 450, 1000, 800]);
%subplot(2,1,1)
plot(x_pd_LM_L,y_pd_LM_L,'ro','linewidth',1);daspect([1 1 1]); hold on; 
plot(x_lpV,y_lpV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
ellipse(vetor_ip_LM_L, covarianceMatrixC_LM_L, confidenceLevelC,'k'); daspect([1 1 1]);hold on; 
ellipse(vetor_ip_LM_L, covarianceMatrix_LM_L, confidenceLevel,'k'); 
fill(x_TLE_ML_L,y_TLE_ML_L, 'g', 'FaceAlpha', 0.1);
title('Pontos de lan\c{c}amento estimados pelo crit\''{e}rio ML','Interpreter','latex','FontSize',15)
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14);
Legend_TLE = num2str(r_TLE_final_ML_L,'TLE 85 (r = %.2fm)');
legend('PI$_\mathrm{E}$ ','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); %%'95\% ellipse (intervalo de confian\c{c}a)',
%axis([2600 3400 990 1310])
exportgraphics(fig_PLP_TRAJ, 'figPLP_TRAJ_nominal_4.pdf', 'ContentType', 'vector');
%set(gcf,'outerposition',get(0,'screensize'))
%export_fig ('test_full_screen_export_macos','-png','-r300', '-transparent',gcf);

% subplot(2,1,2)
% plot(x_pd_LM_L,y_pd_LM_L,'ro','linewidth',1);daspect([1 1 1]);hold on; 
% scatter3(xLP,yLP,zLP, 50, 'm', 'filled');hold on;
% plot(x_lpV,y_lpV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
% ellipse(vetor_ip_LM_L, covarianceMatrixC_LM_L, confidenceLevelC,'k'); daspect([1 1 1]);hold on; 
% ellipse(vetor_ip_LM_L, covarianceMatrix_LM_L, confidenceLevel,'k'); 
% %fill(x_TLE,y_TLE, 'g', 'FaceAlpha', 0.1);
% title('Estimativa inicial do PL e estimativas pelo crit\''{e}rio ML','Interpreter','latex','FontSize',15)
% xlabel('x (m)','Interpreter','latex','FontSize',14)
% ylabel('y (m)','Interpreter','latex','FontSize',14);
% Legend_TLE = num2str(r_TLE_final,'TLE 99 (r = %.2fm)');
% legend('PI$_\mathrm{E}$ ','Estimativa inicial de PL','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',12); %%'95\% ellipse (intervalo de confian\c{c}a)',
% %axis([2600 3400 990 1310])
% %exportgraphics(fig_PIP_TRAJ, 'figPIP_TRAJ_nominal_1.0.pdf', 'ContentType', 'vector');
% %set(gcf,'outerposition',get(0,'screensize'))
% %export_fig ('test_full_screen_export_macos','-png','-r300', '-transparent',gcf);

%=====Algorítmos Genéticos==================================================


% % % Variações dos parametros do vetor X_0 (LIMITES)
% % delta_x = 50;
% % delta_y = 50;
% % delta_z = 10;
% % delta_vx = 100;
% % delta_vy = 100;
% % delta_vz = 100;
% % 
% % 
% % % Limites dos parâmetros [x_0, y_0, z_0, vx_0, vy_0, vz_0]
% % lb = [0, 0,0, Vx(1), Vy(1), Vz(1), 0.5*(Valpha(end)/2) ]; % Limite inferior (exemplo: [-10, -10, -10, -5, -5, -5]) 
% % %lb = [0, 0,0, 50, 50, 50, 10^-4 ]; % Limite inferior (exemplo: [-10, -10, -10, -5, -5, -5]) 
% % ub = [delta_x, delta_y, delta_z, Vx(1)+delta_vx, Vy(1)+delta_vy, Vz(1)+delta_vz, 2*(Valpha(end)/2) ]; % Limite superior (exemplo: [10, 10, 10, 5, 5, 5])2*(Valpha(end)/2)
%ub = [delta_x, delta_y, delta_z, 50+delta_vx, 50+delta_vy, 50+delta_vz, 10^-3]; % Limite superior (exemplo: [10, 10, 10, 5, 5, 5])2*(Valpha(end)/2)

C_D_aux2 = C_D(end); %Valpha(end)/2;
%C_D_aux = C_D(end); 
% % Executando o algoritmo genético
% [x_min, fval] = ga(@funcaoObjetivo, 6, [], [], [], [], lb, ub, [], opcoes);
% % Modifique a chamada do algoritmo genético para usar uma função anônima

% 
% 
% %x = ga(fun,2,A,b)
% % Opções do algoritmo genético
% opcoes = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf, 'PopulationSize', 100, 'MaxGenerations',50, 'Display', 'iter');
% %options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
% %[x_min, fval] = ga(@(x_AG) funcaoObjetivo(x_AG, xE2A', yN2A', zU2A',C_D_aux2,tr2A), 7, [], [], [], [], [], [], [], opcoes);
% [x_min,fval,exitflag,output,population,scores] = ga(@(x) funcaoObjetivo(x, xE2A', yN2A', zU2A',C_D_aux2,tr2A), 7, [], [], [], [], lb, ub, [], opcoes);

% %[x_min, fval] = ga(@(x) funcaoObjetivo(x, x1e, y1e, z1e,C_D_aux,tr2A), 7, [], [], [], [], lb, ub, [], opcoes);
% %[x_min, fval] = ga(@(x) funcaoObjetivo(x, x1e, y1e, z1e,C_D_aux,tr2A), 7, [], [], [], [], [], [], [], opcoes);
% 
% % x = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
% % x = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon)
% % x = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options)
% 
% 
% X0_opt = [x1(1) y1(1) z1(1) vxr(1) vyr(1) vzr(1)];
% 
% % Resultado
% disp('Solução ótima:');
% disp(X0_opt);
% disp('Solução encontrada AG:');
% disp(x_min);
% disp('Valor da função objetivo:');
% disp(fval);


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
    tamanhoMaximo = max([length(xm), length(ym), length(zm), length(xc), length(yc), length(zc)]);%, length(sigma_x), length(sigma_y), length(sigma_z)]);
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



    
    
function erro = funcaoObjetivo(x_AG, xm, ym, zm,C_D_aux,tr2A2)
% x é o vetor [x_0, y_0, z_0, vx_0, vy_0, vz_0]
% trajetoriaMedida = [x1e; y1e; z1e];
%trajetoriaMedida = [xE2A'; yN2A'; zU2A'];

% xm = x1e;
% ym = y1e;
% zm = z1e;
% 
% xm = xE2A';
% ym = yN2A';
% zm = zU2A';


% Desempacotar o vetor x
    xo = x_AG(1);
    yo = x_AG(2);
    zo = x_AG(3);
    vxro = x_AG(4);
    vyro = x_AG(5);
    vzro = x_AG(6);
    C_D = x_AG(7);

 % Chamar RK4 com os parâmetros desempacotados
    [xc, yc, zc, vxc,vyc, vzc, vMc, tc] = RK4(xo, yo, zo, vxro, vyro, vzro, C_D); % function[x1, y1, z1, vxr,vyr, vzr, vMr, t] = RK4(xo, yo, zo, vxro, vyro, vzro, C_D)


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
    tamanhoMaximo = max([length(xm), length(ym), length(zm), length(xc), length(yc), length(zc)]);%, length(sigma_x), length(sigma_y), length(sigma_z)]);
    xm = padArrayToLength(xm, tamanhoMaximo);
    ym = padArrayToLength(ym, tamanhoMaximo);
    zm = padArrayToLength(zm, tamanhoMaximo);
    xc = padArrayToLength(xc, tamanhoMaximo);
    yc = padArrayToLength(yc, tamanhoMaximo);
    zc = padArrayToLength(zc, tamanhoMaximo);



 % Calcula a função objetivo
    

    % Inicializa a função objetivo
    f = 0;  
    for k = 1:size(xm, 1)
        % dx = (xm(k) - xc(k)) / sigma_x(k).^2;
        % dy = (ym(k) - yc(k)) / sigma_y(k).^2;
        % dz = (zm(k) - zc(k)) / sigma_z(k).^2;

         % dx = (xm(k) - xc(k)); 
         % dy = (ym(k) - yc(k)); 
         % dz = (zm(k) - zc(k)); 


    sigma_x = wgn(length(xm), 1, 1); 
    sigma_y = wgn(length(xm), 1, 1);
    sigma_z = wgn(length(xm), 1, 1);

    sigma_x = std(sigma_x);
    sigma_y = std(sigma_y);
    sigma_z = std(sigma_z);
   
     % 
      dx  = (xm(k) - xc(k)) / (sqrt(2)*sigma_x);
      dy  = (ym(k) - yc(k)) / (sqrt(2)*sigma_y);
     dz  = (zm(k) - zc(k)) / (sqrt(2)*sigma_z);

        f = f + dx^2 + dy^2 + dz^2;        
    end

    erro = f;
end









