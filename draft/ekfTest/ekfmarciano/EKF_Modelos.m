close all
clear all
clc
%==========================================================================
% Simulações de Monte Carlo
ENS=1;
auxENS = zeros(ENS,1);
delta_Tracking = 10; % 35.95  //4seg. de tracking\\TOTAL 23.85 Atenção: o tempo total de rastreio é 23,8s //10%= 2.3850   //20%= 4.7700 //25%= 5.9625// 30%=7.1550  // 35%=8.3475 //40%=9.5400 //45%= 10.7325 // 50%= 11.9250
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
erro=zeros(1000,1);
erro_rms = zeros(ENS,1);
erro_rse = 0;
tempo_processamento = 0;
for ens=1:ENS
tic
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
tam = 1500;%850; % tamanho amostra ESTAVA 1000!
N = tam*hr; % deta t

t = 0:hr:N;                  % tamanho do intervalo de t
% t = t + tr2A(1);
vxr = zeros(1,tam);
vyr = zeros(1,tam);
vzr = zeros(1,tam);
vMr= zeros(1,tam);

theta_0 = 70;%260.6385; 
phi_0 = 60; 
v_0 = 200;
range_0 = 800;
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
v_0 = 200;%veloc(1);
% vx1(1) = v_0*cosd(phi_0)*cosd(theta_0); % condições iniciais
% vy1(1) = v_0*cosd(phi_0)*sind(theta_0); %vy(1) = -0.268;
% vz1(1) = v_0*sind(phi_0);             % condições iniciais


vMr(1)= v_0;
x1 = zeros(1,tam);
y1 = zeros(1,tam);
z1 = zeros(1,tam);
% x1(1) = 0;     % condições iniciais  -1.4953   -0.2465    1.2056
% y1(1) = 0;      % condições iniciais
% z1(1) = 0;      % condições iniciais

x1(1) = range_0*cosd(phi_0)*cosd(theta_0); % condições iniciais
y1(1) = range_0*cosd(phi_0)*sind(theta_0); %vy(1) = -0.268;
z1(1) = range_0*sind(phi_0);             % condições iniciais
%[x1(1) y1(1) z1(1)]

%===================== Processo de Wiener: Método para adicionar aleatóridade a trajetória ========================
Tw = 0.05; Nw = tam;
dt = Tw/Nw;
dW = zeros(1,Nw); % preallocate arrays ...
W = zeros(1,Nw); % for efficiency
dW(1) = (sqrt(dt)*randn); % first approximation outside the loop ...
W(1) = dW(1); % since W(0) = 0 is not allowed
for jw = 2:Nw
dW(jw) = sqrt(0.5*dt)*randn; % general increment
W(jw) = W(jw-1) + dW(jw);
end
%load('W_wiener.mat');
% ===== Método II Desvio a partir de um valor nominal a priori conhecido Cd* 
C_d0 = 0.1; 
%delta = wgn(tam,1,0);

%
C_d = C_d0*exp(W(1:tam));
%C_d = C_d0*ones(1,tam);
C_D = p*S*C_d/(2*m);
C_Dx = C_D;
C_Dy = C_D;
C_Dz = C_D;
%[x1,y1,z1,t] = RK4(x1(1), y1(1), z1(1), vxr(1), vyr(1), vzr(1), C_D);
%==================== Equações Diferenciais================================
% 
F_vvx = @(C_Dx,vMr,vxr) -C_Dx.*vMr.*vxr;       % ax
F_vvy = @(C_Dy,vMr,vyr) -C_Dy.*vMr.*vyr ;      % ay
F_vvz = @(C_Dz,vMr,vzr) -C_Dz.*vMr.*vzr - g ;  % az           % função EDO 

F_vx = @(vxr) vxr; % vx
F_vy = @(vyr) vyr; % vy
F_vz = @(vzr) vzr; % vz

for i=1:(tam-1)                                 % acelerações - velocidades comp X
    k1 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i));
    k2 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i)+0.5*k1);
    k3 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i)+0.25*(k1+k2));
    k4 = hr*F_vvx(C_Dx(i),vMr(i),vxr(i)-(k2+2*k3));
    vxr(i+1) = vxr(i) + (1/6)*(k1+4*k3+k4);        

    k1_ = hr*F_vx(vxr(i));                      % velocidades - posições comp X
    k2_ = hr*F_vx(vxr(i)+0.5*k1_);
    k3_ = hr*F_vx(vxr(i)+0.25*(k1_+k2_));
    k4_ = hr*F_vx(vxr(i)-(k2_+2*k3_));

    x1(i+1) = x1(i) + (1/6)*(k1_+4*k3_+k4_);        
%==========================================================================

    Q1 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i));             % acelerações - velocidades comp Y
    Q2 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i)+0.5*Q1);
    Q3 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i)+0.25*(Q1+Q2));
    Q4 = hr*F_vvy(C_Dy(i),vMr(i),vyr(i)-(Q2+2*Q3));
    vyr(i+1) = (vyr(i) + (1/6)*(Q1+4*Q3+Q4));      

    Q1_ = hr*F_vy(vyr(i));                            % velocidades - posições comp Y
    Q2_ = hr*F_vy(vyr(i)+0.5*Q1_);
    Q3_ = hr*F_vy(vyr(i)+0.25*(Q1_+Q2_));
    Q4_ = hr*F_vy(vyr(i)-(Q2_+2*Q3_));

    y1(i+1) = y1(i) + (1/6)*(Q1_+4*Q3_+Q4_);     
%==========================================================================

    J1 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i));               % acelerações - velocidades comp Z
    J2 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i)+0.5*J1);
    J3 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i)+0.25*(J1+J2));
    J4 = hr*F_vvz(C_Dz(i),vMr(i),vzr(i)-(J2+2*J3));
    vzr(i+1) = vzr(i) + (1/6)*(J1+4*J3+J4);    

    J1_ = hr*F_vz(vzr(i));                              % velocidades - posições comp Z
    J2_ = hr*F_vz(vzr(i)+0.5*J1_);
    J3_ = hr*F_vz(vzr(i)+0.25*(J1_+J2_));
    J4_ = hr*F_vz(vzr(i)-(J2_+2*J3_));

    z1(i+1) = z1(i) + (1/6)*(J1_+4*J3_+J4_);


vMr(i+1) = sqrt((vxr(i+1)^2)+(vyr(i+1)^2)+(vzr(i+1)^2)); 

    if z1(i)<0
        z1(i) = z1(i-1);
        y1(i) = y1(i-1);
        x1(i) = x1(i-1);
        iend = i-1;
        break
    end

end
% vm = sqrt((vx.^2)+(vy.^2)+(vz.^2));
x1 = x1(1:iend);  
y1 = y1(1:iend); 
z_aux = z1; 
z1 = z1(1:iend); 
t_aux = t; 

t = t(1:iend);
%==========================================================================
% xMM2A = x1(80:end)';
% yMM2A = y1(80:end)';
% zMM2A = z1(80:end)';
% t=t(80:end);

 xMM2A = x1(1:end)';
 yMM2A = y1(1:end)';
 zMM2A = z1(1:end)';
 t=t(1:end);


tr2A = t';
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

wx = wgn(length(xMM2A), 1, 1);
wy = wgn(length(xMM2A), 1, 1);
wz = wgn(length(xMM2A), 1, 1);

wx_ = wgn(length(xMM2A), 1, 1);
wy_ = wgn(length(xMM2A), 1, 1);
wz_ = wgn(length(xMM2A), 1, 1);

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
vv =  0*(10^-6)*wgn(4, tam, 1);  % Ruído do processo v 0.00001
w  = [wx wy wz]';%wgn(3, tam, 1); % Ruído das medidas w
% 
% sig2x = 1000*ones(1,tam);
% sig2y = sig2x;
% sig2z = sig2x;
% sigxy = 1*rand(1,tam);
% sigxz = 1*rand(1,tam);
% sigyz = 1*rand(1,tam);
 %load cov_medidas
sig2x = (30*wgn(1, tam,1));    %*rand(1,tam);
sig2x = var(sig2x)*ones(1,tam);
sig2y = sig2x;
sig2z = sig2x;
sigxy = 0*rand(1,tam);
sigxz = 0*rand(1,tam);
sigyz = 0*rand(1,tam);

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

sig2v= 150^2;       %
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

q_x =10^-3;   %0.9*10^3;
q_y =q_x;
q_z =q_x;
%q_lamb =(2.2*10^-3)*exp(-h*xk_k(3,k));% 10^-6
%q_lamb = 7.25*(10^-11)*exp(-0.00005*xk_k(3,k));

q_lamb = 10^-12;%  10^-3

%q_lamb =  (10^-12)*exp(-xk_k(3,k));%10^-6; %  1.5600e-04;   1*10^-4

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
f66 = -xk_k(7,k)*((V^2)+vz^2)/(2*V+ee);%f66 = -xk_k(7,k)*(vz^2)/(2*V+ee);
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
tk(k) = t(1) + k*T;
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
   %  if xkm1_km1(3)<0
   %  k_ext=k-1;
   % % break;
   %  end

% % %==========================================================================
% %MÉTODO DE RUNGE KUTTA - PARA EXTRAPOLAR TRAJETÓRIA
% h1e=0.05;    % h = t(i+1)-t(i)                                         % tamanho do passo
% tame=tam;
% Ne = tame*h1e;
% tam = length(v);
% tam = 1000;
% t1e = 0:h1e:Ne-1;                  % tamanho do intervalo de t
% t = t1e + tr2A(1);
% v =v(1:end-1);
% vx1e = zeros(1,tame);
% vy1e = zeros(1,tame);
% vz1e = zeros(1,tame);
% vMe= zeros(1,tame);
% 
% % vMe(1)= sqrt((Vx(1)^2)+(Vy(1)^2)+(Vz(1)^2));
% % deltaX = px(1);
% % deltaY = py(1);
% % deltaZ = pz(1);
% 
% % vMe(1)= sqrt((Vx(1)^2)+(Vy(1)^2)+(Vz(1)^2));
% % deltaX = (px(3)-px(2));
% % deltaY = (py(3)-py(2));
% % deltaZ = (pz(3)-pz(2));
% % 
% vMe(1)= sqrt((vxr(1)^2)+(vyr(1)^2)+(vzr(1)^2));
% deltaX = (x1(2)-x1(1));
% deltaY = (y1(2)-y1(1));
% deltaZ = (z1(2)-z1(1));
% % 
% % vMe(1)= sqrt((Vx(1)^2)+(Vy(1)^2)+(Vz(1)^2));
% % deltaX = (xE2A(2)-xE2A(1));
% % deltaY = (yN2A(2)-yN2A(1));
% % deltaZ = (zU2A(2)-zU2A(1));
% 
% ANG1e = deltaZ/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2)));
% theta_0e = asin(ANG1e); 
% ANG2 = deltaY/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2))*cos(theta_0e));
% phi_0 = acos(ANG2); 
% vx1e(1) = vMe(1)*cos(theta_0e)*sin(phi_0); % condições iniciais
% vy1e(1) = vMe(1)*cos(theta_0e)*cos(phi_0); %vy(1) = -0.268;
% vz1e(1) = vMe(1)*sin(theta_0e);             % condições iniciais
% 
% % vx1e(1) = Vx(2); % condições iniciais da velocidade
% % vy1e(1) = Vy(2); %vy(1) = -0.268;
% % vz1e(1) = Vz(2);             % condições iniciais
% 
% x1e = zeros(1,tame);
% y1e = zeros(1,tame);
% z1e = zeros(1,tame);
% % x1e(1) = px(1) ;     % condições iniciais  -1.4953   -0.2465    1.2056
% % y1e(1) = py(1);      % condições iniciais
% % z1e(1) = pz(1);      % condições iniciais
% 
% x1e(1) = x1(1) ;     % condições iniciais  -1.4953   -0.2465    1.2056
% y1e(1) = y1(1);      % condições iniciais
% z1e(1) = z1(1);      % condições iniciais
% 
% % g = 9.80665;
% % Equações Diferenciais (EDO) - Acelerações para encontrar as velocidades===
% % pz(1) = pz(2);% ??????
% % Valpha(1) = Valpha(2);
% % 
% % 
% % ===================== Processo de Wiener: Método para adicionar aleatóridade a trajetória EXTRAPOLAÇÃO========================
% % Twe = 0.05; Nwe = 1000;
% % dte = Twe/Nwe;
% % dWe = zeros(1,Nwe); % preallocate arrays ...
% % We = zeros(1,Nwe); % for efficiency
% % 
% % dWe(1) = (sqrt(dte)*randn); % first approximation outside the loop ...
% % We(1) = dWe(1); % since W(0) = 0 is not allowed
% % for jwe = 2:Nwe
% %     alpha_dwe = var(Valpha(1:40))/dte;
% %     alpha_dwe=1;
% % dWe(jwe) = sqrt(alpha_dwe*dte)*randn; % general increment
% % We(jwe) = We(jwe-1) + dWe(jwe);
% %end
% %Cd_est(30)=C_d(1);
% %Cd_e = Cd_est(40)*exp(We);
% %C_Dxe = p*S*Cd_e/(2*m);
% C_Dxe = p*S*Cd_est(40)*ones(1,tam)/(2*m);
% C_Dye = C_Dxe;
% C_Dze = C_Dye;
% 
% 
% 
% 
% 
% F_vvx = @(C_Dxe,vMe,vx1e) -C_Dxe.*vMe.*vx1e';       % ax
% F_vvy = @(C_Dye,vMe,vy1e) -C_Dye.*vMe.*vy1e' ;      % ay
% F_vvz = @(C_Dze,vMe,vz1e) -C_Dze.*vMe.*vz1e' - g ;    % az           % função EDO 
% 
% F_vx = @(vx1e) vx1e; % vx
% F_vy = @(vy1e) vy1e; % vy
% F_vz = @(vz1e) vz1e; % vz
% 
% 
% for i=1:(tam)                              % acelerações - velocidades comp X
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
%         % z1e(i) = z1e(i-1);
%         % y1e(i) = y1e(i-1);
%         % x1e(i) = x1e(i-1);
%         % 
%         x1e = x1e(1:iende);    y1e = y1e(1:iende);  z1e = z1e(1:iende);  vMe= vMe(1:iende); 
% break
% 
%     end
% 
% end

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

erro = abs((C_d(1:length(Cd_est))'-Cd_est')./C_d(1:length(Cd_est))')*100 ; %+erro(1:length(Cd_est))
%erro_rse=sqrt((C_d(1:length(Cd_est))' - Cd_est').^2) + erro_rse; %sqrt(mean((observed - predicted).^2));
%erro_aux = erro+erro_aux
%erroX = abs((xMM2A(1:length(x1))-x1')./xMM2A(1:length(x1)))*100 + erroX;
tempo_processamento = toc 
end
%Fim do Loop Monte Carlo
%erro = erro/ENS;
erro_rms=sum(erro_rms)/ENS;
 %tempo_processamento = tempo_processamento/ENS

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

% figPIP_EKF_completo = figure('Position', [100, 100, 900, 700]);
% plot(xaux2,yaux2, 'r.','linewidth',5); hold on; grid
% plot(cord_ponto_central(1,1),cord_ponto_central(1,2),'b+','linewidth',7)
% xlabel('x (m)','Interpreter','latex','FontSize',15)
% ylabel('y (m)','Interpreter','latex','FontSize',15)
% legend('Pontos estimados EKF','Centroide (algoritmo \textit{K-Means})','Interpreter','latex','FontSize',16);
% exportgraphics(figPIP_EKF_completo, 'figPIP_EKF_completo_nom.pdf', 'ContentType', 'vector');






% fig2=figure(2); % Gráfico dos Estados
% 
% % subplot(511);
% % plot(t,xE2A(1:tam),'k','Linewidth',1); grid
% % hold on
% % plot(t(2:end),px(2:end),'k--','Linewidth',1);
% % hold on
% % plot(t,yN2A(1:tam),'r','Linewidth',1);
% % hold on
% % plot(t(2:end),py(2:end),'r--','Linewidth',1);
% % hold on
% % plot(t,zU2A(1:tam),'b','Linewidth',1);
% % hold on
% % plot(t(2:end),pz(2:end),'b--','Linewidth',1);
% % title('Posições');
% % xlabel('tempo(s)');
% % legend('x_{Medido}','x_{Real}');
% % 
% % subplot(512);
% % plot(t,Vx,'k','Linewidth',1); grid
% % hold on;
% % plot(t,Vy,'r','Linewidth',1);
% % hold on;
% % plot(t,Vz,'b','Linewidth',1);
% % title('Velocidades Estimadas');
% % xlabel('tempo(s)');
% % legend('v_x','v_y','v_z' );
% 
% subplot(311);
% %plot(t,alpha,'k','Linewidth',1);
% %hold on;
% plot(t(5:end),Valpha(5:end),'g','Linewidth',1); grid
% title('Parâmetro de Arrasto');
% xlabel('tempo(s)');
% legend('$\hat{\alpha}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %exportgraphics(fig2, 'fig_kalman0.pdf', 'ContentType', 'vector');
% 
% subplot(312);
% %plot(t,alpha,'k','Linewidth',1);
% %hold on;
% plot(t(5:end),1./(Valpha(5:end)+10^-5),'b','Linewidth',1); grid
% title('Parâmetro Balístico');
% xlabel('tempo(s)');
% legend('$\hat{\beta}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %exportgraphics(fig2, 'fig_kalman0.pdf', 'ContentType', 'vector');
% 
% subplot(313);
% 
% plot(t(5:end),Valpha(5:end).*(m*S/p),'k-','Linewidth',1); grid
% title('Coeficiente de Arrasto');
% xlabel('tempo(s)');
% legend('$\hat{C}_{D}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
%legend('\hat{C}_{D}');


% fig3=figure(3);
% 
% %plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% plot3(xE2A(1:tam),yN2A(1:tam),zU2A(1:tam),'k','Linewidth',1); hold on;
% %plot3(xE2A2(1:tam),yN2A2(1:tam),zU2A2(1:tam),'k--','Linewidth',1); 
% %hold on;
% %plot3(px(2:end),py(2:end),pz(2:end),'b--','Linewidth',2);
% % grid
% % xlabel('x (m)')
% % ylabel('y (m)')
% % zlabel('z (m)')
% % AZIMUTE=-74.8675; ELEVACAO=15.9026;    
% %  view(AZIMUTE, ELEVACAO);
% title('Trajetória 3D Nominal');
% legend('Nominal','Saída EKF'); grid
% %legend('Lan\c cadora','Medidas $P_1$','Traj. Estimada $P_1$','Location','northeast','Interpreter','latex','FontSize',12);grid
% %legend('Modelo $P_1$','Medidas','','Lan\c cadora','Ponto de impacto','Interpreter','latex','FontSize',12);grid 
% xlabel('x (m)','Interpreter','latex','FontSize',12)
% ylabel('y (m)','Interpreter','latex','FontSize',12)
% zlabel('z (m)','Interpreter','latex','FontSize',12)
% title('Trajet\''orias 3D ','Interpreter','latex','FontSize',13)
% AZIMUTE=41.8635; ELEVACAO=31.2000; 
% view(AZIMUTE, ELEVACAO);
% exportgraphics(fig3, 'fig_kalmanCd_0.9.pdf', 'ContentType', 'vector');



% fig4=figure(4); 
% %subplot(5,1,[1,3]); % subplot(2,2,[3,4]);
% %plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k','Linewidth',1); hold on;
% plot3(px(1:end),py(1:end),pz(1:end),'b--','Linewidth',2);
% legend('Teórica','Estimada'); grid
% %legend('Lan\c cadora','Medidas $P_1$','Estimativa $P_1$','Location','northeast','Interpreter','latex','FontSize',12);grid
% 
% xlabel('x (m)','Interpreter','latex','FontSize',12)
% ylabel('y (m)','Interpreter','latex','FontSize',12)
% zlabel('z (m)','Interpreter','latex','FontSize',12)
% title('Trajet\''oria Estimada com EKF','Interpreter','latex','FontSize',13)
% %AZIMUTE=-60.6000; ELEVACAO=8.1364; 
% %AZIMUTE=-61.5000; ELEVACAO=4.2000; 
%  AZIMUTE=41.8635; ELEVACAO=31.2000; 
%  view(AZIMUTE, ELEVACAO);
% %exportgraphics(fig4, 'fig_EKF_Cd_0.1TW.pdf', 'ContentType', 'vector');


% subplot(5,1,4);
% plot(t(5:end),Valpha(5:end),'k-','Linewidth',1); grid 
% title('Coeficiente de Arrasto - $P_{w(k)}$= 60 dBm','Interpreter','latex','FontSize',13);
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% legend('$\hat{C}_{D}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %exportgraphics(fig4, 'fig_EKF3.pdf', 'ContentType', 'vector');

fig5=figure('Position', [100, 100, 800, 600]);
subplot(3,2,1:2);
plot(tk(1:end),Cd_est(1:end),'r-','Linewidth',1); grid ; hold on;
esc_min = C_d-0.05;
esc_max = C_d+0.05;
plot(t(1:length(t)),C_d(1:length(t)),'k-.','Linewidth',1); %axis([t(5) t(end) esc_min esc_max])
%title('Estimativa do Coeficiente de Arrasto','Interpreter','latex','FontSize',16);
xlabel('tempo(s)','Interpreter','latex','FontSize',16)
legend('$\mathrm{\hat{C}_{d_{ESTIMADO}}}$','$\mathrm{C_{d_{REAL}}}$','Location','northeast','Interpreter','latex','FontSize',16);


subplot(3,2,3:4);
plot(tk(1:end),Cd_est(1:end),'r-','Linewidth',1); grid ; hold on;
title('Estimativa a partir de 2 seg.','Interpreter','latex','FontSize',14);
plot(t(1:end),C_d(1:length(t)),'k-.','Linewidth',1); axis([t(40) t(length(t)) 0 0.2])
%title('Coeficiente de Arrasto','Interpreter','latex','FontSize',13);
xlabel('tempo(s)','Interpreter','latex','FontSize',16)
legend('$\mathrm{\hat{C}_{d_{ESTIMADO}}}$','$\mathrm{C_{d_{REAL}}}$','Location','northeast','Interpreter','latex','FontSize',16); %Southeast  NorthEast

subplot(3,2,5);
plot(tk(1:end),erro(1:end),'k-','Linewidth',1); grid ; hold on;
title('Erro total da estimativa ','Interpreter','latex','FontSize',16);
%title('Erro da estimativa','Interpreter','latex','FontSize',13);
xlabel('tempo(s)','Interpreter','latex','FontSize',16); ylabel('erro (\%)','Interpreter','latex','FontSize',16)
legend('$\mathrm{e(t)}$','Location','North East','Interpreter','latex','FontSize',16); %Southeast  NorthEast
%exportgraphics(fig5, 'fig_EKF_Cd_0.1W.pdf', 'ContentType', 'vector');

subplot(3,2,6);
plot(tk(40:end),erro(40:end),'k-','Linewidth',1); grid ; hold on;
stem(tk(40),erro(40),'r','Linewidth',1)
title('Erro da estimativa a partir de 2 seg.','Interpreter','latex','FontSize',16);
xlabel('tempo(s)','Interpreter','latex','FontSize',16); ylabel('erro (\%)','Interpreter','latex','FontSize',16)
legend('$\mathrm{e(t))}$','$\mathrm{e(2s)}$','Location','northwest','Interpreter','latex','FontSize',16); %Southeast  NorthEast
exportgraphics(fig5, 'fig_erro_LAMB_nominal10Kalman.pdf', 'ContentType', 'vector');
%rms=rms(C_d(1:length(Cd_est))-Cd_est)
% subplot(3,1,3);
% plot(t(100:end),Cd_est(100:end),'r-','Linewidth',1); grid ; hold on;
% plot(t(100:end),C_d(100:length(t)),'k--','Linewidth',1.5); axis([t(1) t(end) 0.09 0.12])
% title('Coeficiente de Arrasto','Interpreter','latex','FontSize',13);
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% legend('$\hat{C}_{D}(k)$','$C_{D_{REAL}}(k)$','Location','NorthEast','Interpreter','latex','FontSize',14);
%exportgraphics(fig5, 'fig_EKF_Cd_0.9.pdf', 'ContentType', 'vector');








%======================== Plotando TLE  ===================================
% xaux1 = randn(100, 1);
% yaux1 = 2*xaux1 + 0.5*randn(100, 1);
% 
x_ipV = x1(end);
y_ipV = y1(end);
x_ip = cord_ponto_central(1,1); 
y_ip = cord_ponto_central(1,2);
% x_pd = rmmissing(xaux1);
% y_pd = rmmissing(yaux1);

x_pd = (xaux1);
y_pd =(yaux1);

% x_pd = nonzeros(xaux1);
% y_pd = nonzeros(yaux1);

% x_ip = mean(x_pd);
% y_ip = mean(y_pd);

%mu_x ==x_ip;
%mu_y ==y_ip;
var_x_pd = var(x_pd);
cov_x_pd_y_pd = sum((x_pd - mean(x_pd)) .* (y_pd - mean( y_pd))) / (length(x_pd) - 1);
var_y_pd =  var(y_pd);


vetor_ip = [x_ip, y_ip];  % Vetor de médias %meanVector = [mu_x,mu_y];


erro_x =  (x_ip - x_pd);
erro_y =  (y_ip - y_pd);

%covarianceMatrix = cov(x_pd,y_pd); 
covarianceMatrix = cov(erro_x,erro_y); 



%covarianceMatrix = [var_x_pd, cov_x_pd_y_pd; cov_x_pd_y_pd, var_y_pd];  % Matriz de covariância
%covarianceMatrix = cov(x_pd,y_pd); 
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

%======================== Cáculo do TLE  ==================================
r_TLE_s = sqrt((x_pd-x_ipV).^2+(y_pd-y_ipV).^2);
r_TLE_s_ord = sort(r_TLE_s);
TLE_porcent = 99;
%TLE_amostra_porcent = (TLE_porcent*ens)/100;
TLE_amostra_porcent = (TLE_porcent*length(x_pd))/100;TLE_amostra_porcent=round(TLE_amostra_porcent);
r_TLE_final = r_TLE_s_ord(TLE_amostra_porcent );

theta_TLE = linspace(0,2*pi);
x_TLE = r_TLE_final*cos(theta_TLE) + x_ipV(end);
y_TLE = r_TLE_final*sin(theta_TLE) + y_ipV(end);

%============
% Calculo dos vetores de erro e medidas de dispersão

% x_pd = xaux1;
% y_pd = yaux1;
% vet_erroX = x_ipV - x_pd;
% vet_erroY = y_ipV - y_pd;
vet_erroX = x_ipV - x_pd;
vet_erroY = y_ipV - y_pd;
Mod_Erro = sqrt(vet_erroX.^2 + vet_erroY.^2); % módulo dos vetores de erro
% Calcula o desvio padrão dos erros
desv_erros = std(Mod_Erro)
% Calcula o IQR dos erros
% Q1_erro = prctile(Mod_Erro, 25);  % Primeiro quartil (25%)
% Q3_erro = prctile(Mod_Erro, 75);  % Terceiro quartil (75%)
% iqr_erro = Q3 - Q1;

% Calcula o MSE dos erros
mse = mean(Mod_Erro.^2)
desv_mse_erros = std(mse)
rmse=sqrt(mse);

% %Plota os pontos da elipse
% figure;
% plot(points(:,1),points(:,2), 'g-.'); axis equal; hold on;
% plot(ellipsePoints(:,1),ellipsePoints(:,2), 'k-.'); axis equal;
% plot(x_pd,y_pd, 'ro'); hold on;
% ellipse(vetor_ip, covarianceMatrix, confidenceLevel);
% 
% 
% % 
% figScatter=figure(7);
% % Plota o centro da elipse
% scatter(vetor_ip(1), vetor_ip(2), 'r', 'filled');
% 
% title('Elipse de Confiança para Distribuição Normal Multivariada');
% xlabel('X');
% ylabel('Y');
% grid on;
% %axis equal;

% figElipse=figure(8);
% % Use a função ellipse do File Exchange para desenhar a elipse
% ellipse(vetor_ip, covarianceMatrix, confidenceLevel); hold on; 
% ellipse(vetor_ip, covarianceMatrixC, confidenceLevelC); hold on; 
% plot(x_ip,y_ip,'rx','Linewidth',2);


% fig6=figure('Position', [100, 100, 800, 600]); 
% %subplot(5,1,[1,3]); % subplot(2,2,[3,4]);
% %plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% %plot3(xE2A(1:tam),yN2A(1:tam),zU2A(1:tam),'k','Linewidth',1); hold on;
% %stem3(xE2A(1),yN2A(1),zU2A(1),'b','Linewidth',1); hold on;
% plot3(xE2A(1),yN2A(1),zU2A(1),'m+','Linewidth',10); hold on;
% plot3(xE2A(1:end),yN2A(1:end),zU2A(1:end),'k','Linewidth',1); hold on;
% plot3(xE2A(end),yN2A(end),0,'b+','Linewidth',5); hold on;
% %plot3(px(1:end),py(1:end),pz(1:end),'b--','Linewidth',2); hold on;
% plot3(x1(1:20),y1(1:20),z1(1:20),'go','Linewidth',1); hold on;
% plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r','Linewidth',1); hold on;
% plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
% 
% legend('Ponto de detec\c{c}\~ao','Trajet\''oria real','Ponto de impacto real','Medidas radar ($\Delta$t = 2s)','Trajet\''oria extrapolada','Ponto de impacto estimado','Location','North East','Interpreter','latex','FontSize',16); grid
% %legend('Lan\c cadora','Medidas $P_1$','Estimativa $P_1$','Location','northeast','Interpreter','latex','FontSize',12);grid
% 
% xlabel('x (m)','Interpreter','latex','FontSize',16)
% ylabel('y (m)','Interpreter','latex','FontSize',16)
% zlabel('z (m)','Interpreter','latex','FontSize',16)
% %title('Trajet\''oria Estimada com $\lambda$ EKF','Interpreter','latex','FontSize',13); 
% %AZIMUTE=-60.6000; ELEVACAO=8.1364; 
% %AZIMUTE=-61.5000; ELEVACAO=4.2000; 
% %AZIMUTE=41.8635; ELEVACAO=31.2000; 
% AZIMUTE=105.9135; ELEVACAO=21.5202;
%  view(AZIMUTE, ELEVACAO);
% exportgraphics(fig6, 'fig_EKF_traj_nominal1W.pdf', 'ContentType', 'vector');


% figPIP=figure('Position', [100, 100, 700, 500]);
% %plot_ellipse(xaux1,yaux1,x1(end),y1(end))    ; hold on
% plot(x_pd,y_pd,'ro','linewidth',1); daspect([1 1 1]);
% %legend('PIP Previsto com $\Delta t = 40$ \quad seg','Interpreter','latex','FontSize',14);
% %grid
% hold on; 
% plot(x_ipV,y_ipV,'b+','linewidth',10); hold on; 
% 
% %error_ellipse(covarianceMatrix,vetor_ip); hold on
% %plot(px(end),py(end),'y+','linewidth',3);daspect([1 1 1]); hold on;
% 
% ellipse(vetor_ip, covarianceMatrixC, confidenceLevelC);hold on; 
% ellipse(vetor_ip, covarianceMatrix, confidenceLevel); %daspect([1 1 1]); hold on;
% 
% 
% 
% %plot_ellipse(x_pd,y_pd, x_pd,y_pd); %plot_ellipse(x_pd,y_pd, x_ip,y_ip); %plot_ellipse( x,y,X0,Y0 )  
% 
% %plot(x_TLE,y_TLE, 'v-.')
% fill(x_TLE,y_TLE, 'g', 'FaceAlpha', 0.1);
% %legend('PIP Estimado com $\Delta t = 15,55$ \quad seg ');
% 
% %legend('PIP Estimado com $\Delta t = 4$seg','PIP Estimado com $\Delta t = 15,55$seg (todo o rastreio)','Interpreter','latex','FontSize',12);
% %title('Pontos de Impacto Estimados','Interpreter','latex','FontSize',15)
% xlabel('x (m)','Interpreter','latex','FontSize',16)
% ylabel('y (m)','Interpreter','latex','FontSize',16);
% Legend_TLE = num2str(r_TLE_final,'TLE 90 (r = %.2fm)');
% legend('PI$_\mathrm{E}$ ($\Delta t = 2$ seg.)','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',16); %%'95\% ellipse (intervalo de confian\c{c}a)',
% 
% exportgraphics(figPIP, 'figPIP_nominal20W_.pdf', 'ContentType', 'vector');

%[x2,y2,z2,t2] = RK4(x1e(end)-x1e(1), y1e(end)-y1e(1), z1e(1),-Vx(end),-Vy(end), -Vz(end), C_D,tam);   %Valpha/2

%[x2,y2,z2,t2] = RK4(x1(end)-x1(1), y1(end)-y1(1), z1(1),-vxr(1),-vyr(1), vzr(1), C_D,tam);
tamJanela = 60; 

bm = (1/tamJanela)*ones(1,tamJanela);
x_filtMM = filtfilt(bm,1,xE2A);
y_filtMM = filtfilt(bm,1,yN2A);
z_filtMM = filtfilt(bm,1,zU2A);


%taux = zeros(length(x1e));
tz0 = -z_filtMM(1)/(z_filtMM(2)-z_filtMM(1)); 
xLP = x_filtMM(1) +tz0*(x_filtMM(2)-x_filtMM(1));
yLP = y_filtMM(1) +tz0*(y_filtMM(2)-y_filtMM(1));

%===========================================================================
% 
[z_max, index] = max(z_aux); % Encontra o valor máximo em z e o índice correspondente
 t_max = t_aux(index); % Usa o índice para encontrar o tempo correspondente
% vz_00 = g*t_max;
% vx_00 = x1(end)/(2*t_max); % testar com t(end) 2*t_max
% vy_00 = y1(end)/(2*t_max);
% x00 = x1(end) - vx_00*2*t_max;
% y00 = y1(end) - vy_00*2*t_max;
% %z00 = z_max - 0.5*g*(t(end)/2)^2;
% z00 = 0;

 % Considerando uma trajetória parabólicas para aproximar o primeiro
 % "chute" do ponto de lançamento
 [z_max, index] = max(z1e); % identificar o máximo da trajetória
  x_max = x1e(index);       % identificar as componente x e y correspondentes a esse máximo
  y_max = y1e(index);

x00 = 2*x_max - x1e(end);
y00 = 2*y_max - y1e(end);
z00 = 0;
%===========================================================================
fig_PIP_TRAJ = figure('Position', [150, 150, 750, 550]);
subplot(2,1,1)
plot3(xE2A(1),yN2A(1),zU2A(1),'kx','Linewidth',2); hold on;
%plot3([xLP x1e(5)],[yLP y1e(5)],[0 z1e(5)],'r','Linewidth',1); hold on;
%plot3([x00 x1e(5)],[y00 y1e(5)],[ z00 z1e(5)],'b','Linewidth',1); hold on;
plot3(x1,y1,z1,'k','Linewidth',1); hold on;
tr2A = tr2A - tr2A(1);
scatter3(xE2A(1:amostra_parada), yN2A(1:amostra_parada), zU2A(1:amostra_parada), 3, tr2A(1:amostra_parada), 'filled'); hold on;
colormap(jet);  % Utiliza o mapa de cores 'jet' para representar o tempo
% Supondo que você queira representar faixas específicas do mapa de cores
hold on;
scatter3(xE2A(amostra_parada),yN2A(amostra_parada),zU2A(amostra_parada), 30, 'g', 'filled');hold on;

%plot3(x2,y2,z2,'r','linewidth',2); hold on;
%colorbar;       % Adiciona uma barra de cores para indicar a relação com o tempo
cb = colorbar;%('northoutside');  % Coloca a barra de cores na parte superior e na horizontal;  % Adiciona uma barra de cores
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Tempo (s)';  % Legenda da barra de cores
% Ajustando a espessura da barra de cores
cb.Position = cb.Position .* [1.1 1 0.25 1];  % Reduz a espessura pela metade
cb.Label.FontSize = 14;
%plot3(xE2A(end),yN2A(end),0,'b+','Linewidth',5); hold on;
%plot3(x1(1:20),y1(1:20),z1(1:20),'go','Linewidth',1); hold on;
plot3(x1e(1:end),y1e(1:end),z1e(1:end),'r-.','Linewidth',1); hold on;
plot3(x1e(end),y1e(end),z1e(end),'rO','Linewidth',2); hold on;
legend('Ponto de detec\c{c}\~ao','Medida do radar','Fim da detec\c{c}\~ao','Trajet\''oria extrapolada','Ponto de impacto estimado','Location','North East','Interpreter','latex','FontSize',12); grid
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14)
zlabel('z (m)','Interpreter','latex','FontSize',14)
title('Trajet\''oria 100$\%$ do rastreio ','Interpreter','latex','FontSize',15); 
AZIMUTE=105.9135; ELEVACAO=21.5202;
 view(AZIMUTE, ELEVACAO);

subplot(2,1,2)
plot(x_pd,y_pd,'ro','linewidth',1);daspect([1 1 1]);hold on; 
plot(x_ipV,y_ipV,'b+','linewidth',4);daspect([1 1 1]); hold on; 
ellipse(vetor_ip, covarianceMatrixC, confidenceLevelC,'r'); daspect([1 1 1]);hold on; 
ellipse(vetor_ip, covarianceMatrix, confidenceLevel,'k'); 
fill(x_TLE,y_TLE, 'g', 'FaceAlpha', 0.1);
%title('Pontos de Impacto Estimados','Interpreter','latex','FontSize',15)
xlabel('x (m)','Interpreter','latex','FontSize',14)
ylabel('y (m)','Interpreter','latex','FontSize',14);
Legend_TLE = num2str(r_TLE_final,'TLE 99 (r = %.2fm)');
legend('PI$_\mathrm{E}$ ','PI$_\mathrm{V}$','C\''{i}rculo de Alerta 99\% $ \chi^2_2$','Elipse 99\% $\chi^2_2$ ', Legend_TLE, 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',12); %%'95\% ellipse (intervalo de confian\c{c}a)',
%axis([2600 3400 990 1310])
%exportgraphics(fig_PIP_TRAJ, 'figPIP_TRAJ_nominal_1.0.pdf', 'ContentType', 'vector');
%set(gcf,'outerposition',get(0,'screensize'))
%export_fig ('test_full_screen_export_macos','-png','-r300', '-transparent',gcf);

exportgraphics(fig_PIP_TRAJ, 'figPIP_TRAJ_nominal_TOTAL2.pdf', 'ContentType', 'vector');

%export_fig (fig_PIP_TRAJ,'figPIP_TRAJ_nominal_1.0');
%export_fig('figPIP_TRAJ_nominal_2.0.pdf', 'ContentType', 'vector'); %-painters   -opengl
% fig_Ruidos_Medidas = figure('Position', [100, 100, 700, 500]);
% wc1 = wgn(length(t), 1, 1, 'real'); 
% wc2 = wgn(length(t), 1, 10, 'real');
% wc3 = wgn(length(t), 1, 20, 'real');
% subplot(3,1,1)
% plot(t,wc1(1:length(t)), 'k' ); hold on;grid
% xlabel('tempo (s)','Interpreter','latex','FontSize',16)
% ylabel('Intensidade (m)','Interpreter','latex','FontSize',16)
% legend('Cen\''{a}rio I', 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); axis([t(1) t(end) -30 40])%%'95\% ellipse (intervalo de confian\c{c}a)', 
% subplot(3,1,2)
% plot(t,wc2(1:length(t)), 'k');  hold on;grid
% xlabel('tempo (s)','Interpreter','latex','FontSize',16)
% ylabel('Intensidade (m)','Interpreter','latex','FontSize',16)
% legend('Cen\''{a}rio II','Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); axis([t(5) t(end) -30 40])%%'95\% ellipse (intervalo de confian\c{c}a)',
% subplot(3,1,3)
% plot(t,wc3(1:length(t)), 'k'); grid
% xlabel('tempo (s)','Interpreter','latex','FontSize',16)
% ylabel('Intensidade (m)','Interpreter','latex','FontSize',16)
% legend('Cen\''{a}rio III', 'Interpreter', 'latex','Location','northwest','Interpreter','latex','FontSize',15); axis([t(5) t(end) -30 40])%%'95\% ellipse (intervalo de confian\c{c}a)',
% exportgraphics(fig_Ruidos_Medidas, 'fig_Ruidos_Medidas.pdf', 'ContentType', 'vector');


% figModelos=figure('Position', [450, 450, 1000, 800]);
% subplot(4,1,1)
% load('Cd_est_BETA.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'r-.','Linewidth',1); grid ; hold on;
% plot(t(1:end),C_d(1:1:length(t)),'k-','Linewidth',1); 
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% ylabel('$\mathrm{\hat{C}_d}$','Interpreter','latex','FontSize',14)
% set(get(gca, 'ylabel'), 'Rotation', 0);
% title('Modelo $\beta$-EKF','Interpreter','latex','FontSize',15);
% ylim([0.05 0.15]);
% xlim([0 33]);
% legend('$\mathrm{\hat{C}_d}$','$\mathrm{{C}_d}$','Location','North East','Interpreter','latex','FontSize',12); %Southeast  NorthEast
% 
% subplot(4,1,2)
% load('Cd_est_BETA_DELTA.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'r-.','Linewidth',1); grid ; hold on;
% plot(t(1:end),C_d(1:1:length(t)),'k-','Linewidth',1); 
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% ylabel('$\mathrm{\hat{C}_d}$','Interpreter','latex','FontSize',14)
% set(get(gca, 'ylabel'), 'Rotation', 0);
% title('Modelo $\beta_0$-EKF','Interpreter','latex','FontSize',15);
% ylim([0.05 0.15]);
% xlim([0 33]);
% legend('$\mathrm{\hat{C}_d}$','$\mathrm{{C}_d}$','Location','North East','Interpreter','latex','FontSize',12); %Southeast  NorthEast
% 
% subplot(4,1,3)
% load('Cd_est_LN.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'r-.','Linewidth',1); grid ; hold on;
% plot(t(1:end),C_d(1:1:length(t)),'k-','Linewidth',1); 
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% ylabel('$\mathrm{\hat{C}_d}$','Interpreter','latex','FontSize',14)
% set(get(gca, 'ylabel'), 'Rotation', 0);
% title('Modelo $\Delta$-EKF','Interpreter','latex','FontSize',15);
% ylim([0.05 0.15]);
% xlim([0 33]);
% legend('$\mathrm{\hat{C}_d}$','$\mathrm{{C}_d}$','Location','North East','Interpreter','latex','FontSize',12); %Southeast  NorthEast
% 
% subplot(4,1,4)
% load('Cd_est_GAMA.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'r-.','Linewidth',1); grid ; hold on;
% plot(t(1:end),C_d(1:1:length(t)),'k-','Linewidth',1); 
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% ylabel('$\mathrm{\hat{C}_d}$','Interpreter','latex','FontSize',14)
% set(get(gca, 'ylabel'), 'Rotation', 0);
% title('Modelo $\gamma$-EKF','Interpreter','latex','FontSize',15);
% ylim([0.05 0.15]);
% xlim([0 33]);
% legend('$\mathrm{\hat{C}_d}$','$\mathrm{{C}_d}$','Location','North East','Interpreter','latex','FontSize',12); %Southeast  NorthEast
% 


% figModelos2=figure('Position', [400, 400, 900, 700]);
% subplot(2,1,1)
% load('Cd_est_BETA.mat')
% plot(t(1:end),C_d(1:length(t)),'k-','Linewidth',2); hold on;
% plot(t(1:length(Cd_est)),Cd_est(1:end),'r-.','Linewidth',1); grid ; hold on;
% load('Cd_est_BETA_DELTA.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'b-.','Linewidth',1); hold on;
% load('Cd_est_LN.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'g-','Linewidth',1.5); hold on;
% load('Cd_est_GAMA.mat')
% plot(t(1:length(Cd_est)),Cd_est(1:end),'m-','Linewidth',1.5); 
% 
% xlabel('tempo(s)','Interpreter','latex','FontSize',18)
% ylabel('Coeficiente de arrasto','Interpreter','latex','FontSize',18)
% %set(get(gca, 'ylabel'), 'Rotation', 0);
% title('Modelos para Estima\c{c}\~ao do Coeficiente de Arrasto','Interpreter','latex','FontSize',20);
% ylim([0.05 0.18]);
% xlim([0 33]);
% legend('$\mathrm{{C}_d\quad (real)}$','$\mathrm{\hat{C}_d}\quad(\beta$-EKF)','$\mathrm{\hat{C}_d}\quad(\beta_0$-EKF)','$\mathrm{\hat{C}_d}\quad(e$-EKF)','$\mathrm{\hat{C}_d}\quad(\gamma$-EKF)' ,'Location','best','Interpreter','latex','FontSize',16); %Southeast  NorthEast
% 
% 
% 
% subplot(2,1,2)
% load('erro_BETA.mat')
% load('Cd_est_BETA.mat')
% rms_BETA=rms((C_d(1:length(Cd_est)) - Cd_est))
% plot(t(1:length(erro)),erro(1:end),'r-.','Linewidth',1); grid ; hold on;
% load('erro_BETA_DELTA.mat')
% load('Cd_est_BETA_DELTA.mat')
% rms_BETA_DELTA=rms(C_d(1:length(Cd_est))-Cd_est)
% plot(t(1:length(erro)),erro(1:end),'b-.','Linewidth',1); hold on;
% load('erro_LN.mat')
% load('Cd_est_LN.mat')
% rms_LN=rms(C_d(1:length(Cd_est))-Cd_est)
% plot(t(1:length(erro)),erro(1:end),'g-','Linewidth',1.5); hold on;
% load('erro_GAMA.mat')
% load('Cd_est_GAMA.mat')
% rms_GAMA=rms(C_d(1:length(Cd_est))-Cd_est)
% plot(t(1:length(erro)),erro(1:end),'m-','Linewidth',1.5); 
% 
% xlabel('tempo(s)','Interpreter','latex','FontSize',18)
% ylabel('erro$\mathrm{(\%)}$','Interpreter','latex','FontSize',18)
% %set(get(gca, 'ylabel'), 'Rotation', 0);
% title('Erro de Estima\c{c}\~ao ','Interpreter','latex','FontSize',20);
% %ylim([0.05 0.18]);
% xlim([0 33]);
% %legend('$\mathrm{erro(\%)}\quad(\beta$-EKF)\quad RMSE:',num2str(rms_BETA),'$\mathrm{erro(\%)}\quad(\beta_0$-EKF)\quad RMSE:',num2str(rms_BETA_DELTA),'$\mathrm{erro(\%)}\quad(\Delta$-EKF)\quad RMSE:',num2str(rms_LN),'$\mathrm{erro(\%)}\quad(\gamma$-EKF)\quad RMSE:',num2str(rms_GAMA) ,'Location','North East','Interpreter','latex','FontSize',16); %Southeast  NorthEast
% 
% legenda1 = ['$\mathrm{erro(\%)}\quad(\beta$-EKF)\quad RMSE: ', num2str(rms_BETA, '%.4f')];
% legenda2 = ['$\mathrm{erro(\%)}\quad(\beta_0$-EKF)\quad RMSE: ', num2str(rms_BETA_DELTA, '%.4f')];
% legenda3 = ['$\mathrm{erro(\%)}\quad(e$-EKF)\quad RMSE: ', num2str(rms_LN, '%.4f')];
% legenda4 = ['$\mathrm{erro(\%)}\quad(\gamma$-EKF)\quad RMSE: ', num2str(rms_GAMA, '%.4f')];
% 
% % Aplicando as legendas
% legend(legenda1, legenda2, legenda3, legenda4, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
% exportgraphics(figModelos2, 'figModelos.pdf', 'ContentType', 'vector');


figDensidade0=figure(1);%figure('Position', [400, 400, 900, 700]);
subplot(2,1,1)
% plot(t(1:end), pz(1:length(t)), 'k','LineWidth', 2);
% xlabel('tempo (s)', 'Interpreter', 'latex', 'FontSize',13);
% ylabel('Altitude $a$ (m)', 'Interpreter', 'latex', 'FontSize',13);
% %title('Varia\c{c}\~{a}o da Densidade do Ar com a Altitude', 'Interpreter', 'latex','FontSize',16);
% grid on;
% Encontrando o ponto de máximo
[maxVal, maxIdx] = max(pz); % maxVal é o valor máximo, maxIdx é o índice correspondente no vetor
maxTime = t(maxIdx); % Tempo correspondente ao valor máximo

% Arredondando o valor máximo para não ter casas decimais
maxValRounded = round(maxVal);
maxTimeRounded = round(maxTime);
plot(maxTime, maxVal, 'ro','LineWidth', 2); hold on
plot(t(1:length(pz)), pz(1:end), 'k','LineWidth', 2);
xlabel('tempo (s)', 'Interpreter', 'latex', 'FontSize',15);
ylabel('Altitude $z$(m)', 'Interpreter', 'latex', 'FontSize',15);
grid on;


% Adicionando uma anotação ou legenda com o valor máximo arredondado
legend(['$z_{MAX}:$ ', num2str(maxValRounded), ' m em ', num2str(round(maxTime)), ' s'], 'Interpreter', 'latex','Location', 'best', 'FontSize',13);
% Note que também arredondei o tempo para não ter casas decimais


% Destacando o ponto de máximo
hold on;
 % 'ro' cria um marcador vermelho no ponto de máximo
hold off;


subplot(2,1,2)
rho0 = 1.225; % Densidade do ar ao nível do mar em kg/m³
H = 8500; % Escala de altura atmosférica em metros
z = pz; % Faixa de altitude de 0 a 40.000 metros
rho = rho0 * exp(-z/H); % Cálculo da densidade do ar
[minVal_rho, minId_rho] = min(rho); % maxVal é o valor máximo, maxIdx é o índice correspondente no vetor
minZ = z(minId_rho); % Tempo correspondente ao valor máximo
plot(minZ, minVal_rho, 'ko','LineWidth', 2); hold on 
% Arredondando o valor máximo para não ter casas decimais
minVal_rhoRounded = (minVal_rho);
minZRounded = round(minZ);
plot(z, rho, 'r','LineWidth', 2);
xlabel('Altitude $z$(m)', 'Interpreter', 'latex', 'FontSize',15);
ylabel('Densidade do ar $\rho(z)$ (kg/m$^3$)', 'Interpreter', 'latex', 'FontSize',15);
%legend(['%\rho:%', num2str(minVal_rhoRounded), ' Kg/m$^3$ em ', num2str(minZRounded), ' m'], 'Interpreter', 'latex','Location', 'best', 'FontSize',13);
legend(['$\rho= ', num2str(minVal_rhoRounded), '\,kg/m^3$ em ', num2str(minZRounded), 'm'], 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 13);
%title('Varia\c{c}\~{a}o da Densidade do Ar com a Altitude', 'Interpreter', 'latex','FontSize',16);
grid on;
% hold on
% % Definindo a área que queremos ampliar
% x_zoom_min = 958.58415;
% x_zoom_max = 958.5845;
% y_zoom_min = 1.09435620;
% y_zoom_max = 1.09435625;
% % Criando um conjunto de eixos para o recorte
% axes('Position',[.45 .17 .15 .2]) % Posição do recorte [left, bottom, width, height] em proporção da figura
% box on % Cria uma borda ao redor do recorte
% % Plotando o recorte
% plot(z, rho, 'r','LineWidth', 2);
% title('Zoom', 'Interpreter', 'latex','FontSize',13);
% xlim([x_zoom_min x_zoom_max]);
% ylim([y_zoom_min y_zoom_max]);
% grid on;
exportgraphics(figDensidade0, 'figDensidade1.pdf', 'ContentType', 'vector');
% 
% 
% 


% % 
% figDensidade = figure%figure('Position', [400, 400, 900, 700]);
% rho0 = 1.225; % Densidade do ar ao nível do mar em kg/m³
% H = 8500; % Escala de altura atmosférica em metros
% 
% z = 0:100:30000; % Faixa de altitude de 0 a 40.000 metros
% rho = rho0 * exp(-z/H); % Cálculo da densidade do ar
% 
% plot(z, rho, 'r','LineWidth', 2);
% xlabel('Altitude $z$ (m)', 'Interpreter', 'latex', 'FontSize',15);
% ylabel('Densidade do ar $\rho(z)$ (kg/m$^3$)', 'Interpreter', 'latex', 'FontSize',15);
% %title('Varia\c{c}\~{a}o da Densidade do Ar com a Altitude', 'Interpreter', 'latex','FontSize',16);
% grid on;
% exportgraphics(figDensidade, 'figDensidade.pdf', 'ContentType', 'vector');






