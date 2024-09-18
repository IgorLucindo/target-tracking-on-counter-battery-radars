
close all
clear 
clc

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

%MÉTODO DE RUNGE KUTTA - Validar Modelo
%load v.mat;   % carregando vetor v do 81mm PRODAS 
hr=0.05;    % h = t(i+1)-t(i)  % tamanho do passo
tam = 900; % tamanho amostra
N = tam*hr; % deta t

t = 0:hr:N;                  % tamanho do intervalo de t
% t = t + tr2A(1);
%vx = zeros(1,tam);
%vy = zeros(1,tam);
%vz = zeros(1,tam);
%vM= zeros(1,tam);

theta_0 = 60;%260.6385; 
phi_0 = 70; 
v_0 = 200;

% deltaX = (xE2A(3)-xE2A(1));
% deltaY = (yN2A(3)-yN2A(1));
% deltaZ = (zU2A(3)-zU2A(1));
% ANG1 = deltaZ/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2)));
% theta_0 = asin(ANG1); 
% ANG2 = deltaY/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2))*cos(theta_0));
% phi_0 = -acos(ANG2); 
% v_0 = 200;%veloc(1);
vx(1) = v_0*cosd(theta_0)*sind(phi_0); % condições iniciais
vy(1) = v_0*cosd(theta_0)*cosd(phi_0); %vy(1) = -0.268;
vz(1) = v_0*sind(theta_0);             % condições iniciais
vM(1)= v_0;
%x1 = zeros(1,tam);
%y1 = zeros(1,tam);
%z1 = zeros(1,tam);
x1(1) = 0;     % condições iniciais  -1.4953   -0.2465    1.2056
y1(1) = 0;      % condições iniciais
z1(1) = 0;      % condições iniciais


%===================== Processo de Wiener Método I ========================
% Tw = 0.05; Nw = tam;
% dt = Tw/Nw;
% dW = zeros(1,Nw); % preallocate arrays ...
% W = zeros(1,Nw); % for efficiency
% dW(1) = abs(sqrt(dt)*randn); % first approximation outside the loop ...
% W(1) = dW(1); % since W(0) = 0 is not allowed
% for jw = 2:Nw
% dW(jw) = sqrt(0.5*dt)*randn; % general increment
% W(jw) = W(jw-1) + dW(jw);
% end
load('W_wiener.mat');
% ===== Método II Desvio a partir de um valor nominal a priori conhecido Cd* 
C_d0 = 0.1;%0.0873; 
%delta = wgn(tam,1,0);

C_d = C_d0*exp(W);
%C_d = C_d0*ones(length(C_d),1);
C_D = p*S*C_d/(2*m);
C_Dx = C_D;
C_Dy = C_D;
C_Dz = C_D;
%==================== Equações Diferenciais================================

F_vvx = @(C_Dx,vM,vx) -C_Dx.*vM.*vx';       % ax
F_vvy = @(C_Dy,vM,vy) -C_Dy.*vM.*vy' ;      % ay
F_vvz = @(C_Dz,vM,vz) -C_Dz.*vM.*vz' - g ;  % az           % função EDO 

F_vx = @(vx) vx; % vx
F_vy = @(vy) vy; % vy
F_vz = @(vz) vz; % vz

for i=1:(tam-1)                              % acelerações - velocidades comp X
    k1 = hr*F_vvx(C_Dx(i),vM(i),vx(i));
    k2 = hr*F_vvx(C_Dx(i),vM(i),vx(i)+0.5*k1);
    k3 = hr*F_vvx(C_Dx(i),vM(i),vx(i)+0.25*(k1+k2));
    k4 = hr*F_vvx(C_Dx(i),vM(i),vx(i)-(k2+2*k3));
    vx(i+1) = vx(i) + (1/6)*(k1+4*k3+k4);        

    k1_ = hr*F_vx(vx(i));                      % velocidades - posições comp X
    k2_ = hr*F_vx(vx(i)+0.5*k1_);
    k3_ = hr*F_vx(vx(i)+0.25*(k1_+k2_));
    k4_ = hr*F_vx(vx(i)-(k2_+2*k3_));
    
    x1(i+1) = x1(i) + (1/6)*(k1_+4*k3_+k4_);        
%==========================================================================

    Q1 = hr*F_vvy(C_Dy(i),vM(i),vy(i));             % acelerações - velocidades comp Y
    Q2 = hr*F_vvy(C_Dy(i),vM(i),vy(i)+0.5*Q1);
    Q3 = hr*F_vvy(C_Dy(i),vM(i),vy(i)+0.25*(Q1+Q2));
    Q4 = hr*F_vvy(C_Dy(i),vM(i),vy(i)-(Q2+2*Q3));
    vy(i+1) = (vy(i) + (1/6)*(Q1+4*Q3+Q4));      

    Q1_ = hr*F_vy(vy(i));                            % velocidades - posições comp Y
    Q2_ = hr*F_vy(vy(i)+0.5*Q1_);
    Q3_ = hr*F_vy(vy(i)+0.25*(Q1_+Q2_));
    Q4_ = hr*F_vy(vy(i)-(Q2_+2*Q3_));

    y1(i+1) = y1(i) + (1/6)*(Q1_+4*Q3_+Q4_);     
%==========================================================================

    J1 = hr*F_vvz(C_Dz(i),vM(i),vz(i));                % acelerações - velocidades comp Z
    J2 = hr*F_vvz(C_Dz(i),vM(i),vz(i)+0.5*J1);
    J3 = hr*F_vvz(C_Dz(i),vM(i),vz(i)+0.25*(J1+J2));
    J4 = hr*F_vvz(C_Dz(i),vM(i),vz(i)-(J2+2*J3));
    vz(i+1) = vz(i) + (1/6)*(J1+4*J3+J4);    

    J1_ = hr*F_vz(vz(i));                              % velocidades - posições comp Z
    J2_ = hr*F_vz(vz(i)+0.5*J1_);
    J3_ = hr*F_vz(vz(i)+0.25*(J1_+J2_));
    J4_ = hr*F_vz(vz(i)-(J2_+2*J3_));

    z1(i+1) = z1(i) + (1/6)*(J1_+4*J3_+J4_);
   

vM(i+1) = sqrt((vx(i+1)^2)+(vy(i+1)^2)+(vz(i+1)^2)); 

    if z1(i)<0
        z(i+1) = z1(i);
        iend = i-1;
        break
    end
iend=i;
end
% vm = sqrt((vx.^2)+(vy.^2)+(vz.^2));
x1 = x1(1:iend);    y1 = y1(1:iend);  z1 = z1(1:iend);  t = t(1:iend);
%==========================================================================
xMM2A = x1';
yMM2A = y1';
zMM2A = z1';
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

%==========================================================================

% Filtro de Kalman EKF=====================================================

% xE2A = xMM2A + 1*rand(length(xMM2A),1);
% yN2A = yMM2A + 1*rand(length(yMM2A),1);
% zU2A = zMM2A + 1*rand(length(zMM2A),1);


wx = wgn(length(xMM2A), 1, 1)
wy = wgn(length(xMM2A), 1, 1);
wz = wgn(length(xMM2A), 1, 1);

wx2 = 0*wgn(length(xMM2A), 1, 1);
wy2 = 0*wgn(length(xMM2A), 1, 1);
wz2 = 0*wgn(length(xMM2A), 1, 1);
% wx = zeros(length(xMM2A),1);
% wy = zeros(length(xMM2A),1);
% wz = zeros(length(xMM2A),1);

xE2A = xMM2A + wx;
yN2A = yMM2A + wy;
zU2A = zMM2A + wz;



T = 0.05;                    % Período de Amostragem Radar (T =50ms)
t =tr2A(1:end-2);            %0: T: total_time; % vetor de tempo
n_estados = 7;               % Quantidade de estados / 7 estados e 3 medidas
po =  [xE2A(1); yN2A(1); zU2A(1)]; %  posição inicial do projétil
v_ = v_(:,1);          % velocidade inicial do projétil
%v_ = v_(:,1:end-1);
%v_=[10;10;10];
tam=length(t);

% Inicialização das variáveis==============================================
%alpha = zeros(1,tam);        % parâmetro de arrasto - alpha = S*Cd/m (m^2/kg)
%alpha(1) = 10^-4;%S*C_d0/m;% 5.7256*10^-4;%3.9000e-04;% 10^-6;           % inicialização do parâmetro de arrasto
%lamb = zeros(1,tam);
p0 = 1.203411;                % densidade local
%p0 = 0.002378;
h = 1.0361*10^-4;    %3.158*10^-5 (1/feet) = 1.0361*10^-4 (1/m) ---- ft = m x 3.28084 
Bheta = zeros(1,tam); 
Bheta(1) = 10^4;%0.8*8.7328e+03;

%lamb(1) =alpha(1)*p0;
%lamb(1) =  2.7483e-04;
%lamb(1) =0;
% m = 4.036
% d = 81.03e-3;
% S = pi*(d^2)/4;
%xk_k = zeros(n_estados,tam);   % vetor de estados x
%f_xk = zeros(n_estados,tam);   % vetor de estados x' (derivada do vetor x)

%px = zeros(tam,1);
%py = zeros(tam,1);
%pz = zeros(tam,1);
%Vx = zeros(tam,1);
%Vy = zeros(tam,1);
%Vz = zeros(tam,1);
VBheta = zeros(tam,1);
%Cd_est = zeros(tam,1);
% 


w  = [wx wy wz]';%wgn(3, tam, 1); % Ruído das medidas w
w2  = [wx2 wy2 wz2]'; 
% sig2x = 100*ones(1,tam);
% sig2y = sig2x;
% sig2z = sig2x;
% sigxy = 0*rand(1,tam);
% sigxz = 0*rand(1,tam);
% sigyz = 0*rand(1,tam);



sig2x = (30*wgn(1, tam,1));    %*rand(1,tam);
sig2x = var(sig2x)*ones(1,tam);
sig2y = sig2x;
sig2z = sig2x;
sigxy = 0*rand(1,tam);
sigxz = 0*rand(1,tam);
sigyz = 0*rand(1,tam);



 %load cov_medidas


%R = zeros(3,3); % Matriz de Covariância do erro das medidas (inclui os erros da transformação de coordenadas e ruídos/incertezas dos sensores de medição AER)

% R_o = [sig2x(1) sigxy(1) sigxz(1); 
%        sigxy(1) sig2y(1) sigyz(1);
%        sigxz(1) sigyz(1) sig2z(1)]; 

R_o = eye(3)*sig2x(1);
% =================== EKF =================================================
% Inicialização do vetor de estados, 7 estados: x(0|0) = [x(0) y(0) z(0) vx(0) vy(0) vz(0) CD(0)]
%v_(:,1) = [100;100;100];

xk_k(:,1)  = [po; v_(:,1) ; Bheta(1)];    % Definir estimativa inicial do estado
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
%sig2lamb = 10^-6;
sig2bheta = (10^3)^2;%10^-3;

% Inicialização da Matriz de covariância de estimativa do estado Re_o(0|0) ou P (0|0) - Definir covariância de erro inicial 
P_k_k = [R_o zeros(3,3) zeros(3,1); %(7x7)
      zeros(3,3) sig2v*eye(3,3) zeros(3,1);
      zeros(1,3) zeros(1,3) sig2bheta];

 
%P_k_k = diag([10^6 10^6 10^6 10^6 10^6 10^6 12.86*(10^-14)*exp(-7.38*(10^-5)*xk_k(3,1))]);



% Vetor das medidas reais==================================================
z_med = [xE2A(1:tam) yN2A(1:tam) zU2A(1:tam)]'; % Medidas
H=[eye(3) zeros(3,4)];  % -  Dim:(3x7) Matriz de saída H
% Matriz de Ganho T  -  Dim: (7 x 4)
T7 = [((T^2)/2)*eye(3)    zeros(3,1);
       T*eye(3)         zeros(3,1); 
       zeros(1,3)         1       ];


% Matriz F de transição de estados=========================================
%==========================================================================
vv = 0*wgn(4, tam,1); % Ruído do processo v
%w  = 0*wgn(3, tam, 1); % Ruído das medidas w
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
%g = 9.80665;                  % aceleração da gravidade m/s^2
p0 = 1.203411;                % densidade local
%p0 = 0.002378;
%h = 2.926*10^-5 + (10^-10*xk_k(3,k-1))*3.28084;
p = p0*exp(-h*xk_k(3,k));      % modelo exponencial da densidade do ar
paux(k) = p;
Pp = -(1/2)*p*V;            % variável P para auxilar nos cálculos


R = [sig2x(k) sigxy(k) sigxz(k);
     sigxy(k) sig2y(k) sigyz(k);
     sigxz(k) sigyz(k) sig2z(k)]

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

q_x =10^3;%(5)^2;
q_y =q_x;
q_z =q_x;
%q_lamb =(2.2*10^-3)*exp(-h*xk_k(3,k));% 10^-6
%q_lamb = 7.25*(10^-13)*exp(-0.00005*xk_k(3,k));
q_Bheta =(8*10^3)^2;%(4*10^2)^2; 6*10^2

%q_Bheta =  (7.25*10^-6)*exp(-0.00005/xk_k(3,k));%10^-6; %  1.5600e-04;   1*10^-4

Q_k = [q_x        0          0        0; %10^-8
      0           q_y        0        0;
      0           0          q_z      0;
      0           0          0        q_Bheta ]; %5*10^-8 %(7.25*10^-13)*exp(-0.00005*xk_k(3,k-1))

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
Bheta = xk_k(7,k);
%alpha(k-1) = xk_k(7,k-1);
ee =0;
%===== Elementos da Matriz Jacobiana===================================================================================

f43 = -vx*Pp*h/xk_k(7,k);   % xk_k(7,k) = alpha(k)
%f43 = 0;
f44 = (Pp/xk_k(7,k))*(1+((vx^2))/V^2);
f45 = vx*vy*(1/xk_k(7,k))*Pp/V^2;
f46 = vx*vz*Pp/(xk_k(7,k)*V^2);
f47 = -Pp*vx/xk_k(7,k)^2;

f53 = -vy*Pp*h/xk_k(7,k);   % xk_k(7,k) = alpha(k)
%f53 = 0;
f54 = vy*vx*Pp/(xk_k(7,k)*V^2);
f55 = (Pp/xk_k(7,k))*(1+((vy^2))/V^2);
f56 = vy*vz*Pp/(xk_k(7,k)*V^2);
f57 = -vy*Pp/(xk_k(7,k)^2);

f63 = -vz*Pp*h/xk_k(7,k);   % xk_k(7,k) = alpha(k)
%f63 = 0;
f64 = vx*vz*Pp/(xk_k(7,k)*V^2);
f65 = vy*vz*Pp/(xk_k(7,k)*V^2);
f66 = (Pp/xk_k(7,k))*(1+((vz^2)/V^2));
f67 = -Pp*vz/xk_k(7,k)^2;

f73 = 0;%(h^2)*vz*xk_k(7,k);   % xk_k(7,k) = alpha(k)
%f73 = 0;
f74 = 0;
f75 = 0;
f76 = 0;%-xk_k(7,k)*h;
f77 = 0; %-h*vz;



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
x_4 = (Pp/xk_k(7,k))*x_1;
x_5 = (Pp/xk_k(7,k))*x_2;
x_6 = (Pp/xk_k(7,k))*x_3 - g;
x_7 = 0;


f_xk = [x_1; x_2; x_3; x_4; x_5; x_6; x_7];

%xkm1_k = xk_k(:,k) + f_xk*T +  F*f_xk*(T^2)/2 + T7*vv(:,k);
%xkm1_k = xk_k(:,k) + f_xk*T + F*f_xk*(T^2)/2 ;
%xkm1_k = xk_k(:,k) + F*T ;
  

fi_k = eye(7) + F*T;
P_km1_k = fi_k*(P_k_k)*fi_k'+Q;
xkm1_k = xk_k(:,k) + f_xk*T + F*f_xk*(T^2)/2 ;
%P_km1_k = fi_k*(P_k_k+Q)*fi_k'; 
%xkm1_k = xk_k(:,k) + f_xk*T;% + F*f_xk*(T^2)/2 ;

%P_km1_k = F*P_k_k*F'+Q; 
%P_km1_k = F*(P_k_k+Q)*F';
% Observações
%z_km1 = H*xk_k(:,k) + w(:,k-1); % z_km1(:,k) = H*xkm1_k;
z_km1 = H*xkm1_k+ w2(:,k) ; %z_km1(:,k) = H*xkm1_k;   
% Atualização medidas
K = P_km1_k*H'/(H*P_km1_k*H' + R);     % K = P_km1_k*H'*inv(H*P_km1_k*H' + R);               % Cálculo do ganho de Kalman
norm(K)
%xk_k(:,k) = xk_k(:,k) + K*(z_med(:,k) - z_km1);
xkm1_km1 = xkm1_k + K*(z_med(:,k) - z_km1);    % Atualiza os estados estimados
%P_km1_k = (eye(n_estados)-K*H)*P_km1_k; 
P_km1_km1 = (eye(n_estados)-K*H)*P_km1_k;           % Atualiza a matriz de Covariancia do erro

xk_k(:,k+1) = xkm1_km1;
P_k_k = P_km1_km1;

% px(k)  = xk_k(1);
% py(k)  = xk_k(2);
% pz(k)  = xk_k(3);
% Vx(k)  = xk_k(4);
% Vy(k)  = xk_k(5);
% Vz(k)  = xk_k(6);
% Valpha(k-1)  = xk_k(7);

px(k)  = xkm1_km1(1);
py(k)  = xkm1_km1(2);
pz(k)  = xkm1_km1(3);
Vx(k)  = xkm1_km1(4);
Vy(k)  = xkm1_km1(5);
Vz(k)  = xkm1_km1(6);

VBheta(k)  = xkm1_km1(7);
Cd_est(k) = m/(S*VBheta(k));%paux(k)

% px(k)  = xkm1_k(1);
% py(k)  = xkm1_k(2);
% pz(k)  = xkm1_k(3);
% Vx(k)  = xkm1_k(4);
% Vy(k)  = xkm1_k(5);
% Vz(k)  = xkm1_k(6);
% Valpha(k)  = xkm1_k(7);
%alpha(k) = xk_k(7);
end
erro = abs((C_d(1:length(Cd_est))-Cd_est)./C_d(1:length(Cd_est)))*100;
%Vest = sqrt((Vx.^2)+((Vy.^2)+((Vz.^2));





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



fig4=figure(4); 
%subplot(5,1,[1,3]); % subplot(2,2,[3,4]);
%plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
plot3(xE2A(1:tam),yN2A(1:tam),zU2A(1:tam),'k','Linewidth',1); hold on;
plot3(px(1:end),py(1:end),pz(1:end),'b--','Linewidth',2);
legend('Teórica','Estimada'); grid
%legend('Lan\c cadora','Medidas $P_1$','Estimativa $P_1$','Location','northeast','Interpreter','latex','FontSize',12);grid

xlabel('x (m)','Interpreter','latex','FontSize',12)
ylabel('y (m)','Interpreter','latex','FontSize',12)
zlabel('z (m)','Interpreter','latex','FontSize',12)
title('Trajet\''oria Estimada com EKF','Interpreter','latex','FontSize',13)
%AZIMUTE=-60.6000; ELEVACAO=8.1364; 
%AZIMUTE=-61.5000; ELEVACAO=4.2000; 
 AZIMUTE=41.8635; ELEVACAO=31.2000; 
 view(AZIMUTE, ELEVACAO);
%exportgraphics(fig4, 'fig_EKF_Cd_0.1TW.pdf', 'ContentType', 'vector');


% subplot(5,1,4);
% plot(t(5:end),Valpha(5:end),'k-','Linewidth',1); grid 
% title('Coeficiente de Arrasto - $P_{w(k)}$= 60 dBm','Interpreter','latex','FontSize',13);
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% legend('$\hat{C}_{D}(k)$','Location','northwest','Interpreter','latex','FontSize',14);
% %exportgraphics(fig4, 'fig_EKF3.pdf', 'ContentType', 'vector');

fig5=figure(5);
subplot(3,2,1:2);
plot(t(1:end),Cd_est(1:end),'r-','Linewidth',1); grid ; hold on;
esc_min = C_d-0.05;
esc_max = C_d+0.05;
plot(t(1:end),C_d(1:length(t)),'k-','Linewidth',1); %axis([t(5) t(end) esc_min esc_max])
title('Coeficiente de Arrasto','Interpreter','latex','FontSize',13);
xlabel('tempo(s)','Interpreter','latex','FontSize',12)
legend('$\hat{C}_{D}(k)$','$C_{D_{REAL}}(k)$','Location','Southeast','Interpreter','latex','FontSize',14);


subplot(3,2,3:4);
plot(t(1:end),Cd_est(1:end),'r-','Linewidth',1); grid ; hold on;
plot(t(1:end),C_d(1:length(t)),'k-','Linewidth',1); %axis([t(1) t(end) 0.08 0.115])
title('Coeficiente de Arrasto','Interpreter','latex','FontSize',13);
xlabel('tempo(s)','Interpreter','latex','FontSize',12)
legend('$\hat{C}_{D}(k)$','$C_{D_{REAL}}(k)$','Location','North East','Interpreter','latex','FontSize',14); %Southeast  NorthEast

subplot(3,2,5);
plot(t(1:end),erro,'k-','Linewidth',1); grid ; hold on;

title('Erro da estimativa','Interpreter','latex','FontSize',13);
xlabel('tempo(s)','Interpreter','latex','FontSize',12)
legend('$e(t)$','Location','North East','Interpreter','latex','FontSize',14); %Southeast  NorthEast

subplot(3,2,6);
plot(t(60:end),erro(60:end),'k-','Linewidth',1); grid ; hold on;
stem(t(60),erro(60),'r','Linewidth',1)
title('Erro da estimativa','Interpreter','latex','FontSize',13);
xlabel('tempo(s)','Interpreter','latex','FontSize',12)
legend('$e(t))$','$e(3)$','Location','northwest','Interpreter','latex','FontSize',14); %Southeast  NorthEast
%exportgraphics(fig5, 'fig_EKF_Cd_0.1W.pdf', 'ContentType', 'vector');
rms=rms(C_d(1:length(Cd_est))-Cd_est)
exportgraphics(fig5, 'fig_EKF_BETA2.pdf', 'ContentType', 'vector');


%exportgraphics(fig5, 'fig_EKF_Cd_0.1W.pdf', 'ContentType', 'vector');

% subplot(3,1,3);
% plot(t(100:end),Cd_est(100:end),'r-','Linewidth',1); grid ; hold on;
% plot(t(100:end),C_d(100:length(t)),'k--','Linewidth',1.5); axis([t(1) t(end) 0.09 0.12])
% title('Coeficiente de Arrasto','Interpreter','latex','FontSize',13);
% xlabel('tempo(s)','Interpreter','latex','FontSize',12)
% legend('$\hat{C}_{D}(k)$','$C_{D_{REAL}}(k)$','Location','NorthEast','Interpreter','latex','FontSize',14);
%exportgraphics(fig5, 'fig_EKF_Cd_0.9.pdf', 'ContentType', 'vector');

% % %==========================================================================
% %MÉTODO DE RUNGE KUTTA - Validar Modelo
% %load v.mat;   % carregando vetor v do 81mm PRODAS 
% h1=0.05;    % h = t(i+1)-t(i)                                         % tamanho do passo
% N = tam*h1;
% %tam = length(v);
% %tam = 1000;
% t1 = 0:h1:N;                  % tamanho do intervalo de t
% t = t1 + tr2A(1);
% %v =v(1:end-1);
% vx1 = zeros(1,tam);
% vy1 = zeros(1,tam);
% vz1 = zeros(1,tam);
%  vM= zeros(1,tam);
% % deltaX = (xMM2A(2)-xMM2A(1));
% % deltaY = (yMM2A(2)-yMM2A(1));
% % deltaZ = (zMM2A(2)-zMM2A(1));
% 
% deltaX = (xE2A(2)-xE2A(1));
% deltaY = (yN2A(2)-yN2A(1));
% deltaZ = (zU2A(2)-zU2A(1));
% ANG1 = deltaZ/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2)));
% theta_0 = asin(ANG1); 
% ANG2 = deltaY/(sqrt((deltaX^2)+(deltaY^2)+(deltaZ^2))*cos(theta_0));
% phi_0 = -acos(ANG2); 
% vx1(1) = v(1)*cos(theta_0)*sin(phi_0); % condições iniciais
% vy1(1) = v(1)*cos(theta_0)*cos(phi_0); %vy(1) = -0.268;
% vz1(1) = v(1)*sin(theta_0);             % condições iniciais
% vM(1)= v(1);
% x1 = zeros(1,tam);
% y1 = zeros(1,tam);
% z1 = zeros(1,tam);
% x1(1) = xE2A(1) ;     % condições iniciais  -1.4953   -0.2465    1.2056
% y1(1) = yN2A(1);      % condições iniciais
% z1(1) = zU2A(1);      % condições iniciais
% 
% 
% g = 9.80665;
% %Equações Diferenciais (EDO) - Acelerações para encontrar as velocidades===
% pz(1) = pz(2);
% Valpha(1) = Valpha(2);
% C_Dx = p0*exp(-h*pz).*Valpha/(2);
% C_Dy = C_Dx;
% C_Dz = C_Dy/2;
% 
% 
% 
% F_vvx = @(C_Dx,vM,vx1) -C_Dx.*vM.*vx1';       % ax
% F_vvy = @(C_Dy,vM,vy1) -C_Dy.*vM.*vy1' ;      % ay
% F_vvz = @(C_Dz,vM,vz1) -C_Dz.*vM.*vz1' - g ;    % az           % função EDO 
% 
% F_vx = @(vx1) vx1; % vx
% F_vy = @(vy1) vy1; % vy
% F_vz = @(vz1) vz1; % vz
% 
% 
% for i=1:(tam-1)                              % acelerações - velocidades comp X
%     k1 = h1*F_vvx(C_Dx(i),vM(i),vx1(i));
%     k2 = h1*F_vvx(C_Dx(i),vM(i),vx1(i)+0.5*k1);
%     k3 = h1*F_vvx(C_Dx(i),vM(i),vx1(i)+0.25*(k1+k2));
%     k4 = h1*F_vvx(C_Dx(i),vM(i),vx1(i)-(k2+2*k3));
%     vx1(i+1) = vx1(i) + (1/6)*(k1+4*k3+k4);        
% 
%     k1_ = h1*F_vx(vx1(i));                             % velocidades - posições comp X
%     k2_ = h1*F_vx(vx1(i)+0.5*k1_);
%     k3_ = h1*F_vx(vx1(i)+0.25*(k1_+k2_));
%     k4_ = h1*F_vx(vx1(i)-(k2_+2*k3_));
%     
%     x1(i+1) = x1(i) + (1/6)*(k1_+4*k3_+k4_);  
%       
% %==========================================================================
% 
%     Q1 = h1*F_vvy(C_Dy(i),vM(i),vy1(i));                % acelerações - velocidades comp Y
%     Q2 = h1*F_vvy(C_Dy(i),vM(i),vy1(i)+0.5*Q1);
%     Q3 = h1*F_vvy(C_Dy(i),vM(i),vy1(i)+0.25*(Q1+Q2));
%     Q4 = h1*F_vvy(C_Dy(i),vM(i),vy1(i)-(Q2+2*Q3));
%     vy1(i+1) = (vy1(i) + (1/6)*(Q1+4*Q3+Q4));      
% 
%     Q1_ = h1*F_vy(vy1(i));                            % velocidades - posições comp Y
%     Q2_ = h1*F_vy(vy1(i)+0.5*Q1_);
%     Q3_ = h1*F_vy(vy1(i)+0.25*(Q1_+Q2_));
%     Q4_ = h1*F_vy(vy1(i)-(Q2_+2*Q3_));
% 
%     y1(i+1) = y1(i) + (1/6)*(Q1_+4*Q3_+Q4_);
%      
% %==========================================================================
% 
%     J1 = h1*F_vvz(C_Dz(i),vM(i),vz1(i));                % acelerações - velocidades comp Z
%     J2 = h1*F_vvz(C_Dz(i),vM(i),vz1(i)+0.5*J1);
%     J3 = h1*F_vvz(C_Dz(i),vM(i),vz1(i)+0.25*(J1+J2));
%     J4 = h1*F_vvz(C_Dz(i),vM(i),vz1(i)-(J2+2*J3));
%     vz1(i+1) = vz1(i) + (1/6)*(J1+4*J3+J4);    
% 
%     J1_ = h1*F_vz(vz1(i));                             % velocidades - posições comp Z
%     J2_ = h1*F_vz(vz1(i)+0.5*J1_);
%     J3_ = h1*F_vz(vz1(i)+0.25*(J1_+J2_));
%     J4_ = h1*F_vz(vz1(i)-(J2_+2*J3_));
% 
%     z1(i+1) = z1(i) + (1/6)*(J1_+4*J3_+J4_);
% 
% vM(i+1) = sqrt((vx1(i+1)^2)+(vy1(i+1)^2)+(vz1(i+1)^2)); 
% 
%     if z1(i)<0
%         %z(i+1) = z1(i);
%         iend = i-1;
%         break
%     end
% % iend=i;
% end
% 
% vm = sqrt((vx1.^2)+(vy1.^2)+(vz1.^2));
% 
% figure
% plot3(x1,y1,z1,'g--','linewidth',2); title('Validação do Modelo'); hold on; 
% plot3(xlanc,ylanc,zlanc,'gx','linewidth',3); hold on; 
% plot3(xE2A(1:tam),yN2A(1:tam),zU2A(1:tam),'k','Linewidth',1); hold on;
% %plot3(xE2A2(1:tam),yN2A2(1:tam),zU2A2(1:tam),'k--','Linewidth',1); 
% %hold on;
% plot3(px(2:end),py(2:end),pz(2:end),'b--','Linewidth',2);
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

%  w_k = [(tau/2) * Delta_v_o_k; Delta_v_o_k];
%     x_k_km1 = A * x_km1_km1 - w_k;
%     Pe_k_km1 = A * Pe_km1_km1 * A';
%     aux_den = x_k_km1(1)^2 + x_k_km1(2)^2;
%     c_k = [x_k_km1(2)/aux_den; -x_k_km1(1)/aux_den; 0; 0];
%     g_k = Pe_k_km1 * c_k / (c_k' * Pe_k_km1 * c_k + sigma2_k);
%     x_k_k = x_k_km1 + g_k * (beta_k - arctan_BOT(x_k_km1(1:2,:),1,1));
%     Pe_k_k = (eye(4) - g_k * c_k') * Pe_k_km1;


