clear
clc

passo = 0.0001;
tf = 10;

HR = 60;
tc = 60/HR;

t = 0:passo:tf;
N = tf/passo;

% PARAMETERS

% Sistemic arterial section
Rcs = 0.0495;
Ls = 0.0000736;
Cas = 3.47;
Ras = 0.7905;
% Ras = 1.6;

% Pulmonary arterial section
Rcp = 0.00225;
Lp = 0.000035;
Cap = 4;
Rap = 0.075;

% Sistemic venous section
Cvs = 84;
Rvs = 0.06;

% Pulmonary venous section
Cvp = 5;
Rvp = 0.00075;

% Left heart
Rm = 0.003;
Ra = 0.0075;

% Right heart
Rt = 0.00525;
Rp = 0.00225;

% Constants
Vlv0 = 4.0;
Vla0 = 4.0;
% Vrv0 = 10.0;
Vrv0 = 20.0;
Vra0 = 4.0;
Plv0 = 0.0;
Pla0 = 0.0;
Prv0 = 0.0;
Pra0 = 0.0;

% VAD parameters
Ri = 0.0677;
Ro = 0.0677;
Li = 0.0127;
Lo = 0.0127;
Rk = 0;
alpha = -3.5;
B0 = 0.17070;
B1 = 0.02177;
B2 = -0.000093;
x1 = 1.0;
w = 6000*(2*pi)/60;

% DIODES
Dm = 0;
Da = 0;
Dt = 0;
Dp = 0;

Dm_v = [Dm];
Da_v = [Da];
Dt_v = [Dt];
Dp_v = [Dp];

% ELASTANCES 
[Elv, Erv, Ea] = elastance(passo,tf,tc);

% VARIABLES

Vla = zeros(1,N);
Vlv = zeros(1,N);
Qvad = zeros(1,N);
Qa = zeros(1,N);
Pas = zeros(1,N);
Pvs = zeros(1,N);
Vra = zeros(1,N);
Vrv = zeros(1,N);
Pap = zeros(1,N);
Qp = zeros(1,N);
Pvp = zeros(1,N);

Pla = zeros(1,N);
Plv = zeros(1,N);
Pra = zeros(1,N);
Prv = zeros(1,N);

Pao = zeros(1,N);

LEDV = [];
LESV = [];
REDV = [];
RESV = [];

% INITIAL CONDITIONS

Vla(1) = 48;
Vlv(1) = 140; 
Qvad(1) = 0;
Qa(1) = 0; 
Pas(1) = 75; 
Pvs(1) = 9.25; 
Vra(1) = 17; 
Vrv(1) = 120; 
Pap(1) = 9; 
Qp(1) = 0; 
Pvp(1) = 7.4;

Pla(1) = Ea(1)*(Vla(1) - Vla0);
Plv(1) = Elv(1)*(Vlv(1) - Vlv0);
Pra(1) = Ea(1)*(Vra(1) - Vra0);
Prv(1) = Erv(1)*(Vrv(1) - Vrv0);

Pao(1) = Rcs*Qa(1) + Pas(1);

x = [Vla(1), Vlv(1), Qvad(1), Qa(1), Pas(1), Pvs(1), Vra(1), Vrv(1), Pap(1), Qp(1), Pvp(1)]';

for i=1:N
    
    tn = rem(t(i),tc);
    
    if Pla(i) >= Plv(i)
        Dm = 1;
    else
        Dm = 0;
    end
    
    if Plv(i) >= Pas(i)
        Da = 1;
    else
        Da = 0;
    end
    
    if Pra(i) >= Prv(i)
        Dt = 1;
    else
        Dt = 0;
    end
    
    if Prv(i) >= Pap(i)
        Dp = 1;
    else
        Dp = 0;
    end
    
%     if rem(t(i),tc) == 0
%         LEDV = [LEDV Vlv(i)];
%         REDV = [REDV Vrv(i)];
%     elseif rem(t(i)/round(0.4*tc,2)) == 0
%         LESV = [LESV Vlv(i)];
%         RESV = [REDV Vrv(i)];
%     end

    if Plv(i) <= x1
%         Rk = alpha*(Plv(i) - x1);
        Rk = 0;
    else
        Rk = 0;
    end
    
    Dm_v = [Dm_v Dm];
    Da_v = [Da_v Da];
    Dt_v = [Dt_v Dt];
    Dp_v = [Dp_v Dp];
    
    % EQUATIONS
    
    RR = Ri + Ro + Rk + B0;
    LL = Li + Lo + B1;
    
    a11 = -((1/Rvp)+(Dm/Rm))*Ea(i);
    a12 = (Dm/Rm)*Elv(i);
    a111 = 1/Rvp;
    
    a21 = (Dm/Rm)*Ea(i);
    a22 = -((Dm/Rm)+(Da/Ra))*Elv(i);
    
    a32 = Elv(i);
    a33 = -RR;
    
    a44 = -Rcs;
    
    a55 = -1/Ras;
    a56 = 1/Ras;
    
    a65 = 1/Ras;
    a66 = -((1/Ras)+(1/Rvs));
    a67 = Ea(i)/Rvs;
   
    a76 = 1/Rvs;
    a77 = -((1/Rvs)+(Dt/Rt))*Ea(i);
    a78 = (Dt/Rt)*Erv(i);
    
    a87 = (Dt/Rt)*Ea(i);
    a88 = -(Dt/Rt)*Erv(i);
    
    a99 = -1/Rap;
    a911 = 1/Rap;
    
    a108 = Erv(i);
    a1010 = -(Rp + Rcp);
    
    a11_1 = Ea(i)/Rvp;
    a119 = 1/Rap;
    a1111 = -((1/Rap)+(1/Rvp));
    
    % A =    [Vla   Vlv  Qvad   Qa  Pas  Pvs  Vra   Vrv    Pap  Qp   Pvp]
    A = zeros(11,11);
    A(1,:) = [a11  , a12,   0,   0,   0,   0,   0,    0,    0,     0,  a111];         % Vla  
    A(2,:) = [a21  , a22,  -1,   0,   0,   0,   0,    0,    0,     0,     0];         % Vlv
    A(3,:) = [ 0   , a32, a33,   0,   0,   0,   0,    0,    0,     0,     0]/LL;      % Qvad
    A(4,:) = [ 0   ,   0,   0, a44,  -1,   0,   0,    0,    0,     0,     0]*(Da/Ls); % Qa
    A(5,:) = [ 0   ,  0 ,   0,   1, a55, a56,   0,    0,    0,     0,     0]/Cas;     % Pas
    A(6,:) = [ 0   ,  0 ,   0,   0, a65, a66, a67,    0,    0,     0,     0]/Cvs;     % Pvs
    A(7,:) = [ 0   ,  0 ,   0,   0,   0, a76, a77,  a78,    0,     0,     0];         % Vra
    A(8,:) = [ 0   ,  0 ,   0,   0,   0,   0, a87,  a88,    0,    -1,     0];         % Vrv
    A(9,:) = [ 0   ,  0 ,   0,   0,   0,   0,   0,    0,  a99,     1,  a911]/Cap;     % Pap
    A(10,:) =[ 0   ,  0 ,   0,   0,   0,   0,   0, a108,   -1, a1010,     0]*(Dp/Lp); % Qp
    A(11,:)= [a11_1,  0 ,   0,   0,   0,   0,   0,    0, a119,     0, a1111]/Cvp;     % Pvp
 
%     B = zeros(11,1);
    B = [       (((1/Rvp)+(Dm/Rm))*(Ea(i)*Vla0)) - ((Dm/Rm)*(Elv(i)*Vlv0));
                   -(Dm/Rm)*(Ea(i)*Vla0) + ((Dm/Rm)+(Da/Ra))*(Elv(i)*Vlv0);
                             -(Pao(i)/LL) - (B2*(w^2))/LL - (Elv(i)*Vlv0)/LL;
                                                                  Pao(i)/Ls;
                                                                          0;
                                               (1/(Rvs*Cvs))*(- Ea(i)*Vra0);
                 (((1/Rvs)+(Dt/Rt))*(Ea(i)*Vra0)) - ((Dt/Rt)*(Erv(i)*Vrv0));
                                       (Dt/Rt)*(- Ea(i)*Vra0 + Erv(i)*Vrv0);
                                                                          0;
                                                     (1/Lp)*(- Erv(i)*Vrv0);
                                               (1/(Rvp*Cvp))*(- Ea(i)*Vla0)];
    
    if Da == 0
        x(4) = 0;
        B(4) = 0;
    end
    if Dp == 0
        x(10) = 0;
        B(10) = 0;
    end

    x = runkut4(passo, A, x, B);
    
    Vla(i+1) = x(1);
    Vlv(i+1) = x(2);
    Qvad(i+1) = x(3);
    Qa(i+1) = x(4);
    Pas(i+1) = x(5);
    Pvs(i+1) = x(6);
    Vra(i+1) = x(7);
    Vrv(i+1) = x(8);
    Pap(i+1) = x(9);
    Qp(i+1) = x(10);    
    Pvp(i+1) = x(11);
    
    Pla(i+1) = Ea(i+1)*(Vla(i+1) - Vla0);
    Plv(i+1) = Elv(i+1)*(Vlv(i+1) - Vlv0);
    Pra(i+1) = Ea(i+1)*(Vra(i+1) - Vra0);
    Prv(i+1) = Erv(i+1)*(Vrv(i+1) - Vrv0);
    
    Qa_dot = [0 diff(Qa)];
    Pao(i+1) = Rcs*Qa(i+1) + Ls*Qa_dot(end) + Pas(i+1);
       
end

% REDV(end)
% RESV(end)
% LEDV(end) 
% LESV(end)

% SV_right = (REDV(end) - RESV(end));
% SV_left = (LEDV(end) - LESV(end));
% 
% CO_right = HR*0.001*SV_right;
% CO_left = HR*0.001*SV_left;
% 
% ejectionFraction_right = 100*(SV_right/REDV(end));
% ejectionFraction_left = 100*(SV_left/LEDV(end));

figure()
subplot(9,2,1)
plot(t,Vla)
title('Vla')
subplot(9,2,2)
plot(t,Vlv)
title('Vlv')
subplot(9,2,3)
plot(t,Qa)
title('Qa')
subplot(9,2,4)
plot(t,Pas)
title('Pas')
subplot(9,2,5)
plot(t,Pvs)
title('Pvs')
subplot(9,2,6)
plot(t,Vra)
title('Vra')
subplot(9,2,7)
plot(t,Vrv)
title('Vrv')
subplot(9,2,8)
plot(t,Pap)
title('Pap')
subplot(9,2,9)
plot(t,Qp)
title('Qp')
subplot(9,2,10)
plot(t,Pvp)
title('Pvp')

subplot(9,2,11)
plot(t,Pla)
title('Pla')
subplot(9,2,12)
plot(t,Plv)
title('Plv')
subplot(9,2,13)
plot(t,Pra)
title('Pra')
subplot(9,2,14)
plot(t,Prv)
title('Prv')

subplot(9,2,15)
plot(t,Pla,t,Plv,t,Pao)
title('Left heart')
subplot(9,2,16)
plot(t,Pra,t,Prv,t,Pap)
title('Right heart')

subplot(9,2,17)
plot(t,Qvad)
title('Qvad')

