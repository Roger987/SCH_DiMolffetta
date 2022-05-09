clear
clc

passo = 0.0001;
tf = 8;

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
% Rm = 0.05;
Ra = 0.0075;

% Right heart
Rt = 0.00525;
Rp = 0.00225;

% Constants
Vlv0 = 5.0;
Vla0 = 4.0;
% Vrv0 = 10.0;
Vrv0 = 20.0;
Vra0 = 4.0;

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

% Vla(1) = 48;
% Vlv(1) = 140; 
% Qa(1) = 0; 
% Pas(1) = 75; 
% Pvs(1) = 9.25; 
% Vra(1) = 17; 
% Vrv(1) = 120; 
% Pap(1) = 9; 
% Qp(1) = 0; 
% Pvp(1) = 7.4;

Vla(1) = 68;
Vlv(1) = 125; 
Qa(1) = 0; 
Pas(1) = 72; 
Pvs(1) = 8.75; 
Vra(1) = 20; 
Vrv(1) = 125; 
Pap(1) = 12; 
Qp(1) = 0; 
Pvp(1) = 10;

Pla(1) = Ea(1)*(Vla(1) - Vla0);
Plv(1) = Elv(1)*(Vlv(1) - Vlv0);
Pra(1) = Ea(1)*(Vra(1) - Vra0);
Prv(1) = Erv(1)*(Vrv(1) - Vrv0);

Pao(1) = Rcs*Qa(1) + Pas(1);

x = [Vla(1), Vlv(1), Qa(1), Pas(1), Pvs(1), Vra(1), Vrv(1), Pap(1), Qp(1), Pvp(1)]';

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
    
    if rem(t(i),tc) == 0
        LEDV = [LEDV Vlv(i)];
        REDV = [REDV Vrv(i)];
%     elseif rem(t(i),tc) == 0.4
    elseif (rem(t(i),tc) >= 0.4 - 0.0000000001) && (rem(t(i),tc) <= 0.4 + 0.0000000001)
        LESV = [LESV Vlv(i)];
        RESV = [RESV Vrv(i)];
    end
    
    Dm_v = [Dm_v Dm];
    Da_v = [Da_v Da];
    Dt_v = [Dt_v Dt];
    Dp_v = [Dp_v Dp];
    
    % EQUATIONS

    
    a11 = -((1/Rvp)+(Dm/Rm))*Ea(i);
    a12 = (Dm/Rm)*Elv(i);
    a110 = 1/Rvp;
    
    a21 = (Dm/Rm)*Ea(i);
    a22 = -(Dm/Rm)*Elv(i);
    
    a32 = Elv(i);
    a33 = -(Ra + Rcs);
    
    a44 = -1/Ras;
    a45 = 1/Ras;
    
    a54 = 1/Ras;
    a55 = -((1/Ras)+(1/Rvs));
    a56 = Ea(i)/Rvs;
    
    a65 = 1/Rvs;
    a66 = -((1/Rvs)+(Dt/Rt))*Ea(i);
    a67 = (Dt/Rt)*Erv(i);
    
    a76 = (Dt/Rt)*Ea(i);
    a77 = -(Dt/Rt)*Erv(i);
    
    a88 = -1/Rap;
    a810 = 1/Rap;
    
    a97 = Erv(i);
    a99 = -(Rp + Rcp);
    
    a101 = Ea(i)/Rvp;
    a108 = 1/Rap;
    a1010 = -((1/Rap)+(1/Rvp));
    
    % A =    [Vla   Vlv  Qa  Pas  Pvs  Vra   Vrv    Pap  Qp   Pvp]
    A = zeros(10,10);
    A(1,:) = [a11 , a12,   0,   0,   0,   0,    0,    0,   0,  a110];         % Vla  
    A(2,:) = [a21 , a22,  -1,   0,   0,   0,    0,    0,   0,     0];         % Vlv
    A(3,:) = [ 0  , a32, a33,  -1,   0,   0,    0,    0,   0,     0]*(Da/Ls); % Qa
    A(4,:) = [ 0  ,  0 ,   1, a44, a45,   0,    0,    0,   0,     0]/Cas;     % Pas
    A(5,:) = [ 0  ,  0 ,   0, a54, a55, a56,    0,    0,   0,     0]/Cvs;     % Pvs
    A(6,:) = [ 0  ,  0 ,   0,   0, a65, a66,  a67,    0,   0,     0];         % Vra
    A(7,:) = [ 0  ,  0 ,   0,   0,   0, a76,  a77,    0,  -1,     0];         % Vrv
    A(8,:) = [ 0  ,  0 ,   0,   0,   0,   0,    0,  a88,   1,  a810]/Cap;     % Pap
    A(9,:) = [ 0  ,  0 ,   0,   0,   0,   0,  a97,   -1, a99,     0]*(Dp/Lp); % Qp
    A(10,:)= [a101,  0 ,   0,   0,   0,   0,    0, a108,   0, a1010]/Cvp;     % Pvp
 
%     B = zeros(10,1);
    B = [    (-((1/Rvp)+(Dm/Rm))*(-Ea(i)*Vla0)) + ((Dm/Rm)*(- Elv(i)*Vlv0));
                                        (Dm/Rm)*(- Ea(i)*Vla0 + Elv(i)*Vlv0);
                                                       (1/Ls)*(- Elv(i)*Vlv0);
                                                                            0;
                                                 (1/(Rvs*Cvs))*(- Ea(i)*Vra0);
            (-((1/Rvs)+(Dt/Rt))*(- Ea(i)*Vra0)) + ((Dt/Rt)*(- Erv(i)*Vrv0));
                                        (Dt/Rt)*(- Ea(i)*Vra0 + Erv(i)*Vrv0);
                                                                            0;
                                                       (1/Lp)*(- Erv(i)*Vrv0);
                                                 (1/(Rvp*Cvp))*(- Ea(i)*Vla0)];
    
    if Da == 0
        x(3) = 0;
        B(3) = 0;
    end
    if Dp == 0
        x(9) = 0;
        B(9) = 0;
    end

    x = runkut4(passo, A, x, B);
    
    Vla(i+1) = x(1);
    Vlv(i+1) = x(2);
    Qa(i+1) = x(3);
    Pas(i+1) = x(4);
    Pvs(i+1) = x(5);
    Vra(i+1) = x(6);
    Vrv(i+1) = x(7);
    Pap(i+1) = x(8);
    Qp(i+1) = x(9);    
    Pvp(i+1) = x(10);
    
    Pla(i+1) = Ea(i+1)*(Vla(i+1) - Vla0);
    Plv(i+1) = Elv(i+1)*(Vlv(i+1) - Vlv0);
    Pra(i+1) = Ea(i+1)*(Vra(i+1) - Vra0);
    Prv(i+1) = Erv(i+1)*(Vrv(i+1) - Vrv0);
    
    Qa_dot = [0 diff(Qa)];
    Pao(i+1) = Rcs*Qa(i+1) + Qa_dot(end) + Pas(i+1);
       
end

% REDV(end)
% RESV(end)
% LEDV(end) 
% LESV(end)

SV_right = (REDV(end) - RESV(end));
SV_left = (LEDV(end) - LESV(end));

CO_right = HR*0.001*SV_right;
CO_left = HR*0.001*SV_left;

ejectionFraction_right = 100*(SV_right/REDV(end));
ejectionFraction_left = 100*(SV_left/LEDV(end));

figure()
subplot(8,2,1)
plot(t,Vla)
title('Vla')
subplot(8,2,2)
plot(t,Vlv)
title('Vlv')
subplot(8,2,3)
plot(t,Qa)
title('Qa')
subplot(8,2,4)
plot(t,Pas)
title('Pas')
subplot(8,2,5)
plot(t,Pvs)
title('Pvs')
subplot(8,2,6)
plot(t,Vra)
title('Vra')
subplot(8,2,7)
plot(t,Vrv)
title('Vrv')
subplot(8,2,8)
plot(t,Pap)
title('Pap')
subplot(8,2,9)
plot(t,Qp)
title('Qp')
subplot(8,2,10)
plot(t,Pvp)
title('Pvp')

subplot(8,2,11)
plot(t,Pla)
title('Pla')
subplot(8,2,12)
plot(t,Plv)
title('Plv')
subplot(8,2,13)
plot(t,Pra)
title('Pra')
subplot(8,2,14)
plot(t,Prv)
title('Prv')

subplot(8,2,15)
plot(t,Pla,t,Plv,t,Pao)
title('Left heart')
subplot(8,2,16)
plot(t,Pra,t,Prv,t,Pap)
title('Right heart')
