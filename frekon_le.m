function f=frekon_le(t,X);

% Frenkel-Kontorova equation
% phi_i_dot = Lambda*K*P_i
% P_i_dot = (lambda/K)*(sin(phi_{i-1}-phi_i)-sin(phi_i-phi_{i-1}))
%              - ((g_0+g_1*cos(gamma*t))/K*Lambda)*sin(phi_i)

% Values of parameters

nn=100; mm=nn/2;
Lam=50; KK=12;
g0=1; g1=15; gam=10.0;

for i=1:nn
    for j=1:nn
    Y(i,j) = X(nn*j+i);
    end
end

f=zeros(nn^2,1);

gt = g0+g1*cos(gam*t);

% Flow equation

for i=1:mm
    f(i)=Lam*KK*X(mm+i);
end

f(mm+1)=(Lam/KK)*(sin(X(mm)-X(1))-sin(X(1)-X(2)))-gt*sin(X(1))/(KK*Lam);

for i=2:mm-1
    f(mm+i)=(Lam/KK)*(sin(X(i-1)-X(i))-sin(X(i)-X(i+1)))-gt*sin(X(i))/(KK*Lam);
end

f(nn) =(Lam/KK)*(sin(X(mm-1)-X(mm))-sin(X(mm)-X(1)))-gt*sin(X(mm))/(KK*Lam);

% Linearized system

Jac = 0.0;

for i=1:mm
    Jac(i,mm+i) = Lam*KK;
end

Jac(mm+1,1) = -(Lam/KK)*(cos(X(mm)-X(1))+cos(X(1)-X(2)))-(gt/(KK*Lam))*cos(X(1));
Jac(mm+1,2) = (Lam/KK)*cos(X(1)-X(2));
Jac(mm+1,mm) = (Lam/KK)*cos(X(mm)-X(1));

for i=2:mm-1
    Jac(mm+i,i-1) = (Lam/KK)*cos(X(i-1)-X(i));
    Jac(mm+i,i) = -(Lam/KK)*(cos(X(i-1)-X(i))+cos(X(i)-X(i+1)))-(gt/(KK*Lam))*cos(X(i));
    Jac(mm+i,i+1) = (Lam/KK)*cos(X(i)-X(i+1));
end

Jac(nn,1) = (Lam/KK)*cos(X(mm)-X(1));
Jac(nn,mm-1) = (Lam/KK)*cos(X(mm-1)-X(mm));
Jac(nn,mm) = -(Lam/KK)*(cos(X(mm-1)-X(mm))+cos(X(mm)-X(1)))-(gt/(KK*Lam))*cos(X(mm));
  
% Variational equation

f(nn+1:nn+nn^2)=Jac*Y;

% To run: "kuramoto_le_run.m"