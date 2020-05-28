% Keep the file HMF.m in the same folder

clear all; clc;

N=5; % No of coupled pendulum
combination=1; average=100;
tspan=0:0.01:1000; Uspan=.5;

K=1;

H=zeros(length(tspan),1);
T=zeros(length(tspan),1);
M=zeros(length(tspan),1);

for l=1:length(Uspan)

U=Uspan(l);

for k=1:combination
    
%% Initial positions and momenta

pos = 0.01*(rand(N,1)-0.5); % Random initial positions

PE=0;
for i=1:N
    for j=1:N
        PE = PE + cos(pos(i)-pos(j));
    end
end
PE = PE*(K/(2*N));

mom = sqrt(2*(U*N-PE)/N)*(-1).^(floor(2*rand(N,1))); % Initial momenta

%% ODE45

[t,x]=ode45(@(t,y) HMF(t,y,N,K),tspan,[pos,mom]);

%% Calculating potential energy, kinetic energy and total energy in each time step
for m=1:length(tspan)

KE=0; PE=0;

for j=(N+1):2*N
    KE = KE + x(m,j)^2;
end
KE = KE/2;

for i=1:N
    for j=1:N
        PE = PE + cos(x(m,i)-x(m,j));
    end
end
PE = PE*(K/(2*N));

H(m) = KE + PE;
T(m) = 2*KE;

%% Calculating Magnetization (Order parameter)

Mx=0; My=0;
for j=1:N
Mx = Mx + cos(x(m,j));
My = My + sin(x(m,j));
end
M(m) = sqrt(Mx^2+My^2)/N;
end

end

end