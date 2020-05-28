%Keep the file DrivenHMF_rk4.m in the same folder

clear all; close all; clc;

N=80; % No of coupled pendulum
combination=10; average=1000;
tspan=0:0.01:500; Kspan=0:1:4;

r=zeros(length(Kspan),1);

for l=1:length(Kspan)

K=Kspan(l);

for k=1:combination

w = random('Normal',0,1,1,N); % Random initial conditions
[t,x]=ode45(@(t,y) DrivenHMF_rk4(t,y,N,K),tspan,[w,w]);

for i=(length(tspan)-average):length(tspan)
    rx=0; ry=0;
    for j=1:N
        rx = rx + cos(x(i,j));
        ry = ry + sin(x(i,j));
    end
    r(l) = r(l) + sqrt(rx^2+ry^2)/N;
end

end

r(l) = r(l)/(combination*average);

end

plot(Kspan,r)