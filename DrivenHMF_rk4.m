% Simulation of Driven HMF model

function dy=DrivenHMF_rk4(t,y,N,K)

%% Setting the parameters
N=80; % No of coupled pendulum
dy=zeros(2*N,1);
g0=0; g1=0; omega=0;
gt = g0+g1*cos(omega*t);

%% equations for the many body Kapitza pendula
for i=1:N
    dy(i) = y(N+i);
end

for i=1:N
    dy(N+i) = gt*sin(y(i));
    for j=1:N
        dy(N+i) = dy(N+i) + (K/N)*sin(y(i)-y(j));
    end
end

end
