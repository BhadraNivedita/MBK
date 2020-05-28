
% Writing down th equation for many body Kapitza pendula



function dy = frekon(t,y)

%% Setting the parameters
nn=100; mm=nn/2;
dy=zeros(nn,1);
Lam=.001;KK=1;
g0=1.0; g1=10; gam=5.0;

%% equations for the many body Kapitzapendula
for i=1:mm
    dy(i)=Lam*KK*y(mm+i);
end


dy(mm+1)=(Lam/KK)*(sin(y(mm)-y(1))-sin(y(1)-y(2)))-(g0+g1*cos(gam*t))*sin(y(1))/(KK*Lam);

for i=2:mm-1
    dy(mm+i)=(Lam/KK)*(sin(y(i-1)-y(i))-sin(y(i)-y(i+1)))-(g0+g1*cos(gam*t))*sin(y(i))/(KK*Lam);
end

dy(nn) =(Lam/KK)*(sin(y(mm-1)-y(mm))-sin(y(mm)-y(1)))-(g0+g1*cos(gam*t))*sin(y(mm))/(KK*Lam);

%for i=1:mm
 %   dy(i)=Lam*y(mm+i);
%end

%dy(mm+1)=Lam*(sin(y(mm)-y(1))-sin(y(1)-y(2)))-(g0+g1*cos(gam*t))*sin(y(1))/Lam;

%for i=2:mm-1
 %   dy(mm+i)=Lam*(sin(y(mm-1)-y(mm))-sin(y(mm)-y(1)))-(g0+g1*cos(gam*t))*sin(y(mm))/Lam;
%end
%dy(nn) =Lam*(sin(y(mm-1)-y(mm))-sin(y(mm)-y(1)))-(g0+g1*cos(gam*t))*sin(y(mm))/Lam;
%end
%run the code manybodykapitzarun.m
