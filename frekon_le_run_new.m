
tic

clear all;% close all; clc
nn=100; mm=nn/2;

fileID=fopen('frekon_le_K=12_Lambda=1new.dat','w');

pos=ones(mm,1)*pi+0.01; mom=ones(mm,1)*0.001;

[T,Res]=lyapunov(nn,@frekon_le,@ode45,0,0.05,5,[pos,mom]);

N=length(Res);

for ii=1:NN
    
 fprintf(fileID,'%f\n',Res(ii));
 
end


toc

fclose(fileID);