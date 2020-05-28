%Keep the file frekon_le.m in the same folder

tic

clear all; close all; clc
mm=2; nn=mm*2;

tt=50;
pos=ones(mm,1)*pi+0.01; mom=ones(mm,1)*0.001;

[T,Res]=lyapunov(nn,@MBK_le,@ode45,0,0.05,tt,[pos,mom]);

le=[T,Res];

NN=length(Res);

save('MBK_le.dat','le','-ascii')






 plot(T,Res);

toc

