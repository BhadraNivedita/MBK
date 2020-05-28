
%Keep the file frekon.m in the sam folder

clear all; close all; clc;
nn=100;mm=nn/2;

%% integrating the manybody kapitza pendula
pos=ones(mm,1)*pi+0.001; mom=ones(mm,1)*0.01;

[t,x]=ode45('frekon',(0:0.01:100),[pos,mom]);

pos=x(:,1:mm); mom=x(:,mm+1:nn);

%% averaging all the oscillators
avpos=mean(pos');%averaging the position
avmom=mean(mom');%averaging the momentum

%% Fourier spectrum
Tstep=0.01;
N=length(avpos);
Fs = 1/Tstep;
xdft = fft(avpos);
xdft = xdft(1:(N-1)/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(avpos):Fs/2;

%% Plotting 

subplot(2,2,1)
plot(avpos)
axis tight

xlabel('time','Interpreter','LaTex','Fontsize',20)
ylabel('$\phi$','Interpreter','LaTex','Fontsize',20)


subplot(2,2,2)
plot(avmom)
axis tight

xlabel('time','Interpreter','LaTex','Fontsize',20)
ylabel('$\dot{\phi}$','Interpreter','LaTex','Fontsize',20)


subplot(2,2,3)
plot(avpos,avmom)
axis tight
xlabel('$\phi$','Interpreter','LaTex','Fontsize',20)
ylabel('$\dot{\phi}$','Interpreter','LaTex','Fontsize',20)


subplot(2,2,4)

plot(2*pi*freq,10*log10(psdx))

axis([0 50 -100 10])

xlabel('time','Interpreter','LaTex','Fontsize',20)
ylabel('10$\log$10(psdx)','Interpreter','LaTex','Fontsize',20)

print -depsc -painters ManybodykapitzaLam80.eps