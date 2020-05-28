% Simulation of massive Kuramoto model/HMF model
tic

clear all ;close all; clc; 
N=80; %No of coupled pendulum 
T=500;% Total time
tau=0.01;% Timestep

combination=10;

%% Defining the coupling constant K
Klow = 0; Kup = 3.5; Kstep = 0.01;
totsteps = round((Kup-Klow)/Kstep);
Ktot=zeros(totsteps,1);
Ktot(1) = Klow;


for i=2:totsteps;
    Ktot(i) =  Ktot(i-1)+Kstep;% Loop for K
end


%Creating vector for theta and thetadot
thetadot=zeros(N,1);
theta=zeros(N,1);
phi=zeros(N,1);
r=zeros(N,1);
rtot=zeros(combination,totsteps);
%g=.1;g1=1;

g=0;g1=0;
omega=0; % Frequency of applied drive

for run=1:combination % Number of combination

  
    
% Initial conditions with normal distribution for w
sigN=1; 
w = random('Normal',0,sigN,1,N);% Random on initial conditions


for Steps=1:totsteps     %Begining of the coupling(K) loop;
    
    K=Klow+ Steps*Kstep;% Setting the value of K each lap
    
 % Initialize theta and thetadot
    for i=1:N
        theta(1,i)=w(i);
        thetadot(1,i)=w(i);
    end
    
   

    % Attempting to calculate r (Phase coherence between 0 and 1) 
    % and phi (Average angle) for the first steps of time
    rx=0;
    ry=0;
    phi(1)=0;
    
    % Calculating average angle/mean phi for N oscillators
    for i=1:N
        phi(1) = phi(1) + (1/N)*theta(1,i); % Calculation mean angle phi
        rx=rx+(1/N)*cos(theta(1,i)); % Sum of mean x-part of theta
        ry=ry+(1/N)*sin(theta(1,i)); % Sum of mean y-part of theta
    end
    
    r(1) = sqrt(rx*rx + ry*ry); %defining the phase coherence between 0 and 1
    r(2) = r(1);
    phi(2)=phi(1);

    %% The main loop
    for t=2:T   % Initial conditions each timestep in the loop
        rx=0;
        ry=0;
        phi(t+1)=0;
       
        
        %% Simulating the main equation
        
        for i=1:N
           
            theta(t,i)=theta(t-1,i)+tau*thetadot(t-1,i);
            thetadot(t,i)=thetadot(t-1,i)-tau*(theta(t-1,i)+K*r(t)*sin(theta(t-1,i)-phi(t))...
                +(g+g1*cos(omega*t))*sin(theta(t-1,i)));
            
            rx=rx+(1/N)*cos(theta(t,i)); % Sum of mean x-part of theta
            ry=ry+(1/N)*sin(theta(t,i)); % Sum of mean y-part of theta
           
            phi(t+1) = phi(t+1) + (1/N)*theta(t,i); % Calculating mean angle phi for next step of time
        end
    
        r(t+1) = sqrt(rx*rx + ry*ry); % Calculating total mean radius
        
    end
    %% Assigning the final value of r at each K to rtot
    rtot(run,Steps) = mean(r(T-100:T+1));
   % storew(Steps,:)=w + K*r(t)*sin(phi(t)-theta(t-1,:));
end

end


%% Plotting the radius r (coherence) over K as bifurcation diagram
figure(1)
hold off
plot(Ktot,mean(rtot),'-o','Linewidth',1)%b
hlx=xlabel('Coupling strength: K');
hly=ylabel('Average Coherence: r');
% Plotting Kuramotos analytical results
% only with symbolic toolbox
%hold on;
%plot(Ktot2,ar,'k:')
axis([Klow Kup -0.05 1])


print -depsc -painters g0=0g1=0w0Kura.eps
%figure(2)


toc
