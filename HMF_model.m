% Simulation of massive Kuramoto model/HMF model
tic

clear all ;close all; clc; 
N=100; %No of coupled pendulum 
T=1000;% Total time
tau=0.01;% Timestep

combination=1;

%% Defining the coupling constant K
K=1;

%Creating vector for theta and thetadot
thetadot=zeros(N,1);
theta=zeros(N,1);
en=zeros(N,1);
phi=zeros(N,1);
r=zeros(N,1);
EE=zeros(T,1);

for run=1:combination % Number of combination
 
    
% Initial conditions with normal distribution for w
sigN=1; 
w = random('Normal',0,sigN,1,N);% Random on initial conditions

    en=10;

    % Attempting to calculate r (Phase coherence between 0 and 1) 
    % and phi (Average angle) for the first steps of time
    rx=0;
    ry=0;
    phi(1)=0;
    
    % Initialize theta and thetadot
    for i=1:N
        theta(1,i)=w(i);
        %thetadot(1,i)=w(i);
    end
    
    % Calculating average angle/mean phi for N oscillators
    for i=1:N
        phi(1) = phi(1) + (1/N)*theta(1,i); % Calculation mean angle phi
        rx=rx+(1/N)*cos(theta(1,i)); % Sum of mean x-part of theta
        ry=ry+(1/N)*sin(theta(1,i)); % Sum of mean y-part of theta
    end
    
    r(1) = sqrt(rx*rx + ry*ry); %defining the phase coherence between 0 and 1
    r(2) = r(1);
    phi(2)=phi(1);
    
    for i=1:N
      thetadot(1,i)=sqrt(2*(en-r(i).*cos(theta(1,i)-phi(i))));
    end
 
    en1=0.5*thetadot(1,i).^2+r(1)*K*cos(theta(1,i)-phi(i));

    %% The main loop
    for t=2:T   % Initial conditions each timestep in the loop
        rx=0;
        ry=0;
        phi(t+1)=0;
       
        entot=0;
        %% Simulating the main equation
        
        for i=1:N
           
            theta(t,i)=theta(t-1,i)+tau*thetadot(t-1,i);
            thetadot(t,i)=thetadot(t-1,i)-tau*(theta(t-1,i)+K*r(t)*sin(theta(t-1,i)-phi(t))...
                );
            
            rx=rx+(1/N)*cos(theta(t,i)); % Sum of mean x-part of theta
            ry=ry+(1/N)*sin(theta(t,i)); % Sum of mean y-part of theta
            
            en(t,i)=0.5*thetadot(t-1,i).^2+r(t)*K*cos(theta(t-1,i)-phi(t));
           
            phi(t+1) = phi(t+1) + (1/N)*theta(t,i); % Calculating mean angle phi for next step of time
            
            entot=entot+en(t,i);
            
            
        end

    
        EE(t+1)=entot;
    
        r(t+1) = sqrt(rx*rx + ry*ry); % Calculating total mean radius
        
    end

   

end





toc
