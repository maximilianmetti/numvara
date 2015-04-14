function [ff,xstag,t,E,D,time] = energyform(f,gamma,phi,dphi,ddphi,T,N,M)
%%   Computes approximate solution to the heat eqn via energy formulation
%   f = @(x) exp(-20*x.*x);
%   N = 201;
%   M = 100;
%   T = 1;
%   gamma = 1;
% phi   = @(x) zeros(size(x,1),1);
% dphi  = @(x) zeros(size(x,1),1);
% ddphi = @(x) zeros(size(x,1),1);
% phi   = @(x) barrierpotential(x,.001,100,.6);
% dphi  = @(x) dbarrierpotential(x,.001,100,.6);
% ddphi = @(x) zeros(size(x,1),1);
% b   = @(x) (1.0000001-x.*x);
% db  = @(x) -2*x;
% ddb = @(x) -2*ones(size(x,1),1);




%% Define discrete flow map space and initial condition
%%  Initialize Grid and variables
tic;
dt = T/M;
x  = (-1:(2/N):1)';
xhat = x;
t  = zeros(M,1);
E  = zeros(M,1);
D  = zeros(M,1);
Moment = zeros(M,3);
xstag = zeros(M,2*N);
ff = zeros(M,2*N);





%%  Evaluate Energy and Dissipation
%E(1) = generalFickEnergy(x,xhat,b,f);
E(1) = porousEnergy(x,xhat,f,gamma)+potentialEnergy(x,xhat,f,phi);%fickEnergy(x,xhat,f);%
D(1) = darcyDiss(x,x,xhat,dt,f);
Moment(1,1) = moment(0,x,xhat,f);
Moment(1,2) = moment(1,x,xhat,f);
Moment(1,3) = moment(2,x,xhat,f);
display(['Moments are ', num2str(Moment(1,1)),' ',num2str(Moment(1,2)),' ',num2str(Moment(1,3))]);




%%  Plot Initial Value
figure(2);clf;
dx     = x(2:N+1)-x(1:N);
xstag  = zeros(2*N,1);
fstag  = zeros(2*N,1);
xstag(2*(1:N)-1) = x(1:N);
xstag(2*(1:N))   = x(2:N+1);
fstag(2*(1:N)-1) = f((-1:(2/N):1-2/N)')./(dx*N/2);
fstag(2*(1:N))   = f((-1+2/N:(2/N):1)')./(dx*N/2);
%xSTAG(1,:) = xstag';
%fSTAG(1,:) = fstag';
subplot(2,1,1);
plot(xstag,fstag,x,-.05*ones(size(x,1),1),'r+',x,zeros(size(x,1),1),':');
title(['Density at t=',num2str(t(1)),' with N=',num2str(length(x))]);
axis([-1 1 -.06 1.1]);
subplot(2,1,2);
plot(t(1),E(1),t(1),D(1));
axis([0 T (E(1)-1) E(1)]);
title('Energy & Dissipation');
legend('Energy','Dissipation');
pause();




%% Solver for energy minimizing flow map
steady = false;   % flag for steady state
for m=1:M
    t(m+1)= t(m)+dt;
    xprev = x;
    
    %W  = @(p,pprev,phat) fisherWrightEnergy(p,phat,f) + .5*dt*darcyDiss(p,pprev,phat,dt,f);
    W  = @(p,pprev,phat) porousEnergy(p,phat,f,gamma) + potentialEnergy(p,phat,f,phi)...
                + .5*dt*darcyDiss(p,pprev,phat,dt,f);
    %DW = @(p,pprev,phat) fisherWrightVar(p,phat,f) + .5*dt*darcyVar(p,pprev,phat,dt,f);
    DW = @(p,pprev,phat) porousVar(p,phat,f,gamma) + potentialVar(p,phat,f,dphi)...
                + .5*dt*darcyVar(p,pprev,phat,dt,f);
    %DDW= @(p,pprev,phat) fisherWrightHess(p,phat,f) + .5*dt*darcyHess(p,pprev,phat,dt,f);
    DDW= @(p,pprev,phat) porousHess(p,phat,f,gamma) + potentialHess(p,phat,f,ddphi)...
                + .5*dt*darcyHess(p,pprev,phat,dt,f);
            
	display('Solving');
    for cycle=1:1
        [x,xprev,xhat] = Newton(x,xprev,xhat,W,DW,DDW,1e-5,30,true);%nonlinconjgrad(x,W,DW,.0001,false); %
        display('Completed Full Newton');
        [x,xprev,xhat] = nonlinearGSnewton(x,xprev,xhat,W,DW,DDW,0.1,0);
        display('Completed G-S Newton');
    end
    N  = length(x)-1;
    
    
    %% Compute energies
    %E(m+1) = fisherWrightEnergy(x,xhat,f);%
    E(m+1) = porousEnergy(x,xhat,f,gamma) + potentialEnergy(x,xhat,f,phi);%fickEnergy(x,xhat,f);%
    D(m+1) = darcyDiss(x,xprev,xhat,dt,f);
    Moment(m+1,1) = moment(0,x,xhat,f);
    Moment(m+1,2) = moment(1,x,xhat,f);
    Moment(m+1,3) = moment(2,x,xhat,f);
    display(['Moments are ', num2str(Moment(m+1,1)),' ',num2str(Moment(m+1,2)),' ',num2str(Moment(m+1,3))]);
%     if abs(D(m+1)) < 1e-5 && ~steady
%         steady=true;
%         display('At steady state');
%     end
    




    %%  Plot Solution
    dx     = x(2:N+1)-x(1:N);
    dxhat  = xhat(2:N+1)-xhat(1:N);
    xstag  = zeros(2*N,1);
    fstag  = zeros(2*N,1);
    xstag(2*(1:N)-1) = x(1:N);
    xstag(2*(1:N))   = x(2:N+1);
    fstag(2*(1:N)-1) = f(xhat(1:N  ))./(dx./dxhat);%f((-1:(2/N):1-2/N)')./(dx*N/2);
    fstag(2*(1:N))   = f(xhat(2:N+1))./(dx./dxhat);%f((-1+2/N:(2/N):1)')./(dx*N/2);
    %xSTAG(m+1,:) = xstag';
    %fSTAG(m+1,:) = fstag';
    yy = (x(1:N)+x(2:N+1)) / 2;
    ff = zeros(2*N,1);
    ff(2*(1:N)-1) = (f(xhat(1:N))./(dx./dxhat) + f(xhat(2:N+1))./(dx./dxhat)) / 2;
    ff(2*(1:N))   = (f(xhat(1:N))./(dx./dxhat) + f(xhat(2:N+1))./(dx./dxhat)) / 2;
    subplot(2,1,1);
    %plot(xstag,fstag,x,-.05*ones(size(x,1),1),'r+',x,zeros(size(x,1),1),':');
    plot(xstag,ff,x,-.05*ones(size(x,1),1),'r+',x,zeros(size(x,1),1),':');
    title(['Density at t=',num2str(t(m+1)),' with N=',num2str(length(x))]);
    axis([-1 1 -.06 1.1]);
    subplot(2,1,2);
    plot(t(1:m+1),E(1:m+1),t(1:m+1),D(1:m+1));
    axis([0 T E(m+1) E(1)]);
    title('Energy & Dissipation');
    legend('Energy','Dissipation');
    pause(.01);
end
clear dx xprev xstag fstag
time=toc;