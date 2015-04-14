function [ff,xstag,t,E,D,time] = radialEnergyForm(f,gamma,phi,dphi,ddphi,T,N,M)
%%   Computes approximate solution to the heat eqn via energy formulation
%   f = @(x) exp(-20*x.*x);
%   N = 200;
%   M = 100;
%   T = 1;
%   gamma = 1;
% phi   = @(x) zeros(size(x,1),1);
% dphi  = @(x) zeros(size(x,1),1);
% ddphi = @(x) zeros(size(x,1),1);
% phi   = @(x) barrierpotential(x,.01,100,.6);
% dphi  = @(x) dbarrierpotential(x,.01,100,.6);
% ddphi = @(x) zeros(size(x,1),1);
% b   = @(x) (1.0000001-x.*x);
% db  = @(x) -2*x;
% ddb = @(x) -2*ones(size(x,1),1);


%%  Initialize Grid and  variables
tic;
dt = T/M;
x  = (1/(N+1):(1/(N+1)):1)';
xhat = x;
xprev= x;
t  = zeros(M,1);
E  = zeros(M,1);
D  = zeros(M,1);
xstag = zeros(M,2*N);
ff = zeros(M,2*N);

%%  Evaluate Energy and Dissipation
E(1) = radialFickEnergy(x,xhat,f);
D(1) = radialDarcyDiss(x,xprev,xhat,dt,f);


%%  Plot Initial Value
clf;
dx     = x(2:N+1)-x(1:N);
xstag  = zeros(2*N,1);
fstag  = zeros(2*N,1);
xstag(2*(1:N)-1) = x(1:N);
xstag(2*(1:N))   = x(2:N+1);
fstag(2*(1:N)-1) = f(x(1:N))./(dx*N);
fstag(2*(1:N))   = f(x(2:N+1))./(dx*N);
subplot(2,1,1);
plot(xstag,fstag,x,-.05*ones(size(x,1),1),'r+',x,zeros(size(x,1),1),':');
title(['Density at t=',num2str(t(1)),' with N=',num2str(length(x))]);
axis([0 1 -.06 1.1]);
subplot(2,1,2);
plot(t(1),E(1),t(1),D(1));
axis([0 T E(1)-1 E(1)]);
title('Energy & Dissipation');
legend('Energy','Dissipation');
pause();



for m=1:M
    t(m+1)= t(m)+dt;
    xprev = x;
    
    W  = @(p,pprev,phat) radialFickEnergy(p,phat,f) + .5*dt*radialDarcyDiss(p,pprev,phat,dt,f);
    DW = @(p,pprev,phat) radialFickVar(p,phat,f) + .5*dt*radialDarcyVar(p,pprev,phat,dt,f);
    DDW= @(p,pprev,phat) radialFickHess(p,phat,f) + .5*dt*radialDarcyHess(p,pprev,phat,dt,f);
            
    %[x,xprev,xhat] = Newton(x,xprev,xhat,W,DW,DDW,1e-3,50,true);%nonlinconjgrad(x,W,DW,.0001,false); %
    [x,xprev,xhat] = nonlinearGSbisection(x,xprev,xhat,W,DW,0.15,100,3);
    %[x,xprev,xhat] = nonlinearGSnewton(x,xprev,xhat,W,DW,DDW,0.1,5);
    N  = length(x)-1;
    
    E(m+1) = radialFickEnergy(x,xhat,f);
    D(m+1) = radialDarcyDiss(x,xprev,xhat,dt,f);


    %%  Plot Solution
    dx     = x(2:N+1)-x(1:N);
    dxhat  = xhat(2:N+1)-xhat(1:N);
    xstag  = zeros(2*N,1);
    fstag  = zeros(2*N,1);
    xstag(2*(1:N)-1) = x(1:N);
    xstag(2*(1:N))   = x(2:N+1);
    fstag(2*(1:N)-1) = f(xhat(1:N  ))./(dx./dxhat);%f((-1:(2/N):1-2/N)')./(dx*N/2);
    fstag(2*(1:N))   = f(xhat(2:N+1))./(dx./dxhat);%f((-1+2/N:(2/N):1)')./(dx*N/2);
    yy = (x(1:N)+x(2:N+1)) / 2;
    ff = zeros(2*N,1);
    ff(2*(1:N)-1) = (f(xhat(1:N))./(dx./dxhat) + f(xhat(2:N+1))./(dx./dxhat)) / 2;
    ff(2*(1:N))   = (f(xhat(1:N))./(dx./dxhat) + f(xhat(2:N+1))./(dx./dxhat)) / 2;
    subplot(2,1,1);
    plot(xstag,ff,x,-.05*ones(size(x,1),1),'r+',x,zeros(size(x,1),1),':');
    title(['Density at t=',num2str(t(m+1)),' with N=',num2str(length(x))]);
    axis([0 1 -.06 1.1]);
    subplot(2,1,2);
    plot(t(1:m+1),E(1:m+1),t(1:m+1),D(1:m+1));
    axis([0 T E(m+1) E(1)]);
    title('Energy & Dissipation');
    legend('Energy','Dissipation');
    pause(.01);
end


clear dx xprev xstag fstag
time=toc;