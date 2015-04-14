function DDD = radialDarcyHess(x,xprev,xhat,dt,f)
%%  Evaluate the Hessian of Darcy's dissipation for the FE flow map
N  = length(x)-1;
ddD = @(x,xp,xhat,y,dxhat,dt,f) (2/(dt*dt))*gallery('tridiag',...
                 (1-y)*y*dxhat(2:N-1).*xhat(2:N-1).*f(2:N-1),...
                     y*y*dxhat(1:N-1).*xhat(1:N-1).*f(1:N-1)...
            +(1-y)*(1-y)*dxhat(2:N).*xhat(2:N).*f(2:N),...
                 (1-y)*y*dxhat(2:N-1).*xhat(2:N-1).*f(2:N-1));

xm    = .5*(x(1:N)+x(2:N+1));
xpm   = .5*(xprev(1:N)+xprev(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dxhat = xhat(2:N+1)-xhat(1:N);

%%  Evaluate dissipation Hessian element-wise via quadrature (Simpson's)
DDD = (  ddD(x(1:N),xprev(1:N),xhat(1:N),0,dxhat,dt,f(xhat(1:N))) ...
      +4*ddD(xm,xpm,xhatm,.5,dxhat,dt,f(xhatm))...
      +  ddD(x(2:N+1),xprev(2:N+1),xhat(2:N+1),1,dxhat,dt,f(xhat(2:N+1))) )/6;      
  
clear xm xpm xhatm dxhat;