function DD = darcyVar(x,xprev,xhat,dt,f)
%%  Evaluate the variation Darcy's dissipation for the FE flow map
N  = length(x)-1;
dD = @(x,xp,y,dxhat,dt,f) 2*(  y*dxhat(1:N-1).*f(1:N-1).*(x(1:N-1)-xp(1:N-1))...
                          +(1-y)*dxhat(2:N).*f(2:N)  .*(x(2:N)-xp(2:N)) )/(dt*dt);

xm    = .5*(x(1:N)+x(2:N+1));
xpm   = .5*(xprev(1:N)+xprev(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate dissipation Variation element-wise via quadrature (Simpson's)
DD = (   dD(x(1:N),xprev(1:N),0,dxhat,dt,f(xhat(1:N))) ...
      +4*dD(xm,xpm,.5,dxhat,dt,f(xhatm))...
        +dD(x(2:N+1),xprev(2:N+1),1,dxhat,dt,f(xhat(2:N+1)))  )/6;   
    
%%  Fixed Boundary
DD = [0;DD;0];
  
clear xm xpm xhatm dxhat ;