function DDE = porousHess(x,xhat,f,gamma)
%%  Evaluate the Hessian of Power energy for the FE flow map

if gamma==1,    DDE = fickHess(x,xhat,f);
else
    
N  = length(x)-1;
ddE = @(dx,dxhat,f) gallery('tridiag',...
      -gamma*(gamma-1)*dxhat(2:N-1).*(f(2:N-1).^gamma)./(dx(2:N-1).^(gamma+1)),...
       gamma*(gamma-1)*dxhat(1:N-1).*(f(1:N-1).^gamma)./(dx(1:N-1).^(gamma+1))...
      +gamma*(gamma-1)*dxhat(2:N)  .*(f(2:N).^gamma)  ./(dx(2:N).^(gamma+1)),...
      -gamma*(gamma-1)*dxhat(2:N-1).*(f(2:N-1).^gamma)./(dx(2:N-1).^(gamma+1)));

xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy Hessian element-wise via quadrature (Simpson's rule)
DDE = (  ddE(dx,dxhat,f(xhat(1:N)))...
      +4*ddE(dx,dxhat,f(xhatm))...
      	+ddE(dx,dxhat,f(xhat(2:N+1)))   )/6;

clear xhatm dx dxhat;

end