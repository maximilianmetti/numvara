function DE = porousVar(x,xhat,f,gamma)
%%  Evaluate the variation of Power energy for the FE flow map

if gamma==1,    DE = fickVar(x,xhat,f);
else
    
N  = length(x)-1;
dE = @(dx,dxhat,f) -(gamma-1)*(f(1:N-1).^gamma).*dxhat(1:N-1)./(dx(1:N-1).^gamma)...
                   +(gamma-1)*(f(2:N).^gamma).*dxhat(2:N)./(dx(2:N).^gamma);

xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy Variation element-wise via quadrature (Simpson's rule)
DE = (   dE(dx,dxhat,f(xhat(1:N) ))...
      +4*dE(dx,dxhat,f(xhatm(1:N)))...
        +dE(dx,dxhat,f(xhat(2:N+1)))   )/6;
    
%%  Fixed Boundary
DE = [0;DE;0];

clear xhatm dx dxhat;

end