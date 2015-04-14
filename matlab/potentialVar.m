function DE = potentialVar(x,xhat,f,dphi)
%%  Evaluate the variation of Fick's energy for the FE flow map
N  = length(x)-1;
dE = @(y,dxhat,f,dphi)    y*f(1:N-1).*dphi(1:N-1).*dxhat(1:N-1)...
                     +(1-y)*f(2:N).*dphi(2:N).*dxhat(2:N);

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(x(1:N)+x(2:N+1));
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy Variation element-wise via quadrature (Simpson's rule)
DE = (   dE(0, dxhat,f(xhat(1:N)),  dphi(x(1:N))) ...
      +4*dE(.5,dxhat,f(xhatm(1:N)), dphi(xm(1:N)))...
        +dE(1, dxhat,f(xhat(2:N+1)),dphi(x(2:N+1)))   )/6;
    
%%  Fixed Boundary
DE = [0;DE;0];

clear xhatm xm dxhat;