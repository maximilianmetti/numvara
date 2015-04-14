function DE = radialFickVar(x,xhat,f)
%%  Evaluate the variation of Fick's energy for the FE flow map
N  = length(x)-1;
dE = @(x,xhat,dx,dxhat,f,bas) xhat(1:N-1).*f(1:N-1).*dxhat(1:N-1).*(-(1/dx(1:N-1))'-(bas/x(1:N-1))')...
                            + xhat(2:N).*f(2:N).*dxhat(2:N).*((1/dx(2:N))'-((1-bas)/x(2:N))');

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy Variation element-wise via quadrature (Simpson's rule)
DE = (   dE(x(1:N),xhat(1:N),dx,dxhat,f(xhat(1:N)),0)...
      +4*dE(xm(1:N),xhatm(1:N),dx,dxhat,f(xhatm(1:N)),.5)...
        +dE(x(2:N+1),xhat(2:N+1),dx,dxhat,f(xhat(2:N+1)),1)   )/6;
    
%%  Fixed Boundary
DE = [0;DE;0];

clear xhatm dx dxhat;