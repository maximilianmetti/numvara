function DE = generalFickVar(x,xhat,b,db,f)
%%  Evaluate the variation of Fick's energy for the FE flow map
N  = length(x)-1;
dE = @(basis,dx,dxhat,b,db,f) dxhat(1:N-1).*f(1:N-1)...
                                .*(-(1/dx(1:N-1))' + basis *db(1:N-1)./b(1:N-1))...
                            + dxhat(2:N).*f(2:N)...
                                .*( (1/dx(2:N))'   + (1-basis)*db(2:N)./b(2:N));

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy Variation element-wise via quadrature (Simpson's rule)
DE = (   dE(0.0,dx,dxhat,b(x(1:N)),  db(x(1:N)),  f(xhat(1:N) ))...
      +4*dE(0.5,dx,dxhat,b(xm(1:N)), db(xm(1:N)), f(xhatm(1:N)))...
        +dE(1.0,dx,dxhat,b(x(2:N+1)),db(x(2:N+1)),f(xhat(2:N+1)))  )/6;
    
%%  Fixed Boundary
DE = [0;DE;0];

clear xm xhatm dx dxhat;