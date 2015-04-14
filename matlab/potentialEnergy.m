function E = potentialEnergy(x,xhat,f,phi)
%%  Evaluate the energy given by Potential energy for the FE flow map
N  = length(x)-1;
H  = @(f,phi) f.*phi;

xm    = .5*(x(1:N)+x(2:N+1));
dx    = -x(1:N)+x(2:N+1);
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy element-wise via quadrature
E = dxhat'*(   H(f(xhat(1:N)),phi(x(1:N)))...
            +4*H(f(xhatm),phi(xm))...
            +  H(f(xhat(2:N+1)),phi(x(2:N+1))) )/6;
      
if min(dx)<=0, E=inf; end

clear xhatm xm dx dxhat;