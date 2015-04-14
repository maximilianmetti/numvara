function E = radialFickEnergy(x,xhat,f)
%%  Evaluate the energy given by Fick's law for the FE flow map
N  = length(x)-1;
H  = @(x,xhat,dx,dxhat,f,basis) xhat.*f.*log(xhat.*f./(x.*dx./dxhat));

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy element-wise via quadrature
E = dxhat'*(    H(x(1:N),xhat(1:N),dx,dxhat,f(xhat(1:N)))...
             +4*H(xm,xhatm,dx,dxhat,f(xhatm))...
             +  H(x(2:N+1),xhat(2:N+1),dx,dxhat,f(xhat(2:N+1))) )/6;
      
if min(dx)<=0, E=inf; end

clear xhatm dx dxhat;