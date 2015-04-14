function E = generalFickEnergy(x,xhat,b,f)
%%  Evaluate the energy given by Fisher-Wrigt Energy for the FE flow map
N  = length(x)-1;
H  = @(dx,dxhat,b,f) f.*log(b.*f./(dx./dxhat));

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy element-wise via quadrature
E = dxhat'*(    H(dx,dxhat,b(x(1:N)),f(xhat(1:N)))...
             +4*H(dx,dxhat,b(xm),f(xhatm))...
             +  H(dx,dxhat,b(x(2:N+1)),f(xhat(2:N+1))) )/6;
      
if min(dx)<=0, E=inf; end

clear xm xhatm dx dxhat;