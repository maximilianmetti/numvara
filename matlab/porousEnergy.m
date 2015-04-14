function E = porousEnergy(x,xhat,f,gamma)
%%  Evaluate the energy given by Power law for the FE flow map

if gamma==1,    E = fickEnergy(x,xhat,f);
else
    
H  = @(dx,dxhat,f) (f.^gamma)./(dx.^(gamma-1));

N  = length(x)-1;
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy element-wise via quadrature
E = dxhat'*(    H(dx,dxhat,f(xhat(1:N)))...
             +4*H(dx,dxhat,f(xhatm))...
             +  H(dx,dxhat,f(xhat(2:N+1))) )/6;
      
if min(dx)<=0, E=inf; end

clear xhatm dx dxhat;

end