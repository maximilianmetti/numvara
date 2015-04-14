function M = moment(m,x,xhat,f)
%%  Evaluate the energy given by Fick's law for the FE flow map
N  = length(x)-1;
H  = @(m,x,f) (x.^m).*f;

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate Energy element-wise via quadrature
M = dxhat'*(    H(m,x(1:N),f(xhat(1:N)))...
             +4*H(m,xm,f(xhatm))...
             +  H(m,x(2:N+1),f(xhat(2:N+1))) )/6;
      
if min(dx)<=0, M=inf; end

clear xm xhatm dx dxhat;