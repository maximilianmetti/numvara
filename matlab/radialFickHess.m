function DDE = radialFickHess(x,xhat,f)
%%  Evaluate the Hessian of Fick's energy for the FE flow map
N  = length(x)-1;
ddE = @(x,xhat,dx,dxhat,f,bas) gallery('tridiag',...
            -dxhat(2:N-1).*xhat(2:N-1).*f(2:N-1).*((1/dx(2:N-1).^2)'-(bas*(1-bas)/x(2:N-1).^2)'),...
             dxhat(1:N-1).*xhat(1:N-1).*f(1:N-1).*((1/dx(1:N-1).^2)'+(bas*bas/x(1:N-1).^2)')...
            +  dxhat(2:N).*xhat(2:N).*f(2:N).*((1/dx(2:N).^2)'+((1-bas)*(1-bas)/x(2:N).^2)'),...
            -dxhat(2:N-1).*xhat(2:N-1).*f(2:N-1).*((1/dx(2:N-1).^2)'-(bas*(1-bas)/x(2:N-1).^2)') );

xm    = .5*(x(1:N)+x(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dx    = -x(1:N)+x(2:N+1);
dxhat = -xhat(1:N)+xhat(2:N+1);


%%  Evaluate Energy Hessian element-wise via quadrature (Simpson's rule)
DDE = (  ddE(x(1:N),xhat(1:N),dx,dxhat,f(xhat(1:N)),0)...
      +4*ddE(xm,xhatm,dx,dxhat,f(xhatm),.5)...
      	+ddE(x(2:N+1),xhat(2:N+1),dx,dxhat,f(xhat(2:N+1)),1)   )/6;
    
clear xhatm dx dxhat;