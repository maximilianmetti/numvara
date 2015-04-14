function D = darcyDiss(x,xprev,xhat,dt,f)
%%  Evaluate the dissipation given by Darcy's law for the FE flow map
N  = length(x)-1;
D  = @(x,xp,dt,f) f.*(x-xp).^2/(dt*dt);

if length(x)-length(xprev)
    display('dimension mismatch in dissipation');
end
if length(x)-length(xhat)
    display('dimension mismatch in dissipation');
end
xm    = .5*(x(1:N)+x(2:N+1));
xpm   = .5*(xprev(1:N)+xprev(2:N+1));
xhatm = .5*(xhat(1:N)+xhat(2:N+1));
dxhat = -xhat(1:N)+xhat(2:N+1);

%%  Evaluate dissipation element-wise via quadrature
D = dxhat'*(   D(x(1:N),xprev(1:N),dt,f(xhat(1:N)))...
            +4*D(xm,xpm,dt,f(xhatm))...
            +  D(x(2:N+1),xprev(2:N+1),dt,f(xhat(2:N+1))) )/6;
             
clear xm xpm xhatm dxhat;