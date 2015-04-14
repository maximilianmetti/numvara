function [x,xprev,xhat,iter,DX] = Newton(x0,x0prev,x0hat,W,DW,DDW,tol,maxit,plotpath)
%% Conjugate Gradient using inexact line search with Armijo backtracking
% Ex1:
% F  = @(x) ((x(1)-2)*(x(2)-1)*(x(1)+3)*(x(2)+3))^2;
% DF = @(x) [ 2*(x(1)-2)*(x(2)-1)*(x(1)+3)*(x(2)+3)*(x(2)-1)*(x(2)+3)*(2*x(1)+1) ;
%             2*(x(1)-2)*(x(2)-1)*(x(1)+3)*(x(2)+3)*(x(1)-2)*(x(1)+3)*(2*x(2)+2)];
% Ex2:
% F  = @(x)  sin(x(1))+sin(1+x(2));
% DF = @(x) [0;cos(x(1));cos(1+x(2));0];
% x0 = [0;1;1;0];
% plotpath=false;
% tol= .0000001

%%  Initialize parameters
if isempty(tol), tol=10^(-5); end
if isempty(plotpath), plotpath=false; end
iter = 0;
N    = length(x0)-1;

%%  Initialize search direction (steepest descent)
x  = x0;
xprev = x0prev;
xhat  = x0hat;
dx = x(2:N+1)-x(1:N);
quasiFactor = 0.0;
hmin = quasiFactor*mean(dx);
df = DW(x,xprev,xhat);
r0 = norm(df);
B  = DDW(x,xprev,xhat);
        
while( max(df) > tol && iter<maxit)        
        %%  Solve Hessian System
        s = zeros(length(x),1);
        s(2:N) = -B\df(2:N);

        %%  (Inexact) Line Search
        G      = W(x,xprev,xhat);
        sDF    = s'*df;
        if sDF>0, s=-s; sDF=-sDF; end
        a = 1; k=0;
        wolfe=false; armijo=false; curve=false;
        while ~wolfe && k<100
            %   Armijo Condition
            if W(x+a*s,xprev,xhat) <= G+0.0001*a*sDF, armijo=true; end
            
            j=0;
            while ~armijo && j<100
                a = .9*a;
                j = j+1;
                if W(x+a*s,xprev,xhat) <= G+0.0001*a*sDF, armijo=true; end
            end
            
            %   Curvature Condition
            if s'*DW(x+a*s,xprev,xhat) >= .9*sDF, curve=true; end
            j=0;
            while ~curve && j<100
                a = 1.1*a;
                j = j+1;
                if s'*DW(x+a*s,xprev,xhat) >= .9*sDF, curve=true; end
            end;
            
            if W(x+a*s,xprev,xhat) > G+0.0001*a*sDF, armijo=false; end
            if curve && armijo, wolfe=true; end
            k = k+1;
        end
        % Ensure Decrease
        j=0;
        while ~armijo && j<20
            a = .9*a;
            j = j+1;
            if W(x+a*s,xprev,xhat) <= G+0.0001*a*sDF, armijo=true; end
        end
        %if k==100, s = -df; display('Wolfe failed'); end
        x = x + a*s;      % update map
        
        %%  Adaptivity
        toosmall = true;
        while toosmall
             dx = x(3:end)-x(1:end-2);
             good = find(dx>=2*hmin);
            if length(good)==length(dx), toosmall=false;
            else
                x     = x([1;(good)+1;end]);
                xprev = xprev([1;(good)+1;end]);
                xhat  = xhat([1;(good)+1;end]);
            end
        end
        k=2;
        while k<length(x)
            if x(k+1)-x(k-1) < 2*hmin
                x     = [x(1:k-1);x((k+1):end)];
                xprev = [xprev(1:k-1);xprev((k+1):end)];
                xhat  = [xhat(1:k-1);xhat((k+1):end)];
                iter  = 0;
            else k =k+1;
            end
        end
        N  = length(x)-1;
        hmin = quasiFactor*mean(dx);
        
        %X = [X x(1)]; Y = [Y x(2)]; plot(X,Y);
        
        %% Compute Gradient and descent direction with BFGS Hessian approx
        %dfp= df;
        %display(norm(df))
        df = DW(x,xprev,xhat);
        %y  = df-dfp;
        B  = DDW(x,xprev,xhat);%B + y*y'/(a*s'*y) - B*(s*s')*B'/(s'*B*s);
        
        
        iter = iter+1;
        if iter==90, display('Still searching with Newton'); end
        
        if plotpath
        	plot(xhat,xhat,xhat,x);
            title(['At iteration ', num2str(iter), ' there are ',num2str(N), ' nodes']);
            pause(.1);
        end
end
if iter==160, display('Newton Failed'); end
clear df dfp y B s;

DX = DW(x,xprev,xhat);
