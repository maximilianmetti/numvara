function [x,iter,DX] = quasiNewton(x0,W,DW,tol,plotpath)
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
if isempty(tol), tol=10^(-10); end
if isempty(plotpath), plotpath=false; end
iter = 0;
N    = length(x0)-1;

%%  Initialize search direction (steepest descent)
x  = x0;
B  = eye(length(x));
df = DW(x);
s  = zeros(N+1,1);
ds = zeros(N+1,1);
        
%X=[];Y=[];
while( norm(df)>tol )        
        %%  Solve Hessian System by Gauss-Seidel
        s(2:N) = -B(2:N,2:N)\df(2:N);
        
%         L = tril(B(2:N,2:N));
%         r  = -df;%-B*s;
%         k = 0;
%         while norm(r)>.001*norm(df) && k<500
%             ds(2:N) = L\r(2:N);
%             s  = s+ds;
%             r  = -df-B*ds;
%             k  = k+1;
%         end
%         if k==500, s = -df; end
%         if sum(isnan(s))>0, s = -df; end

        %%  (Inexact) Line Search
        G      = W(x);
        sDF    = s'*df;
        if sDF>0, s=-s; sDF=-sDF; end
        a = 1; %k=0;
        wolfe=false; armijo=false; curve=false;
        while ~wolfe %&& k<100
            %   Armijo Condition
            if W(x+a*s) <= G+0.0001*a*sDF, armijo=true; end
            
            %j=0;
            while ~armijo %&& j<100
                a = .9*a;
                %j = j+1;
                if W(x+a*s) <= G+0.0001*a*sDF, armijo=true; end
            end
            
            %   Curvature Condition
            if s'*DW(x+a*s) >= .9*sDF, curve=true; end
            %j=0;
            while ~curve %&& j<100
                a = 1.1*a;
                %j = j+1;
                if s'*DW(x+a*s) >= .9*sDF, curve=true; end
            end;
            
            if W(x+a*s) > G+0.0001*a*sDF, armijo=false; end
            if curve && armijo, wolfe=true; end
            %k = k+1;
        end
        % Ensure Decrease
        %j=0;
        while ~armijo %&& j<20
            a = .9*a;
            %j = j+1;
            if W(x+a*s) <= G+0.0001*a*sDF, armijo=true; end
        end
        %if k==100, s = -df; display('Wolfe failed'); end
        x = x + a*s;      % update map
        
        %X = [X x(1)]; Y = [Y x(2)]; plot(X,Y);
        
        %% Compute Gradient and descent direction with BFGS Hessian approx
        dfp= df;
        %display(norm(df))
        df = DW(x);
        y  = df-dfp;
        B  = B + y*y'/(a*s'*y) - B*(s*s')*B'/(s'*B*s);
        iter = iter+1;
        
        if plotpath
        	plot(x0,x);pause(.1);
        end
end
%if iter==100, display('Quasi-Newton Failed'); end
clear df dfp y B s;

DX = DW(x);
