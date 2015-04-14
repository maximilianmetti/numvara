function [x,iter,DX] = nonlinconjgrad(x0,W,DW,tol,plotpath)
%% Conjugate Gradient using inexact line search with Armijo backtracking
% Ex1:
% F  = @(x) ((x(1)-2)*(x(2)-1)*(x(1)+3)*(x(2)+3))^2;
% DF = @(x) [ 2*(x(1)-2)*(x(2)-1)*(x(1)+3)*(x(2)+3)*(x(2)-1)*(x(2)+3)*(2*x(1)+1) ;
%             2*(x(1)-2)*(x(2)-1)*(x(1)+3)*(x(2)+3)*(x(1)-2)*(x(1)+3)*(2*x(2)+2)];
% Ex2:
% F  = @(x)  sin(x(1))+sin(1+x(2));
% DF = @(x) [cos(x(1));cos(1+x(2))];

%%  Initialize parameters
if isempty(tol), tol=10^(-10); end
if isempty(plotpath), plotpath=false; end
iter = 0;

%%  Initialize search direction (steepest descent)
x  = x0;
DX = -DW(x);
G  = W(x);
s  = DX;

%%  (Inexact) Line Search
a = 1;
armijo=false; %curve=false;
while ~armijo %|| ~curve
	a = a/5;                    % argmin_a (S/2+H)(x+a*s)
    if W(x+a*s)     <= G-.0001*a*s'*DX, armijo=true;
    else armijo=false;
    end
    
    %if s'*DF(x+a*s) >= .1*s'*DF(x), curve=true;
    %else curve=false;
    %end
    %display(a)
    %pause();
end

%X=[];Y=[];
while( norm(DX)>tol )
        %display(a)
        %display(norm(DX))
        DY= DX;           % store prev update direction
        x = x + a*s;      % update map
        
        %X = [X x(1)]; Y = [Y x(2)]; plot(X,Y);
        
        %% Compute Gradient and descent direction
        DX = -DW(x);
        G  = W(x);
        b  = max([0,DX'*(DX-DY)/(DY'*DY)]); % Polak-Ribiere with auto-reset
        s  = DX + b*s;
        if s'*DX<=0, s = DX; end
        iter = iter+1;
        
        if plotpath
        	plot(x0,x);pause(.1);
        end
        
        %%  (Inexact) Line Search
        sDX = s'*DX;
        a   = 1;
        wolfe=false; armijo=false; curve=false;
        while ~wolfe
            
            %   Armijo Condition
            if W(x+a*s) <= G+0.0001*a*sDX, armijo=true; end
            while ~armijo
                a = .9*a;
                if W(x+a*s) <= G+0.0001*a*sDX, armijo=true; end
            end
            
            %   Curvature Condition
            if s'*DW(x+a*s) >= .9*sDX, curve=true; end
            while ~curve
                a = 1.1*a;
                if s'*DW(x+a*s) >= .9*sDX, curve=true; end
            end;
            
            if W(x+a*s) > G+0.0001*a*sDX, armijo=false; end
            if curve && armijo, wolfe=true; end
        end
end

x = x + a*s;      % update map
clear s;
DX = DW(x);
