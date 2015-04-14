function [x,xprev,xhat] = nonlinearGSnewton(x0,x0prev,x0hat,W,DW,DDW,nonlintol,maxsweep)

%% Measure Variation and initialize
%fn   = W(x0,x0prev,x0hat);
obj0 = norm(DW(x0,x0prev,x0hat));
obj  = obj0;
plotit=true;
%plotit=false;

N  = length(x0);
x  = x0; xprev = x0prev; xhat = x0hat;
xc = x0; xd = x0;
sweeps = 0;

%% Nonlinear Gauss-Seidel
while (obj/obj0 > nonlintol) && (sweeps<maxsweep)
    sweeps = sweeps+1;
    
    %% Forward sweep
    for i=2:N-1
        %% Newton method
        df  = DW([x(i-1);x(i);x(i+1)],[xprev(i-1);xprev(i);xprev(i+1)],[xhat(i-1);xhat(i);xhat(i+1)]);
        ddf = DDW([0;x(i-1);x(i);x(i+1);0],[0;xprev(i-1);xprev(i);xprev(i+1);0],[0;xhat(i-1);xhat(i);xhat(i+1);0]);
        if ddf(2,2)~=0, dy  = -df(2)/ddf(2,2);
        else dy = -df(2)*(x(i+1)-x(i-1))/abs(df(2)); end
        if x(i)+dy < x(i-1)+1e-16
            x(i) = x(i) - (x(i)-x(i-1))*.9;%tanh(dy/(x(i-1)-x(i)));
        elseif x(i)+dy > x(i+1)-1e-16
            x(i) = x(i) + (x(i+1)-x(i))*.9;%tanh(dy/(x(i+1)-x(i)));
        else
            x(i) = x(i)+dy;
        end
        if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
    end
    
    %% Backward sweep
    for k=2:N-1
        i=N+1-k;
        %% Newton method
        df  = DW([x(i-1);x(i);x(i+1)],[xprev(i-1);xprev(i);xprev(i+1)],[xhat(i-1);xhat(i);xhat(i+1)]);
        ddf = DDW([0;x(i-1);x(i);x(i+1);0],[0;xprev(i-1);xprev(i);xprev(i+1);0],[0;xhat(i-1);xhat(i);xhat(i+1);0]);
        if ddf(2,2)~=0, dy  = -df(2)/ddf(2,2);
        else dy = -df(2)*(x(i+1)-x(i-1))/abs(df(2)); end
        if x(i)+dy < x(i-1)+1e-16
            x(i) = x(i) - (x(i)-x(i-1))*.5;%*tanh(dy/(x(i-1)-x(i)));
        elseif x(i)+dy > x(i+1)-1e-16
            x(i) = x(i) + (x(i+1)-x(i))*.5;%tanh(dy/(x(i+1)-x(i)));
        else
            x(i) = x(i)+dy;
        end
        if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
    end
    
    %% Measure Variation
    obj = norm(DW(x,xprev,xhat));
    display(['Relative Residual after ', num2str(sweeps), ' sweeps: ', num2str(obj/obj0)]);
    if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
end
        

