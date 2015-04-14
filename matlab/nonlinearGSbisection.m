function [x,xprev,xhat] = nonlinearGSbisection(x0,x0prev,x0hat,W,DW,nonlintol,maxsweep,maxbisect)

%% Measure Variation and initialize
fn   = W(x0,x0prev,x0hat);
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
        %if plotit, display(['forsweep ', num2str(i)]); end
        a = x(i-1); b = x(i+1); y = x(i);
        dy= .5*min(y-a,b-y);
        %% Bisction method
        for j=1:maxbisect
            c = y-dy; d = y+dy;
            nod = [ a c y d b ];
            xc(i) = c; xd(i) = d;
            fc = W(xc,xprev,xhat);
            fd = W(xd,xprev,xhat);
            F  = [ fc fn fd ];
            fm = min(F);
            idx= find(F==fm)+1;
            if length(idx)>1, fm = F(2); idx = 3; end %tie-break
            a = nod(idx-1); b = nod(idx+1);
            y = nod(idx); x(i) = y; fn = fm;
            dy= .5*min(y-a,b-y);
            if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
        end
        xc(i) = y; xd(i) = y;
        if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
    end
    
    %% Backward sweep
    for k=2:N-1
        i=N+1-k;
        %if plotit, display(['backsweep ', num2str(k)]); end
        a = x(i-1); b = x(i+1); y = x(i);
        dy= .5*min(y-a,b-y);
        %% Bisction method
        for j=1:maxbisect
            c = y-dy; d = y+dy;
            nod = [ a c y d b ];
            xc(i) = c; xd(i) = d;
            fc = W(xc,xprev,xhat);
            fd = W(xd,xprev,xhat);
            F  = [ fc fn fd ];
            fm = min(F);
            idx= find(F==fm)+1;
            if length(idx)>1, fm = F(2); idx = 3; end %tie-break
            a = nod(idx-1); b = nod(idx+1);
            y = nod(idx); x(i) = y; fn = fm;
            dy= .5*min(y-a,b-y);
            if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
        end
        xc(i) = y; xd(i) = y;
        if plotit, plot(xhat,xhat,xhat,x); pause(.01); end
    end
    
    %% Measure Variation
    obj = norm(DW(x,xprev,xhat));
    display(['Relative Residual after ', num2str(sweeps), ' sweeps: ', num2str(obj/obj0)]);
    if true, plot(xhat,xhat,xhat,x); pause(.01); end
end
        

