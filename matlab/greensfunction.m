f = @(y) 1;
G = @(x,y,M) 1 + 5*(0.5*(x-y)-M.*(x-y));

x = 0:.01:1;
y = 0:.01:1;
by= y(1:length(y)-1);
fy= y(2:length(y));
my= .5*(fy+by);
dy= fy-by;
M = zeros(3,length(my));
u = x;

for i=1:1:length(x)
    u(i) = 0.0;
    for j=1:1:length(y)-1
        if x(i) >= by(j)
            M(1,j) = 1;
        end
        if x(i) >= my(j)
            M(2,j) = 1;
        end
        if x(i) >= fy(j)
            M(3,j) = 1;
        end
    end
    xx   = x(i)*ones(1,length(my));
	GF = ( G(xx,by,M(1,:)).*f(by)...
       + 4*G(xx,my,M(2,:)).*f(my)...
       +   G(xx,fy,M(3,:)).*f(fy) )/6.;
    Qgf  = (GF)*dy';
    u(i) = u(i) + Qgf;
    Gr = ( G(xx,by,M(1,:))...
       + 4*G(xx,my,M(2,:))...
       +   G(xx,fy,M(3,:)) )/6.;
    %plot(my,Gr);axis([0 1 0 1]);pause(.1);
end

plot(x,u);