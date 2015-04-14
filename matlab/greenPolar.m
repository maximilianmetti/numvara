function Q = greenPolar(r,s,delta,M)
%%  Computes Green interaction for Green'€™s function in Polar Coordinates
% M=100;delta = 1.e-14;

lr = length(r);
ls = length(s);
I  = zeros(M+1,M+1);
Q  = zeros(lr,ls);
ee = ones(M+1,1);
dphi = 2*pi/M;
phi  = 0:dphi:(2*pi);

for i=1:lr
for j=1:ls
    I =  (cos(phi)')*(cos(phi)) + (sin(phi)')*(sin(phi));
    I = -log( max(delta, r(i)*r(i)+s(j)*s(j) - 2*r(i)*s(j)*I ) );
    Q(i,j) = (0.25/pi)*dphi*dphi*ee'*I*ee;
end
end
	
Q = max(Q,0);

clear phi ee I;