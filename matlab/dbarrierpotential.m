function dphi = dbarrierpotential(x,w,h,z)
% Returns a barrier function that is zero for x<z
% The barrier reaches a height of h with a continuous
% linear interface of width w
    dphi          = zeros(length(x),1);
    intface       = find(x>z & x<z+w);
    dphi(intface) = h*ones(length(intface),1)/w;
end