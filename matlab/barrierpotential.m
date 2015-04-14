function phi = barrierpotential(x,w,h,z)
% Returns a barrier function that is zero for x<z
% The barrier reaches a height of h with a continuous
% linear interface of width w
    phi         = zeros(length(x),1);
    phi(x>=z)   = h*(x(x>=z)-z)/w;
    phi(x>=z+w) = h*ones(length(find(x>=z+w)),1);
end