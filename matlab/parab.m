function f = parab(x)

f = 4*(-x.*x+.25);
f( abs(x)>.5 ) = 0.0000001;

end

    
    