function [ p,ph,pl ] = Adjust_Price( D,S,p,ph,pl,tol )

if abs(D-S)>tol
    if D>S
        pl=p;
        p=0.5*p+0.5*ph;
    else
        ph=p;
        p=0.5*p+0.5*pl;
    end
end
        



end

