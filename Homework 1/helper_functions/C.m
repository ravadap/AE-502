function C_z = C(z)
    if z > 0
        C_z = (1-cos(sqrt(z)))/z;
    elseif z < 0
        C_z = (cosh(sqrt(-z))-1)/(-z);
    else
        C_z = 1/2;
    end
end