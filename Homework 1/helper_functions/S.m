function S_z = S(z)
    if z > 0
        S_z = (sqrt(z)-sin(sqrt(z)))/(sqrt(z)^3);
    elseif z < 0
        S_z = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z)^3);
    else
        S_z = 1/6;
    end
end