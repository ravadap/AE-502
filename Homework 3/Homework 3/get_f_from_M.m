function f = get_f_from_M(M, e, tol)
    E = M;
    ratio = 1;

    while abs(ratio) > tol
        f = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        ratio = f/fp;
        E = E - ratio;
    end

    f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end