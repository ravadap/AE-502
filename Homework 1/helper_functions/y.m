function y_z = y(r1, r2, A, z)
    y_z = r1 + r2 + A*((z*S(z)-1)/sqrt(C(z)));
end