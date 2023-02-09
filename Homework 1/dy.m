function dy_z = dy(A, z)
    dy_z = (A/(2*C(z)^(3/2)))*((1-z*S(z))*dC(z)+2*(S(z)+z*dS(z))*C(z));
end