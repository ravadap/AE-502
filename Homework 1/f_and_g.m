function [f,g] = f_and_g(mu, alpha, r0, dt, x)
    f = 1 - (x^2/r0)*C(alpha*x^2);
    g = dt - (1/sqrt(mu))*x^3*S(alpha*x^2);
end