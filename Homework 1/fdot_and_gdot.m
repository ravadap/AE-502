function [fdot, gdot] = fdot_and_gdot(mu, alpha, r0, r, x)
    fdot = (sqrt(mu)/(r*r0))*(alpha*x^3*S(alpha*x^2)-x);
    gdot = 1 - (x^2/r)*C(alpha*x^2);
end