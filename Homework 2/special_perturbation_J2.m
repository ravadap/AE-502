function dydt = special_perturbation_J2(~, y, mu, R, J2)
    % position mag
    r = norm([y(1) y(2) y(3)]);

    % compute J2 pertubation
    I = eye(3);
    J2_perturbation = (3*J2*mu*R^2/(2*r^5))* ((5*y(3)^2/r^2-1) * (y(1)*I(1,:)+y(2)*I(2,:)) + y(3)*(5*y(3)^2/r^2-3)*I(3,:));
    
    dydt = [y(4) y(5) y(6) -(mu/r^3)*y(1)+J2_perturbation(1) -(mu/r^3)*y(2)+J2_perturbation(2) -(mu/r^3)*y(3)+J2_perturbation(3)]';
end