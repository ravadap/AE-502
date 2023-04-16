function [r,v] = elm2rv_PR(a,e,i,w,OM,f,mu)
    % compute angular momentum
    h = sqrt(a*mu*(1-e^2));

    % r and v in perifocal frame
    r_xyz = (h^2/mu)*(1/(1+e*cos(f)))*[cos(f) sin(f) 0]';
    v_xyz = (mu/h)*[-sin(f) e+cos(f) 0]';

    % transformation from perifocal to geocentric
    Q_xX = [cos(OM)*cos(w)-sin(OM)*sin(w)*cos(i) -cos(OM)*sin(w)-sin(OM)*cos(i)*cos(w) sin(OM)*sin(i);
            sin(OM)*cos(w)+cos(OM)*cos(i)*sin(w) -sin(OM)*sin(w)+cos(OM)*cos(i)*cos(w) -cos(OM)*sin(i);
            sin(i)*sin(w) sin(i)*cos(w) cos(i)];

    % transform to geocentric frame
    r = Q_xX*r_xyz;
    v = Q_xX*v_xyz;
end