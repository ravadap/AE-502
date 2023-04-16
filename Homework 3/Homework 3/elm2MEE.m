function [h,k,p,q] = elm2MEE(e, i, w, OM)
    h = e.*sin(w+OM);
    k = e.*cos(w+OM);
    p = tan(i./2).*sin(OM);
    q = tan(i./2).*cos(OM);
end