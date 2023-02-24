function retrograde = is_retrograde(r1, v1)
    h = cross(r1,v1);
    if h(3) < 0
        retrograde = 1;
    else
        retrograde = 0;
    end
end