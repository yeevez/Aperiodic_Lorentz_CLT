function s = same_side(p1,p2,a,b)
    cp1 = cross(b-a, p1-a);
    cp2 = cross(b-a, p2-a);
    s = dot(cp1, cp2) >= 0;
end