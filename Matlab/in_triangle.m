function b = in_triangle(p,p1,p2,p3)
    b = same_side(p,p1,p2,p3) && same_side(p,p2,p1,p3) && same_side(p,p3,p1,p2);
end