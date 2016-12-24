function v = uvw(pp,p1,p2)
    r0 = p2 - p1;
    r1 = pp - p1;
    r2 = pp - p2;
    r1xr2 = cross(r1,r2);
    v = 1/(4*pi) * (r0' * (r1/norm(r1) - r2/norm(r2))) *  r1xr2 / (r1xr2' * r1xr2);
end
