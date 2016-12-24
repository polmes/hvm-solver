function v = uvwt(pp,p1,p2)
    r0 = p2 - p1;
    r1 = pp - p1;
    r2 = pp - p2;

    r1xr2 = cross(r1,r2);
    
    v = 1/(4*pi) * bsxfun(@times, dot(r0, bsxfun(@rdivide, r1, sqrt(sum(r1.^2))) - bsxfun(@rdivide, r2, sqrt(sum(r2.^2))), 1), bsxfun(@rdivide, r1xr2, dot(r1xr2, r1xr2, 1)));
end
