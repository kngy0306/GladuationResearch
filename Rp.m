function reflectivity = Rp(a, n_a, n_b)
%radians to degrees
    if a == 0
        a = 0.1;
    end
    a = deg2rad(a);
    b = asin(sin(a) * (n_a / n_b));
    
    reflectivity = ((tan(b-a)/tan(b+a)) ^2) * 100; % 論文(3.2)
end