%%
% Rs(入射角, 侵入側屈折率,　出射側屈折率)
%%
function reflectivity = Rs(a, n_a, n_b)
%radians to degrees
    if a == 0
        a = 0.1;
    end
    a = deg2rad(a);
    b = asin(sin(a) * (n_a / n_b));
    
    reflectivity = ((sin(b-a)/sin(b+a)) ^2) * 100; % 論文(3.1)
end