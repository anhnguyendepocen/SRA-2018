

v = zeros(size_cap,size_omega) + 100000;

for i = 1:size_cap ;
    for j=1:size_omega ;
        vexp = exp(omega_grid(1,j)+15)*2;
        vadd = sqrt(cap_grid(1,i)^(1.75)*sqrt(vexp));
        v(i,j) = v(i,j) + (vadd - v(i,j))/100000;
    end;
end;