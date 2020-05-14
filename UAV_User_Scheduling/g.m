function G = g(total_cell, m, k, n, l, height, M)

%parameters：
%total_cell：cell information
%m，n：UAV index
%k：UAV internal user index
%l：hovering time index
%height：UAV height
%M：UAV location number

    %initial from[12]
    LOS = 0.1;
    NLOS = 21;
    a = 5.0188;
    b = 0.3511;
    f = 2000000000;
    c = 300000000;
    
    %plane distance
    if mod(l,M)==0
        d_2d = norm(total_cell{m}{2}{M} - total_cell{n}{3}{k});
    else
        d_2d = norm(total_cell{m}{2}{mod(l,M)} - total_cell{n}{3}{k});
    end
    
    %space distance
    d = sqrt(d_2d^2 + height^2);
    
    %calculation
    A = LOS - NLOS;
    B = 20*log10(d) + 20*log10(4*pi*f/c) + NLOS;
    theta = 180/pi * asin(height/d);
    
    G = A / (1+a*exp(-b*(theta-a))) + B;
    G = power(10, -G/10); 
end

