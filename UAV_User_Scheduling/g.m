function G = g(total_cell, m, k, n, l, height, M)

%参数说明：
%total_cell：cell结构，包含UAV的cell中心点，UAV位置分布，内部用户分布
%m，n：UAV序号
%k：UAV内部用户序号
%l：hovering time序号
%height：UAV飞行高度
%M：UAV位置个数

    %参数初始化 from[12]
    LOS = 0.1;
    NLOS = 21;
    a = 5.0188;
    b = 0.3511;
    f = 2000000000;
    c = 300000000;
    
    %平面距离
    if mod(l,M)==0
        d_2d = norm(total_cell{m}{2}{M} - total_cell{n}{3}{k});
    else
        d_2d = norm(total_cell{m}{2}{mod(l,M)} - total_cell{n}{3}{k});
    end
    %空间距离
    d = sqrt(d_2d^2 + height^2);
    
    %具体计算
    A = LOS - NLOS;
    B = 20*log10(d) + 20*log10(4*pi*f/c) + NLOS;
    theta = 180/pi * asin(height/d);
    
    G = A / (1+a*exp(-b*(theta-a))) + B;
    G = power(10, -G/10); 
end

