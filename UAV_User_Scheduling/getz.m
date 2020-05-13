function z = getz(zz, order)
    for i = 1 : length(order)
        zz(i, order(i)) = 1;
        z = zz;
    end
end