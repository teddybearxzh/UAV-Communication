function [DAE, MT, d] = AE_MT_generating(cp_total, number, M, r, n)
    original_point = cp_total{number};
    DAE = cell(1, M);
    
    if M == 3
        random_DAElength = unifrnd(0.4, 0.7);
        DAE{1} = original_point + [0, random_DAElength * r];
        random_DAElength = unifrnd(0.4, 0.7);
        DAE{2} = original_point + [-random_DAElength * r * sin(pi/3), -random_DAElength * r * sin(pi/6)];
        random_DAElength = unifrnd(0.4, 0.7);
        DAE{3} = original_point + [random_DAElength * r * sin(pi/3), -random_DAElength * r * sin(pi/6)];
        
    elseif M == 4
        random_DAElength = unifrnd(0.4, 0.7);
        random_DAElength = 0.6;
        DAE{1} = original_point + [random_DAElength * r * sin(pi/4), random_DAElength * r * sin(pi/4)];
        
        random_DAElength = unifrnd(0.4, 0.7);
        random_DAElength = 0.6;
        DAE{2} = original_point + [-random_DAElength * r * sin(pi/4), random_DAElength * r * sin(pi/4)];
        
        random_DAElength = unifrnd(0.4, 0.7);
        random_DAElength = 0.6;
        DAE{3} = original_point + [-random_DAElength * r * sin(pi/4), -random_DAElength * r * sin(pi/4)];
        
        random_DAElength = unifrnd(0.4, 0.7);
        random_DAElength = 0.6;
        DAE{4} = original_point + [random_DAElength * r * sin(pi/4), -random_DAElength * r * sin(pi/4)];
        d = [0.95,0.78,0.81,0.86,0.90,0.93,0.95,0.97,0.99,0.998];
    end
    
    while(1)
        [MT_x, MT_y] = getCoordinate(r);
        MT = [MT_x, MT_y] + original_point;
        flag = 1;
        for i = 1 : M
            if norm(DAE{i} - MT) <= r/n
                flag = 0;
            end
        end
        if flag == 1
            break
        end
    end
end
