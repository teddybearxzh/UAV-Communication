function [UAV_position, user] = userUAVgenerating(centers, R, r1, r2, M, K, n_dis, number)

    original_position = centers{number};
    
    UAV_position = cell(1, M);
    for i=1:M
        %theta = 2*pi/M*(i-1) + (2*pi/M*i-2*pi/M*(i-1)) * rand();
        theta = pi/2 + (i-1)*2*pi/3;
        r = r1 + (r2-r1) * rand();
        UAV_x = r * cos(theta);
        UAV_y = r * sin(theta);
        UAV_position{i} = [UAV_x, UAV_y] + original_position;
    end
    
    user = cell(1, K);
	for i=1:K
        if i==1
            [user_x, user_y] = getCoordinate(R);
            user{i} = [user_x, user_y] + original_position;
        elseif i>1
            while(1)
                [user_x, user_y] = getCoordinate(R);
                flag = 1;
                for j=1:(i-1)
                    if norm([user_x, user_y] + original_position - user{j}) < (R/n_dis)
                        flag = 0;
                        break;
                    else
                        continue;
                    end
                end
                if flag==1
                    user{i} = [user_x, user_y] + original_position;
                    break;
                else
                    continue;
                end
            end
        end
    end
end


