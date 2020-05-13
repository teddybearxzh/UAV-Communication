function [x, y] = getCoordinate(r)
%在以原点为圆心，r为半径的圆内，随机产生点
%要求：点不能离圆心太远，也不能太近

	while(1) %不断产生x和y
		x = (rand() * 2 - 1) * r * 0.6;
		y = (rand() * 2 - 1) * r * 0.6;
        %直到x，y满足一定条件
		if ((x^2 + y^2) <= (0.6 * r)^2) && ((x^2 + y^2) >= (0.3 * r)^2)
			break;
		end
	end
end