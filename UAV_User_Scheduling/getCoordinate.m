function [x, y] = getCoordinate(r)

	while(1) 
		x = (rand() * 2 - 1) * r * 0.6;
		y = (rand() * 2 - 1) * r * 0.6;
		
		if ((x^2 + y^2) <= (0.6 * r)^2) && ((x^2 + y^2) >= (0.3 * r)^2)
			break;
		end
	end
end
