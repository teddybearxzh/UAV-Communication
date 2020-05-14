function [x, y] = getCoordinate(r)

	while(1) 
		x = (rand() * 2 - 1) * r * 0.5;
		y = (rand() * 2 - 1) * r * 0.5;
		
		if ((x^2 + y^2) <= (0.5 * r)^2) && ((x^2 + y^2) >= (0.2 * r)^2)
			break;
		end
	end
end
