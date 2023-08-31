function [boolean, posOfColl, newV] = collWithHyperCube(pos, center, vel, side)
	posOfColl = [];
    newV = [];
    %make it able to handle lists for centerCube????
    boolean = false;
	%vector velocity signs
	vSign = sign(vel);
	%size of dimensions of velocity
	ldim = length(vel);
	%center of cube relative to position
	nwc = center - pos;
	%possible collision side 
	pcsd = nwc - (side/2)*vSign;
	%possible collision side sign
	pcss = sign(pcsd);
	for i = 1:ldim
		if ( vSign(i) == pcss(i) )
			%possible collision point
			pcp = vel*pcsd(i)/vel(i);
			%if the biggest difference between the center of the side and the projection is bigger on some margin of the hypercube then it does not hit
			if( max( abs(pcp - nwc) - 0.000001 ) <= side/2 )
				boolean = true;
				%It only bounces on one side on which it bounces in the same direction
				newV = vel;
				newV(i) = newV(i)*(-1);
				posOfColl = pcp + pos;
				return
			end
		end
	end
end
