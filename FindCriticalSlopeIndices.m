%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Given the x-axes and the wave height, this function returns the range of indices
%	that exceed the critical slope value of 1/7.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mins, maxes] = FindCriticalSlopeIndices( x, wave )
	deltaX = x(2)-x(1);
	slope = diff(wave)/deltaX;
	criticalIndices = find( abs(slope) > 1/7 );


	mins = [];
	maxes = [];
	
	% This algorithm searches to see if this point to neighboring points has a slope
	% that exceeds 1/7. For example, lets say the slope between indices 4 and 5 exceeds
	% the threshold. This looks to see if the slope between 3 and 5 also exceeds the
	% threshold. This looks both backwards and forwards.
	for index=criticalIndices'
		if ( ismember( index, maxes) )
			continue;
		end

		minIndex = index;
		maxIndex = index+1;
	
		mySlope = abs((wave(maxIndex) - wave(minIndex))/(x(maxIndex) - x(minIndex)));	
	
		while ( mySlope > 1/7 )
			minIndex = minIndex-1;
			mySlope = abs((wave(maxIndex) - wave(minIndex))/(x(maxIndex) - x(minIndex)));		
		end
		minIndex = minIndex+1;
	
		mySlope = abs((wave(maxIndex) - wave(index))/(x(maxIndex) - x(index)));
	
		while ( mySlope > 1/7 )
			maxIndex = maxIndex+1;
			mySlope = abs((wave(maxIndex) - wave(index))/(x(maxIndex) - x(index)));	
		end
		maxIndex = maxIndex-1;
	
		mins = [mins; minIndex];
		maxes = [maxes; maxIndex];	
	end
	
	% this part of the algorithm groups overlapping index ranges together.
	i=1;
	while ( i < length(maxes)-1 )
		if ( maxes(i) >= mins(i+1) )
			maxes(i) = maxes(i+1);
			mins(i+1)=[];
			maxes(i+1)=[];
		else
			i=i+1;
		end
	end
	
	% This next part of the algorithm groups together ranges that are part of the same
	% wave.
	
% 	i=1;
% 	while ( i < length(maxes)-1 )
% 		if ( maxes(i) >= mins(i+1) )
% 			maxes(i) = maxes(i+1);
% 			mins(i+1)=[];
% 			maxes(i+1)=[];
% 		else
% 			i=i+1;
% 		end
% 	end
end