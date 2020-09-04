function [x] = output_file(phi, theta, fn)

	fileID = fopen(strcat(fn, ".csv"), "w");
	% * \t Column 1 should contain the residue number indexed from 1\n
	% * \t Column 2 should contain the theta angle (in radians)\n
	% * \t Column 3 should contain the phi angle (in radians)\n
	for a=1:length(phi)
		fprintf(fileID, "%d %f %f\n", a, theta(a), phi(a));
	end
	fclose(fileID);
end

