% https://ask.metafilter.com/49252/rotation-shmotation flagrantly copied from here.
function newr = rotateIt(vectorToBeRotatedAbout, angle, vectorToBeRotated)
	r = vectorToBeRotated;
	n = vectorToBeRotatedAbout;
	newr = r*cos(angle) + n*dot(n,r)*(1-cos(angle)) + cross(r,n)*sin(angle);
end
