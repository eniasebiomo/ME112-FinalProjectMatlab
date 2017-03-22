% Function to create 3x3 affine transformation 
% corresponding do a simple rotation of an object
% by angle theta anticlockwise about the origin 
% in the plane:  p2 = R*p1 where 
% p1 = [px, py 1]' (using homogeneous coordinates)
% and R is the translation operator.
% usage: R = Rotq(theta)
% 2Mar2010 -mrc

function Trans = Rotq (Theta)

ct = cos(Theta);
st = sin(Theta);

Trans = [ ct -st  0.0;
          st  ct  0.0;
          0.0 0.0 1.0 ];

end
