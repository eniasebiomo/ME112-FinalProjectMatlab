% Function to create 3x3 affine transformation 
% corresponding do a simple X,Y translation of
% an object in the plane: p2 = T*p1 where 
% p1 = [px, py 1]' (using homogeneous coordinates)
% and T is the translation operator.
% usage: T = Tranxy(1,1)
% 2Mar2010 -mrc

function Trans = Tranxy (X, Y)

Trans = [1.0  0.0  X;
         0.0  1.0  Y;
         0.0  0.0  1.0];

end
