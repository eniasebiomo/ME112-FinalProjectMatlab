% Simple function to plot lines between a series
% of points. At the end, we wrap around to beginning
% to make a closed polygon.
% 3Mar2010; notes added 30Dec2010 -mrc

function plotlines(points,color1,color2)
%Points should be array with X values in 1st row, Y in 2nd row.
%Any other rows (e.g. ones for homogeneous coords) are ignored
%color1, color2 = 'r' or 'g' etc.
%Possible bug: using length(points) only works if you 
%have at least as many points as rows. So you need
%at least 3 points for homogeneous coords (x,y,1).

for i=1:length(points)-1
 X1 = points(1,i);
 Y1 = points(2,i);
 X2 = points(1,i+1);
 Y2 = points(2,i+1);
 line([X1;X2],[Y1;Y2],'Color',color1,'LineWidth',1); 
end
%Last one where we wrap around to beginning again...
 X1 = points(1,length(points));
 Y1 = points(2,length(points));
 X2 = points(1,1);
 Y2 = points(2,1);
 line([X1;X2],[Y1;Y2],'Color',color2,'LineWidth',1); 

end