% Updated slightly - Feb19 2011 -mrc
% Tutorial file showing how the transformations work
% for a collection of points. 
  
% Assume data points are initially in (x,y) pairs for each row
data = [1.0, 1.0; 0.7, 1.5; 1.0, 2.0; 2.25, 2.0; 2.25, 1.0]

% Make array where we augment data with "1"
% to get homogenous coordinates:
rcol = ones(length(data),1);
augmented = cat(2,data,rcol);
% Transform it to a series of columns:
oldpoints = augmented'

%Define our translation and rotation
%(Opposite if moving coord frame by same amount)
T1 = Tranxy(1.0,1.0);     % Translate by (x,y)
R1 = Rotq(pi/4);          % Rotate about origin by theta

% Check that inverse is same as transforming
% by negative amount:
NegT1 = Tranxy(-1.0,-1.0)
disp('matrix inverse of T1:');
inverse(T1)
NegR1 = Rotq(-pi/4)
disp('matrix inverse of R1:')
inverse(R1)

%%%%%%%%%%%%%%%
%Set up axes for plotting
clf;
figure(1);
%Guess a reasonable plot area (may need to tweak these)
leftlim = -1;
rightlim = 4;
toplim = 4;
bottomlim = -1;
axis([leftlim rightlim bottomlim toplim]);
axis square;
grid;
%%%%%%%%%%%%%%%%%

%Plot the original shape
plotlines(oldpoints,'g','r')

%The transformation multiplies each
%column of oldpoints (x,y) to get new points
%Rotate about origin by theta:
points = R1*oldpoints;
%Plot the transposed shape
plotlines(points,'b','r')

%Concatenate some transformations (right to left)
%inverse(T1) happens to move lower left corner to origin
points = T1*R1*inverse(T1)*oldpoints;
%Plot the transposed shape
plotlines(points,'b','k')