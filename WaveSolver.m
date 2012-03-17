%
% Thomas Miller thlmille
% Final project for AMS 147
% FDwavEq.m
%   This file uses a finite difference numerical method to
%     solve the wave equation in one dimension.
%
L = 10;
time = 5;
c = 1;

tstep = time/256;
Lstep = L/256;

tvals = [0:tstep:time];
xvals = [0:Lstep:L];

vecsize = size(xvals);
xveclength = vecsize(2);
vecsize = size(tvals);
tveclength = vecsize(2);

u = zeros(xveclength, tveclength);

%
% Fill in boundaries
%
% u(0, 0) -> u(L, 0)
for i = 1:xveclength
    u(i, 1) = f(xvals(i));
end

% u(0, 0) -> u(0, time)
for i = 1:tveclength
   u(1, i) = A(tvals(i)); 
end

% u(L, 0) -> u(L, time)
for i = 1:tveclength
   u(xveclength, i) = A(tvals(i));
end


%
% Use finite difference method to fill u matrix that represents
%   solution u(x, t)
%
for k = 2:tveclength-1
    for j = 2:xveclength-1
        numer = c^2*u(j+1, k-1) - 2*u(j, k-1) + u(j-1, k-1);        
        u(j, k) = tstep^2 * numer/(Lstep^2) + 2*u(j, k-1) - u(j-1, k-1);
    end
end

%mesh (xvals, tvals, u)

%plot(xvals, u(1:veclength, 1))

h = 1;
for k = 1:10:tveclength
    plot(xvals, u(1:xveclength, k));
    title (num2str(tvals(k)))
    %axis([0, pi, -2, 2])
    M(h) = getframe;
    h = h + 1;
end

clf
movie(M, 3, 10)



