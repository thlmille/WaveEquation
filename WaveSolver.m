%
% Thomas Miller thlmille
% Final project for AMS 147
% FDwavEq.m
%   This file uses a finite difference numerical method to
%     solve the wave equation in one dimension.
%
L = pi;
time = .01;
c = 1;

tstep = time/512;
Lstep = L/512;

tvals = [0:tstep:time];
xvals = [0:Lstep:L];

vecsize = size(xvals);
veclength = vecsize(2);

u = zeros(veclength, veclength);

%
% Fill in boundaries
%
% u(0, 0) -> u(L, 0)
for i = 1:veclength
    u(i, 1) = f(xvals(i));
end

% u(0, 0) -> u(0, time)
for i = 1:veclength
   u(1, i) = A(tvals(i)); 
end

% u(L, 0) -> u(L, time)
for i = 1:veclength
   u(veclength, i) = A(tvals(i));
end


%
% Use finite difference method to fill u matrix that represents
%   solution u(x, t)
%
for j = 2:veclength-1
    for k = 2:veclength-1
        numer = c^2*u(j+1, k-1) - 2*u(j, k-1) + u(j-1, k-1);        
        u(j, k) = tstep^2 * numer/(Lstep^2) + 2*u(j, k-1) - u(j-1, k-1);
    end
end

%mesh (xvals, tvals, u)

%plot(xvals, u(1:veclength, 200))


for k = 1:10:veclength
    plot(xvals, u(1:veclength, k));
    title (num2str(tvals(k)))
    %axis([0, pi, -2, 2])
    M(k) = getframe;
end

clf
movie(M,1)



