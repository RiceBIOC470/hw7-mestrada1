function [timecourse, t99] = eq_integ(x0, a)
rhs = @(t,x) a.*x.*(1-x);
dt = 0.01; 
interval = [0 10];
nstep = (interval(2)-interval(1))/dt;
sol1(1) = x0;
for ii = 2:nstep
    sol1(ii) = sol1(ii-1)+rhs(ii, sol1(ii-1))*dt;
end
tt = linspace(interval(1), interval(2), nstep);
timecourse = [sol1; tt];
t99 = tt(find(tt>0.99, 1));