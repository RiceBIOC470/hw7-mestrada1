%HW7

% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows down as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% At X = 0, the population is 0, so the rate of growth/division (derivative
% of population growth function) will be equal to 0, and thus, X = 0 remains
% a fixed point. At X = 1, the population is maxed out (at N), so the rate 
% of growth/division (derivative of population growth function) will also be
% equal to 0, and thus, X = 1 remains a fixed point.

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

a = 2; % Smaller a parameter
X = 0:0.1:2;
plot(X, a.*X.*(1-X), 'LineWidth', 2); hold on;
plot([0 2], [0 0], 'k', 'LineWidth', 2);
xlabel('X'); ylabel('g(x)'); set(gca, 'FontSize', 24);
hold off;

figure;
a = 5; % Larger a parameter
X = 0:0.1:2;
plot(X, a.*X.*(1-X), 'LineWidth', 2); hold on;
plot([0 2], [0 0], 'k', 'LineWidth', 2);
xlabel('X'); ylabel('g(x)'); set(gca, 'FontSize', 24);
hold off;

% The stability of the fixed points does not depend on the parameter a;
% the point at X = 0 will always be unstable and the point at X = 1 will
% always be stable. The a parameter simply indicates how quickly the system
% will approach the stable steady state at X = 1. 

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 

% See eq_integ function code in repository

% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

x0 = rand(1,200);
for a = 1:4
    for ii = 1:200
        xf = x0(ii);
        for jj = 1:200
            xf = a * xf * (1-xf);
        end
        plot(a, xf, 'r.', 'MarkerSize', 24);
        hold on;
    end
end
xlabel('a'); ylabel('xf');

% When a = 1, the only steady state value is found at X(t) = 0. When a = 2,
% the only steady state value is found at X(t) = 0.5. When a = 3, the
% steady state value is found between 0.6 and 0.7. When a = 4, the steady
% state value should be 0.75, but at a = 4, the system is pretty unstable
% and hence, shows a range of fixed point values.

% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 

% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 

% dA/dt = (V+B^4)./(1+B^4) - A
% dB/dt = (V+A^4)./(1+A^4) - B

% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 

% V = 5; A0 > B0

V = 5;
A0 = 2;
B0 = 1;
rhs2 = @(t, x) [(V+x(2)^4)./(1+x(2)^4) - x(1) ; (V+x(1)^4)./(1+x(1)^4) - x(2)];
sol2 = ode23(rhs2, [0 10], [A0 B0]); 
plot(sol2.x, sol2.y(1,:), 'r.-', 'LineWidth', 3, 'MarkerSize', 18); hold on;
plot(sol2.x, sol2.y(2,:), 'g.-', 'LineWidth', 3, 'MarkerSize', 18);
xlabel('Time'); ylabel('Expression');
set(gca, 'FontSize', 24);
legend({'A', 'B'});

% V = 5; B0 > A0

V = 5;
A0 = 1;
B0 = 2;
rhs2 = @(t, x) [(V+x(2)^4)./(1+x(2)^4) - x(1) ; (V+x(1)^4)./(1+x(1)^4) - x(2)];
sol2 = ode23(rhs2, [0 10], [A0 B0]); 
plot(sol2.x, sol2.y(1,:), 'r.-', 'LineWidth', 3, 'MarkerSize', 18); hold on;
plot(sol2.x, sol2.y(2,:), 'g.-', 'LineWidth', 3, 'MarkerSize', 18);
xlabel('Time'); ylabel('Expression');
set(gca, 'FontSize', 24);
legend({'A', 'B'});

% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 

gx = @(x, V) (V*x^4)/(1+x^4)-x;
for V = 0:0.1:10
    gx2 = @(x) gx(x,V);
    for x0 = 0:0.1:3
        [rt, ~, exitflag] = fzero(gx2, x0);
        if exitflag == 1
            plot(V, rt, 'k.');
            hold on;
        end 
    end 
end