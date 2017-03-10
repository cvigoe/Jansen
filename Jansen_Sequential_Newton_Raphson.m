function [t1, t2, t3, t4, t5, t6, t7, t8] = Jansen_Sequential_Newton_Raphson(ti, t1_0, t2_0, t3_0, t4_0, t5_0, t6_0, t7_0, t8_0)
%[t1, t2, t3, t4, t5, t6, t7, t8] = JANSEN_SEQUENTIAL_NEWTON_RAPHSON(ti, t1_0, t2_0, t3_0, t4_0, t5_0, t6_0, t7_0, t8_0): 
%function to approximate the vector of unknowns X for the Jansen Linkage Mechanism in
%problem 2 of MP1 using the multi-dimensional Newton-Raphson method with
%initial estimates X_0 and crank angle ti with a sequential approach
%
%Input ti = crank angle (rad)
%Input t1_0 = initial estimate angle theta_1 (rad)
%Input t2_0 = initial estimate angle theta_2 (rad)
%Input t3_0 = initial estimate angle theta_3 (rad)
%Input t4_0 = initial estimate angle theta_4 (rad)
%Input t5_0 = initial estimate angle theta_5 (rad)
%Input t6_0 = initial estimate angle theta_6 (rad)
%Input t7_0 = initial estimate angle theta_7 (rad)
%Input t8_0 = initial estimate angle theta_8 (rad)
%Output t1 = approximate angle theta_1 (rad)
%Output t2 = approximate angle theta_2 (rad)
%Output t3 = approximate angle theta_3 (rad)
%Output t4 = approximate angle theta_4 (rad)
%Output t5 = approximate angle theta_5 (rad)
%Output t6 = approximate angle theta_6 (rad)
%Output t7 = approximate angle theta_7 (rad)
%Output t8 = approximate angle theta_8 (rad)

%   Version 1: created 09/03/2017. Author: Conor Igoe
%   This MATLAB function M-file is not flexible. It works for the Jansen Linkage mechansim
%   in problem 2 of MP1 only.
%
%   The limit on the number of allowable iterations and the acceptable 
%   tolerance for the Newton-Raphson (NR) algorithm are internally 
%   generated.

% -------------------------------------------------------------------------

% Check input and output arguments
if (nargin ~= 9), error('Incorrect number of input arguments.'); end
if (nargout ~= 8), error('Incorrect number of output arguments.'); end

% -------------------------------------------------------------------------

% Internal parameter ITERATION_LIMIT = maximum number of steps permitted.
% Internal parameter TOLERANCE = minimum acceptable value for |f_dash(x)|
%                                evaluated at the approximate roots.
% Internal parameters li, l1, l2, l3, l4, l5, l6, l7, l8 = linkage lengths (unitless)
% Internal parameter a = fixed joint vertical separation
% Internal parameter b = fixed joint horizontal separation

ITERATION_LIMIT = 10;
TOLERANCE = 10^-5;

a = 7.8;
b = 38;

lt = 15;
l1 = 50;
l2 = 41.5;
l3 = 55.8;
l4 = 40.1;
l5 = 39.4;
l6 = 61.9;
l7 = 39.3;
l8 = 36.7;

solutions = [0 0 0 0 0 0 0 0];

% -------------------------------------------------------------------------

% Theta1 and Theta2 pair NLEs
f1 = @(ti, t1, t2) ( lt*sin(ti) + l1*sin(t1) + a - l2*sin(t2) );
f2 = @(ti, t1, t2) ( lt*cos(ti) + l1*cos(t1) + b - l2*cos(t2) );

% Theta1 and Theta2 pair Jacobian
d_f1_t1 = @(t1, t2) ( l1*cos(t1) );
d_f1_t2 = @(t1, t2) ( -l2*cos(t2) );

d_f2_t1 = @(t1, t2) ( -l1*sin(t1) );
d_f2_t2 = @(t1, t2) ( l2*sin(t2) );

% Theta1 and Theta2 pair Newton-Raphson
X = [t1_0; t2_0];                       % Initialize initial guess
iterations = 0;                         % Initialize iterations counter
 
for count = 1:ITERATION_LIMIT + 1
    if count == ITERATION_LIMIT + 1
        error('The Iteration limit has been reached for Theta1 and Theta2 pair and did not converge')
        break
    end
     
    % Update root estimate
    t1 = X(1);
    t2 = X(2);
     
    % Check if we have met the desired tolerance
    if abs(f1(ti, t1, t2)) < TOLERANCE && abs(f2(ti, t1, t2)) < TOLERANCE
        break
    end

    % Calculate f(x)
    f = [f1(ti, t1, t2); f2(ti, t1, t2)];
    
    % Calculate D_f_dash(x)
    D_f =      [[d_f1_t1(t1, t2) d_f1_t2(t1, t2)];
                [d_f2_t1(t1, t2) d_f2_t2(t1, t2)]];
 
    % Calculate next approximation
    X = X - (D_f)\f;
    
    iterations = iterations + 1;
end

% Update solutions
solutions(1) = X(1);
solutions(2) = X(2);


% Theta3 and Theta4 pair NLEs
f3 = @(ti, t1, t3, t4) ( lt*sin(ti) + l1*sin(t1) + l3*sin(t3) + a - l4*sin(t4) );
f4 = @(ti, t1, t3, t4) ( lt*cos(ti) + l1*cos(t1) + l3*cos(t3) + b - l4*cos(t4) );

% Theta3 and Theta4 pair Jacobian
d_f3_t3 = @(t3, t4) ( l3*cos(t3) );
d_f3_t4 = @(t3, t4) ( -l4*cos(t4) );

d_f4_t3 = @(t3, t4) ( -l3*sin(t3) );
d_f4_t4 = @(t3, t4) ( l4*sin(t4) );

% Theta3 and Theta4 pair Newton-Raphson
X = [t3_0; t4_0];                       % Initialize initial guess
iterations = 0;                         % Initialize iterations counter
 
for count = 1:ITERATION_LIMIT + 1
    if count == ITERATION_LIMIT + 1
        error('The Iteration limit has been reached for Theta1 and Theta2 pair and did not converge')
        break
    end
     
    % Update root estimate
    t1 = solutions(1);
    t2 = solutions(2);
    t3 = X(1);
    t4 = X(2);    
     
    % Check if we have met the desired tolerance
    if abs(f3(ti, t1, t3, t4)) < TOLERANCE && abs(f4(ti, t1, t3, t4)) < TOLERANCE
        break
    end

    % Calculate f(x)
    f = [f3(ti, t1, t3, t4); f4(ti, t1, t3, t4)];
    
    % Calculate D_f_dash(x)
    D_f =      [[d_f3_t3(t3, t4) d_f3_t4(t3, t4)];
                [d_f4_t3(t3, t4) d_f4_t4(t3, t4)]];
 
    % Calculate next approximation
    X = X - (D_f)\f;
    
    iterations = iterations + 1;
end

% Update solutions
solutions(3) = X(1);
solutions(4) = X(2);


% Theta6 and Theta7 pair NLEs
f5 = @(ti, t6, t7) ( lt*sin(ti) + l6*sin(t6) + a - l7*sin(t7) );
f6 = @(ti, t6, t7) ( lt*cos(ti) + l6*cos(t6) + b - l7*cos(t7) );

% Theta6 and Theta7 pair Jacobian
d_f5_t6 = @(t6, t7) ( l6*cos(t6) );
d_f5_t7 = @(t6, t7) ( -l7*cos(t7) );

d_f6_t6 = @(t6, t7) ( -l6*sin(t6) );
d_f6_t7 = @(t6, t7) ( l7*sin(t7) );

% Theta6 and Theta7 pair Newton-Raphson
X = [t6_0; t7_0];                       % Initialize initial guess
iterations = 0;                         % Initialize iterations counter
 
for count = 1:ITERATION_LIMIT + 1
    if count == ITERATION_LIMIT + 1
        error('The Iteration limit has been reached for Theta1 and Theta2 pair and did not converge')
        break
    end
     
    % Update root estimate
    t6 = X(1);
    t7 = X(2);    
     
    % Check if we have met the desired tolerance
    if abs(f5(ti, t6, t7)) < TOLERANCE && abs(f6(ti, t6, t7)) < TOLERANCE
        break
    end

    % Calculate f(x)
    f = [f5(ti, t6, t7); f6(ti, t6, t7)];
    
    % Calculate D_f_dash(x)
    D_f =      [[d_f5_t6(t6, t7) d_f5_t7(t6, t7)];
                [d_f6_t6(t6, t7) d_f6_t7(t6, t7)]];
 
    % Calculate next approximation
    X = X - (D_f)\f;
    
    iterations = iterations + 1;
end

% Update solutions
solutions(6) = X(1);
solutions(7) = X(2);


% Theta5 and Theta8 pair NLEs
f7 = @(ti, t1, t3, t5, t7, t8) ( lt*sin(ti) + l1*sin(t1) + l3*sin(t3) + l5*sin(t5) + a - l7*sin(t7) - l8*sin(t8) );
f8 = @(ti, t1, t3, t5, t7, t8) ( lt*cos(ti) + l1*cos(t1) + l3*cos(t3) + l5*cos(t5) + b - l7*cos(t7) - l8*cos(t8) );

% Theta5 and Theta8 pair Jacobian
d_f7_t5 = @(t5, t8) ( l5*cos(t5) );
d_f7_t8 = @(t5, t8) ( -l8*cos(t8) );

d_f8_t5 = @(t5, t8) ( -l5*sin(t5) );
d_f8_t8 = @(t5, t8) ( l8*sin(t8) );

% Theta5 and Theta8 pair Newton-Raphson
X = [t5_0; t8_0];                       % Initialize initial guess
iterations = 0;                         % Initialize iterations counter
 
for count = 1:ITERATION_LIMIT + 1
    if count == ITERATION_LIMIT + 1
        error('The Iteration limit has been reached for Theta1 and Theta2 pair and did not converge')
        break
    end
     
    % Update root estimate
    t1 = solutions(1);
    t3 = solutions(3);
    t5 = X(1);
    t7 = solutions(7);
    t8 = X(2);
     
    % Check if we have met the desired tolerance
    if abs(f7(ti, t1, t3, t5, t7, t8)) < TOLERANCE && abs(f8(ti, t1, t3, t5, t7, t8)) < TOLERANCE
        break
    end

    % Calculate f(x)
    f = [f7(ti, t1, t3, t5, t7, t8); f8(ti, t1, t3, t5, t7, t8)];
    
    % Calculate D_f_dash(x)
    D_f =      [[d_f7_t5(t5, t8) d_f7_t8(t5, t8)];
                [d_f8_t5(t5, t8) d_f8_t8(t5, t8)]];
 
    % Calculate next approximation
    X = X - (D_f)\f;
    
    iterations = iterations + 1;
end
end
