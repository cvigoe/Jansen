function [tpositions, trace] = plot_Jansen(ti, t_estimates, a_estimates, input_angle_step_size, cycles)
%[positions] = PLOT_JANSEN(angles): 
%function to plot the 4-Legged Jansen Mechanism and trace of joint P from a 
%starting crank angle over a number of crank angles, as determined by the
%step size and number of cycles
%
%Input ti = starting crank angle (rad)
%Input t_estimates = [t1_0 t2_0 t3_0 t4_0 t5_0 t6_0 t7_0 t8_0] theta angle initial estimates (rad)
%Input a_estimates = [a1_0 a2_0 a3_0 a4_0 a5_0 a6_0 a7_0 a8_0] alpha angle initial estimates (rad)
%Input input_angle_step_size = step size (rad)
%Input cycles = number of desired cycles to simulate pseudo-statically 
%Output positions = [[tJointA_x tJointA_y] 
%                    [tJointB_x tJointB_y] 
%                    [tJointC_x tJointC_y] 
%                    [tJointD_x tJointD_y] 
%                    [tJointE_x tJointE_y]
%                    [tJointF_x tJointF_y]
%                    [tJoint0_x tJoint0_y] 
%                    [tJoint1_x tJoint1_y]] final positions for theta leg joints A:F , 0 & 1

%   Version 1: created 09/03/2017. Author: Conor Igoe
%   This MATLAB function M-file is not flexible. It works for the Jansen 
%   Linkage mechansim in problem 2 of MP1 only.

% -------------------------------------------------------------------------

% Check input and output arguments
% Argument check
if (nargin ~= 5), error('Incorrect number of input arguments.'); end
if (nargout ~= 2), error('Incorrect number of output arguments.'); end

% Close all open figures
close all

% Assign inital estimates
t1_0 = t_estimates(1);
t2_0 = t_estimates(2);
t3_0 = t_estimates(3);
t4_0 = t_estimates(4);
t5_0 = t_estimates(5);
t6_0 = t_estimates(6);
t7_0 = t_estimates(7);
t8_0 = t_estimates(8);

w1_0 = t1_0;
w2_0 = t2_0;
w3_0 = t3_0;
w4_0 = t4_0;
w5_0 = t5_0;
w6_0 = t6_0;
w7_0 = t7_0;
w8_0 = t8_0;  

a1_0 = a_estimates(1);
a2_0 = a_estimates(2);
a3_0 = a_estimates(3);
a4_0 = a_estimates(4);
a5_0 = a_estimates(5);
a6_0 = a_estimates(6);
a7_0 = a_estimates(7);
a8_0 = a_estimates(8);

b1_0 = a1_0;
b2_0 = a2_0;
b3_0 = a3_0;
b4_0 = a4_0;
b5_0 = a5_0;
b6_0 = a6_0;
b7_0 = a7_0;
b8_0 = a8_0;       

% Preallocate trace matrix
trace = zeros(2, floor(2*cycles*pi/input_angle_step_size));

% Counter for trace indexing
count = 1;

% Assign input angles
winput = ti;
ainput = ti + pi;
binput = ti + pi;

% Loop through input angles, 
for tinput = ti:input_angle_step_size:ti + 2*cycles*pi + input_angle_step_size
    % Calculate angles using NR
    [t1, t2, t3, t4, t5, t6, t7, t8] = Jansen_Sequential_Newton_Raphson(tinput, t1_0, t2_0, t3_0, t4_0, t5_0, t6_0, t7_0, t8_0);
    tangles = [tinput, t1, t2, t3, t4, t5, t6, t7, t8];    
    
    [w1, w2, w3, w4, w5, w6, w7, w8] = Jansen_Sequential_Newton_Raphson(winput, w1_0, w2_0, w3_0, w4_0, w5_0, w6_0, w7_0, w8_0);
    wangles = [winput, w1, w2, w3, w4, w5, w6, w7, w8];
    
    [a1, a2, a3, a4, a5, a6, a7, a8] = Jansen_Sequential_Newton_Raphson(ainput, a1_0, a2_0, a3_0, a4_0, a5_0, a6_0, a7_0, a8_0);
    aangles = [ainput, a1, a2, a3, a4, a5, a6, a7, a8];    
    
    [b1, b2, b3, b4, b5, b6, b7, b8] = Jansen_Sequential_Newton_Raphson(binput, b1_0, b2_0, b3_0, b4_0, b5_0, b6_0, b7_0, b8_0);
    bangles = [binput, b1, b2, b3, b4, b5, b6, b7, b8];        
    
    % Find positions of nodes from angles
    tpositions = find_joint_positions(tangles);
    wpositions = find_joint_positions(wangles);
    apositions = find_joint_positions(aangles);
    bpositions = find_joint_positions(bangles);

    tjointA = [tpositions(1,1) , tpositions(1,2)];
    tjointB = [tpositions(2,1) , tpositions(2,2)];
    tjointC = [tpositions(3,1) , tpositions(3,2)];
    tjointD = [tpositions(4,1) , tpositions(4,2)];
    tjointE = [tpositions(5,1) , tpositions(5,2)];
    tjointF = [tpositions(6,1) , tpositions(6,2)];
    tjoint0 = [tpositions(7,1) , tpositions(7,2)];
    tjoint1 = [tpositions(8,1) , tpositions(8,2)];
    
    wjointA = [-wpositions(1,1) , wpositions(1,2)];
    wjointB = [-wpositions(2,1) , wpositions(2,2)];
    wjointC = [-wpositions(3,1) , wpositions(3,2)];
    wjointD = [-wpositions(4,1) , wpositions(4,2)];
    wjointE = [-wpositions(5,1) , wpositions(5,2)];
    wjointF = [-wpositions(6,1) , wpositions(6,2)];
    wjoint0 = [-wpositions(7,1) , wpositions(7,2)];
    wjoint1 = [-wpositions(8,1) , wpositions(8,2)];    
    
    ajointA = [apositions(1,1) , apositions(1,2)];
    ajointB = [apositions(2,1) , apositions(2,2)];
    ajointC = [apositions(3,1) , apositions(3,2)];
    ajointD = [apositions(4,1) , apositions(4,2)];
    ajointE = [apositions(5,1) , apositions(5,2)];
    ajointF = [apositions(6,1) , apositions(6,2)];
    ajoint0 = [apositions(7,1) , apositions(7,2)];
    ajoint1 = [apositions(8,1) , apositions(8,2)];     
    
    bjointA = [-bpositions(1,1) , bpositions(1,2)];
    bjointB = [-bpositions(2,1) , bpositions(2,2)];
    bjointC = [-bpositions(3,1) , bpositions(3,2)];
    bjointD = [-bpositions(4,1) , bpositions(4,2)];
    bjointE = [-bpositions(5,1) , bpositions(5,2)];
    bjointF = [-bpositions(6,1) , bpositions(6,2)];
    bjoint0 = [-bpositions(7,1) , bpositions(7,2)];
    bjoint1 = [-bpositions(8,1) , bpositions(8,2)];      
    
    % In-loop plotting; clear previous plot; set axis for correct scaling
    drawnow
    clf        
    axis([-115 115 -115 115])
    
    % Plot legs "furthest" from side viewpoint first
    patch('Faces', [1 2 3], 'Vertices', [ajoint1; ajointB; ajointC],'FaceColor',[0.8 0.8 0.8])
    patch('Faces', [1 2 3], 'Vertices', [ajointD; ajointE; ajointF],'FaceColor',[0.8 0.8 0.8])
     
    line([ajoint0(1) ajointA(1)], [ajoint0(2) ajointA(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([ajointA(1) ajointB(1)], [ajointA(2) ajointB(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajointB(1) ajointC(1)], [ajointB(2) ajointC(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajointE(1) ajointD(1)], [ajointE(2) ajointD(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajointA(1) ajointD(1)], [ajointA(2) ajointD(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajointB(1) ajoint1(1)], [ajointB(2) ajoint1(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajoint1(1) ajointD(1)], [ajoint1(2) ajointD(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajoint1(1) ajointC(1)], [ajoint1(2) ajointC(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([ajointE(1) ajointF(1)], [ajointE(2) ajointF(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajointD(1) ajointF(1)], [ajointD(2) ajointF(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([ajointC(1) ajointE(1)], [ajointC(2) ajointE(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)    
    
    % Plot trace inbetween "far" legs and "near" legs
    trace(1,count) = tjointF(1);
    trace(2,count) = tjointF(2);
    line(trace(1,1:count), trace(2,1:count));    
    
    patch('Faces', [1 2 3], 'Vertices', [bjoint1; bjointB; bjointC], 'FaceColor',[0.8 0.8 0.8])
    patch('Faces', [1 2 3], 'Vertices', [bjointD; bjointE; bjointF], 'FaceColor',[0.8 0.8 0.8])
     
    line([bjoint0(1) bjointA(1)], [bjoint0(2) bjointA(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([bjointA(1) bjointB(1)], [bjointA(2) bjointB(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjointB(1) bjointC(1)], [bjointB(2) bjointC(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjointE(1) bjointD(1)], [bjointE(2) bjointD(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjointA(1) bjointD(1)], [bjointA(2) bjointD(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjointB(1) bjoint1(1)], [bjointB(2) bjoint1(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjoint1(1) bjointD(1)], [bjoint1(2) bjointD(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjoint1(1) bjointC(1)], [bjoint1(2) bjointC(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([bjointE(1) bjointF(1)], [bjointE(2) bjointF(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjointD(1) bjointF(1)], [bjointD(2) bjointF(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([bjointC(1) bjointE(1)], [bjointC(2) bjointE(2)], 'Color',[0.8 0.8 0.8], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)        
    
    patch('Faces', [1 2 3], 'Vertices', [tjoint1; tjointB; tjointC])
    patch('Faces', [1 2 3], 'Vertices', [tjointD; tjointE; tjointF])

    line([tjoint0(1) tjointA(1)], [tjoint0(2) tjointA(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([tjointA(1) tjointB(1)], [tjointA(2) tjointB(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjointB(1) tjointC(1)], [tjointB(2) tjointC(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjointE(1) tjointD(1)], [tjointE(2) tjointD(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjointA(1) tjointD(1)], [tjointA(2) tjointD(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjointB(1) tjoint1(1)], [tjointB(2) tjoint1(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjoint1(1) tjointD(1)], [tjoint1(2) tjointD(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjoint1(1) tjointC(1)], [tjoint1(2) tjointC(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([tjointE(1) tjointF(1)], [tjointE(2) tjointF(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjointD(1) tjointF(1)], [tjointD(2) tjointF(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([tjointC(1) tjointE(1)], [tjointC(2) tjointE(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)

    patch('Faces', [1 2 3], 'Vertices', [wjoint1; wjointB; wjointC])
    patch('Faces', [1 2 3], 'Vertices', [wjointD; wjointE; wjointF])
     
    line([wjoint0(1) wjointA(1)], [wjoint0(2) wjointA(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([wjointA(1) wjointB(1)], [wjointA(2) wjointB(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjointB(1) wjointC(1)], [wjointB(2) wjointC(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjointE(1) wjointD(1)], [wjointE(2) wjointD(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjointA(1) wjointD(1)], [wjointA(2) wjointD(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjointB(1) wjoint1(1)], [wjointB(2) wjoint1(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjoint1(1) wjointD(1)], [wjoint1(2) wjointD(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjoint1(1) wjointC(1)], [wjoint1(2) wjointC(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0])
    line([wjointE(1) wjointF(1)], [wjointE(2) wjointF(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjointD(1) wjointF(1)], [wjointD(2) wjointF(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)
    line([wjointC(1) wjointE(1)], [wjointC(2) wjointE(2)], 'Color',[0.5 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 40)

    % Update initial estimates for NR for next input angles
    t1_0 = t1;
    t2_0 = t2;
    t3_0 = t3;
    t4_0 = t4;
    t5_0 = t5;
    t6_0 = t6;
    t7_0 = t7;
    t8_0 = t8;
    
    w1_0 = w1;
    w2_0 = w2;
    w3_0 = w3;
    w4_0 = w4;
    w5_0 = w5;
    w6_0 = w6;
    w7_0 = w7;
    w8_0 = w8;    
    
    a1_0 = a1;
    a2_0 = a2;
    a3_0 = a3;
    a4_0 = a4;
    a5_0 = a5;
    a6_0 = a6;
    a7_0 = a7;
    a8_0 = a8;    

    b1_0 = b1;
    b2_0 = b2;
    b3_0 = b3;
    b4_0 = b4;
    b5_0 = b5;
    b6_0 = b6;
    b7_0 = b7;
    b8_0 = b8;          

    % Incremenet counter for trace indexing
    count = count + 1;
    
    % Update input angles
    winput = winput - input_angle_step_size;
    ainput = ainput + input_angle_step_size;
    binput = binput - input_angle_step_size;
end
end
