function [positions] = find_joint_positions(angles)
%[positions] = FIND_JOINT_POSITIONS(angles): 
%function to calculate the position of joints A through F for given angles
%
%Input angles = [ti t1 t2 t3 t4 t5 t6 t7 t8] given angles of Linkage system (rad)
%Output positions = [[JointA_x JointA_y] 
%                    [JointB_x JointB_y] 
%                    [JointC_x JointC_y] 
%                    [JointD_x JointD_y] 
%                    [JointE_x JointE_y]
%                    [JointF_x JointF_y]
%                    [Joint0_x Joint0_y] 
%                    [Joint1_x Joint1_y]]

%   Version 1: created 09/03/2017. Author: Conor Igoe
%   This MATLAB function M-file is not flexible. It works for the Jansen 
%   Linkage mechansim in problem 2 of MP1 only.

% -------------------------------------------------------------------------

% Check input and output arguments
if (nargin ~= 1), error('Incorrect number of input arguments.'); end
if (nargout ~= 1), error('Incorrect number of output arguments.'); end

% -------------------------------------------------------------------------

% Internal parameters li, l1, l2, l3, l4, l5, l6, l7, l8, l0 = linkage lengths (unitless)
% Internal parameter a = fixed joint vertical separation (unitless)
% Internal parameter b = fixed joint horizontal separation (unitless)
% Internal parameter c = fixed internal angle EDF (rad)

li = 15;
l1 = 50;
l2 = 41.5;
l3 = 55.8;
l4 = 40.1;
l5 = 39.4;
l6 = 61.9;
l7 = 39.3;
l8 = 36.7;
l9 = 49;

a = 7.8;
b = 38;

c = 1.729556;

% -------------------------------------------------------------------------

% Extract angles
ti = angles(1);
t1 = angles(2);
t2 = angles(3);
t3 = angles(4);
t4 = angles(5);
t5 = angles(6);
t6 = angles(7);
t7 = angles(8);
t8 = angles(9);

% Joint A
positions(1, 1)  = li*cos(ti);
positions(1, 2)  = li*sin(ti);

% Joint B
positions(2, 1)  = l1*cos(t1) + positions(1, 1);
positions(2, 2)  = l1*sin(t1) + positions(1, 2);

% Joint C
positions(3, 1)  = l3*cos(t3) + positions(2, 1);
positions(3, 2)  = l3*sin(t3) + positions(2, 2);

% Joint D
positions(4, 1)  = -b + l7*cos(t7);
positions(4, 2)  = -a + l7*sin(t7);

% Joint E
positions(5, 1)  = l8*cos(t8) + positions(4, 1);
positions(5, 2)  = l8*sin(t8) + positions(4, 2);

% Joint F
positions(6, 1)  = l9*cos(t8 + c) + positions(4, 1);
positions(6, 2)  = l9*sin(t8 + c) + positions(4, 2);

% Joint 0
positions(7, 1) = 0;
positions(7, 2) = 0;

% Joint 1
positions(8, 1) = -b;
positions(8, 2) = -a;

end
