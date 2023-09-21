b             = 18  ;    % wingspan in meters
rc            = 4.0 ;    % root chord in meters
tc            = 1.0 ;    % tip chord in meters
sweep         = 5  ;    % in degrees
AOA           = 0   ;    % in degrees
U             = 400  ;    % in metres/second
ro            = 1.225;    % density of fluid in kg/m3
Strength_of_Actuator_Vortex = 100;

number_of_spanwise_divisions    = 40; 
number_of_chordwise_divisions   = 5;

NACA = 4415;
T = rem(NACA,100);
P = rem((NACA-T)/100,10);
M = ((NACA-T)/100 - P)/10;
P = 0.1*P;
M = 0.01*M;
T = 0.01*T;


root_chord_y     = linspace(0, rc, number_of_chordwise_divisions+1);
tip_chord_y      = linspace(tan(sweep*pi/180)*b/2 + rc/4 - 0.25*(rc-((rc-tc)/(b/2))*b/2), tan(sweep*pi/180)*b/2 + rc/4 + 0.75*(rc-((rc-tc)/(b/2))*b/2), number_of_chordwise_divisions+1);
panel_edges_x    = linspace(0, b/2, number_of_spanwise_divisions+1);
panel_midpoint_x = linspace(b/(4*number_of_spanwise_divisions), b/2 - b/(4*number_of_spanwise_divisions), number_of_spanwise_divisions);


panel_control_points = zeros(number_of_spanwise_divisions,number_of_chordwise_divisions, 2);
panel_vortex_points  = zeros(number_of_spanwise_divisions,number_of_chordwise_divisions, 2);
panel_slope          = zeros(number_of_spanwise_divisions,number_of_chordwise_divisions, 2);

%figure
%exxageration_factor = 5;
for i = 1:number_of_spanwise_divisions
    chord_length_for_this_strip = rc - panel_midpoint_x(i)*(rc-tc)*2/b;
    leading_edge_y = tan(sweep*pi/180)*panel_midpoint_x(i) + rc/4 - 0.25*(rc-((rc-tc)/(b/2))*panel_midpoint_x(i));
    for j = 1:number_of_chordwise_divisions

        panel_control_points(i,j,1) = panel_midpoint_x(i);
        panel_vortex_points(i,j,1)  = panel_edges_x(i+1);
        
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_midpoint_x(i)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_midpoint_x(i)*2/b);
        panel_control_points(i,j,2) = (3*trailing_y_edge + leading_y_edge)/4;

        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b);
        panel_vortex_points(i,j,2)  = (3*leading_y_edge + trailing_y_edge)/4;
        
        y = (panel_control_points(i,j,2)-leading_edge_y)/chord_length_for_this_strip;
        retval = func(y, P, M, T);
        panel_slope(i,j,1) = retval(5);
        y = (panel_vortex_points(i,j,2)-leading_edge_y)/chord_length_for_this_strip;
        retval = func(y, P, M, T);
        panel_slope(i,j,2) = retval(5);

        %scatter(panel_control_points(i,j,1), panel_control_points(i,j,2), 5, 'k');
        %scatter3(panel_vortex_points(i,j,1), panel_vortex_points(i,j,2), 0, 10, 'b');
        %hold on
    end
end
%axis equal
%title('Control Points and Vortex Points')

% fplot(@(x) tan(sweep*pi/180)*x + rc/4, [0,b/2], 'b') % Quarter Chord Line
%fplot(@(x) tan(sweep*pi/180)*x + rc/4 + 0.75*(rc-((rc-tc)/(b/2))*x), [0,b/2], 'r')  % Trailing Edge
%fplot(@(x) tan(sweep*pi/180)*x + rc/4 - 0.25*(rc-((rc-tc)/(b/2))*x), [0,b/2], 'r')  % Leading Edge
%axis equal


Coeff_Matrix = zeros(number_of_spanwise_divisions*number_of_chordwise_divisions);
Tip_Actuator_Effect_Matrix = zeros(number_of_spanwise_divisions*number_of_chordwise_divisions,1);
B_Matrix = ones(number_of_spanwise_divisions*number_of_chordwise_divisions,1);


wakelength = rc*100; 
for i = 1:number_of_spanwise_divisions
    for j = 1:number_of_chordwise_divisions
        for k = 1:number_of_spanwise_divisions
            for l = 1:number_of_chordwise_divisions

                if k==1
                    vortex_left_point = [0, 0.25*(3*root_chord_y(l)+root_chord_y(l+1))];
                else
                    vortex_left_point = [panel_vortex_points(k-1,l,1),panel_vortex_points(k-1,l,2),0];
                end

                %{
                delta1 = panel_edges_x(k) - panel_control_points(i,j,1);
                delta2 = panel_edges_x(k+1) - panel_control_points(i,j,1);
                delta3 = panel_vortex_points(k,l,2) - panel_control_points(i,j,2);
                
                delta4 = -panel_edges_x(k) - panel_control_points(i,j,1);
                delta5 = -panel_edges_x(k+1) - panel_control_points(i,j,1);


                theta1 = atan(delta3/delta1);
                theta2 = atan(delta3/delta2);

                theta3 = atan(delta3/delta4);
                theta4 = atan(delta3/delta5);
                %}

                v = Line_Effect([vortex_left_point(1),vortex_left_point(2),0],[vortex_left_point(1),vortex_left_point(2)+wakelength,0],[panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
                v = v - Line_Effect([-vortex_left_point(1),vortex_left_point(2),0],[-vortex_left_point(1),vortex_left_point(2)+wakelength,0],[panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
                v = v - Line_Effect([panel_vortex_points(k,l,1),panel_vortex_points(k,l,2),0],[panel_vortex_points(k,l,1),panel_vortex_points(k,l,2)+wakelength,0],[panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
                v = v + Line_Effect([-panel_vortex_points(k,l,1),panel_vortex_points(k,l,2),0],[-panel_vortex_points(k,l,1),panel_vortex_points(k,l,2)+wakelength,0],[panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
                v = v + Line_Effect([panel_vortex_points(k,l,1),panel_vortex_points(k,l,2),0],[vortex_left_point(1),vortex_left_point(2),0],[panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
                v = v - Line_Effect([-panel_vortex_points(k,l,1),panel_vortex_points(k,l,2),0],[-vortex_left_point(1),vortex_left_point(2),0],[panel_control_points(i,j,1),panel_control_points(i,j,2),0]);

                Coeff_Matrix((j-1)*number_of_spanwise_divisions+i, (l-1)*number_of_spanwise_divisions+k) = v(3)*cos(panel_slope(i,j,1));


                
                
                %{
                if k>i
                    coeff = (1/(4*pi))*((1-sin(theta1))/delta1 - (1-sin(theta2))/delta2 + (cos(theta2) - cos(theta1))/delta3 + (1-sin(theta3))/delta4 - (1-sin(theta4))/delta5 + (cos(theta4) - cos(theta3))/delta3)*cos(panel_slope(i,j,1));
                elseif k<i
                    coeff = (1/(4*pi))*((1-sin(theta1))/delta1 - (1-sin(theta2))/delta2 + (cos(theta2) - cos(theta1))/delta3 - (1-sin(theta3))/delta4 + (1-sin(theta4))/delta5 + (cos(theta4) - cos(theta3))/delta3)*cos(panel_slope(i,j,1));
                elseif k==i
                    coeff = (1/(4*pi))*((1-sin(theta1))/delta1 - (1-sin(theta2))/delta2 + (cos(theta2) - cos(theta1))/delta3 - (1-sin(theta3))/delta4 - (1-sin(theta4))/delta5 + (cos(theta4) - cos(theta3))/delta3)*cos(panel_slope(i,j,1));
                end
                %}
                
                %{
                coeff = (1/(4*pi))*((1-sin(theta1))/delta1 - (1-sin(theta2))/delta2 + (cos(theta2) - cos(theta1))/delta3 + (1-sin(theta3))/delta4 - (1-sin(theta4))/delta5 + (cos(theta4) - cos(theta3))/delta3)*cos(panel_slope(i,j,1));
                Coeff_Matrix((j-1)*number_of_spanwise_divisions+i, (l-1)*number_of_spanwise_divisions+k) = coeff;
                %}
                

            end
        end
        delta4 = panel_edges_x(number_of_spanwise_divisions+1) - panel_control_points(i,j,1);
        delta5 = tip_chord_y(1) - panel_control_points(i,j,2);
        theta3 = atan(delta5/delta4);
        v = Line_Effect([panel_edges_x(number_of_spanwise_divisions+1),tip_chord_y(1),0],[panel_edges_x(number_of_spanwise_divisions+1),tip_chord_y(1)+wakelength,0], [panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
        v = v - Line_Effect([-panel_edges_x(number_of_spanwise_divisions+1),tip_chord_y(1),0],[-panel_edges_x(number_of_spanwise_divisions+1),tip_chord_y(1)+wakelength,0], [panel_control_points(i,j,1),panel_control_points(i,j,2),0]);
        Tip_Actuator_Effect_Matrix((j-1)*number_of_spanwise_divisions+i,1) = cos(panel_slope(i,j,1))*(Strength_of_Actuator_Vortex)*v(3);
        B_Matrix((j-1)*number_of_spanwise_divisions+i,1) = U*(sin(panel_slope(i,j,1)-AOA*pi/180));
    end
end

B2_Matrix = B_Matrix - Tip_Actuator_Effect_Matrix;
Coeff_Matrix;
B_Matrix;
%I = inv(Coeff_Matrix);
Solution1 = Coeff_Matrix\B_Matrix;
Solution2 = Coeff_Matrix\B2_Matrix;
%Solution = I * B_Matrix;

%{
%figure
for i = 1:number_of_spanwise_divisions
    for j  = 1:number_of_chordwise_divisions
        %scatter3(panel_vortex_points(i,j,1), panel_vortex_points(i,j,2), Solution((j-1)*number_of_spanwise_divisions+i,1), 50, "red");
        %scatter3(-panel_vortex_points(i,j,1), panel_vortex_points(i,j,2), Solution((j-1)*number_of_spanwise_divisions+i,1), 50, "red");
        %hold on
    end
end
%}

Panel_Lift1 = ro*U*Solution1*b/(2*number_of_spanwise_divisions);
Panel_Lift2 = ro*U*Solution2*b/(2*number_of_spanwise_divisions);
Cl1 = Panel_Lift1/(0.5*ro*U^2);
Cl2 = Panel_Lift2/(0.5*ro*U^2);
Del_Cp1 = zeros(number_of_spanwise_divisions*number_of_chordwise_divisions,1);
Del_Cp2 = zeros(number_of_spanwise_divisions*number_of_chordwise_divisions,1);
Spanwise_Cl1 = zeros(number_of_spanwise_divisions);
Spanwise_Cl2 = zeros(number_of_spanwise_divisions);
Spanwise_Moment1 = zeros(number_of_spanwise_divisions);
Spanwise_Moment2 = zeros(number_of_spanwise_divisions);
for i  = 1:number_of_spanwise_divisions
    spanCl1 = 0;
    spanCl2 = 0;
    spanM1 = 0;
    spanM2 = 0;
    for j = 1:number_of_chordwise_divisions
        spanCl1 = spanCl1 + Cl1((j-1)*number_of_spanwise_divisions+i);
        spanCl2 = spanCl2 + Cl2((j-1)*number_of_spanwise_divisions+i);
        spanM1 = spanM1 + Cl1((j-1)*number_of_spanwise_divisions+i)*panel_midpoint_x(i)*cos(AOA*pi/180);
        spanM2 = spanM2 + Cl2((j-1)*number_of_spanwise_divisions+i)*panel_midpoint_x(i)*cos(AOA*pi/180);
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b);
        areap = (b/(2*number_of_spanwise_divisions))*(trailing_y_edge-leading_y_edge)/2;
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i)*2/b);
        areap = areap + (b/(2*number_of_spanwise_divisions))*(trailing_y_edge-leading_y_edge)/2;
        Del_Cp1((j-1)*number_of_spanwise_divisions+i) = -Cl1((j-1)*number_of_spanwise_divisions+i)/(areap*cos(panel_slope(i,j,1)));
        Del_Cp2((j-1)*number_of_spanwise_divisions+i) = -Cl2((j-1)*number_of_spanwise_divisions+i)/(areap*cos(panel_slope(i,j,1)));
    end
    Spanwise_Cl1(i) = spanCl1;
    Spanwise_Cl2(i) = spanCl2;
    Spanwise_Moment1(i) = spanM1/b;
    Spanwise_Moment2(i) = spanM2/b;
end
 
%Del_Cp = Del_Cp/(0.5*ro*U^2);
%{
figure
hold on
for i = 1:number_of_spanwise_divisions
    for j = 1:number_of_chordwise_divisions
        %scatter(panel_control_points(i,j,1), panel_control_points(i,j,2), 50, Panel_Lift((j-1)*number_of_spanwise_divisions+i), 'filled');
        %scatter(-panel_control_points(i,j,1), panel_control_points(i,j,2), 50, Panel_Lift((j-1)*number_of_spanwise_divisions+i), 'filled');
        %scatter(panel_control_points(i,j,1), panel_control_points(i,j,2), 50, Del_Cp((j-1)*number_of_spanwise_divisions+i), 'filled');
        %scatter(-panel_control_points(i,j,1), panel_control_points(i,j,2), 50, Del_Cp((j-1)*number_of_spanwise_divisions+i), 'filled');
        a = [panel_edges_x(i), panel_edges_x(i+1), panel_edges_x(i+1), panel_edges_x(i)];
        y = [root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i)*2/b), root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b), root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b), root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i)*2/b)];
        fill(a,y, Cp((j-1)*number_of_spanwise_divisions+i));
        a = [-panel_edges_x(i), -panel_edges_x(i+1), -panel_edges_x(i+1), -panel_edges_x(i)];
        fill(a,y, Cp((j-1)*number_of_spanwise_divisions+i));
    end
end
colormap hot;
colorbar;
axis equal
%}




function y = func(x,p,m,t)
if x <p
    yc = (m/(p^2))*(2*p*x - x^2);
    slope = (m/(p^2))*(2*p - 2*x);
end

if x >=p
    yc = (m/(1-p)^2)*((1-2*p) + 2*p*x - x^2);
    slope = (m/(1-p)^2)*(2*p - 2*x);   
end

yt = (t/0.2)*(0.2969*x^(0.5) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4);

theta = atan(slope);

yu = yc + yt*cos(theta);
yl = yc - yt*cos(theta);

y = [yu,yl,yc,yt,theta];

end

function retval = Line_Effect(v1, v2, P)
% Effect of Line (v1-v2) on Control Point P
r1_c_r2_x = (P(2)-v1(2))*(P(3)-v2(3)) - (P(3)-v1(3))*(P(2)-v2(2));
r1_c_r2_y = -(P(1)-v1(1))*(P(3)-v2(3)) - (P(3)-v1(3))*(P(1)-v2(1));
r1_c_r2_z = (P(1)-v1(1))*(P(2)-v2(2)) - (P(2)-v1(2))*(P(1)-v2(1));
r1_c_r2_square = r1_c_r2_x^2 + r1_c_r2_y^2 + r1_c_r2_z^2;

r1 = ((P(1)-v1(1))^2 + (P(2)-v1(2))^2 + (P(3)-v1(3))^2)^0.5;
r2 = ((P(1)-v2(1))^2 + (P(2)-v2(2))^2 + (P(3)-v2(3))^2)^0.5;

rnot_dot_r1 = (v2(1)-v1(1))*(P(1)-v1(1)) + (v2(2)-v1(2))*(P(2)-v1(2)) + (v2(3)-v1(3))*(P(3)-v1(3)); 
rnot_dot_r2 = (v2(1)-v1(1))*(P(1)-v2(1)) + (v2(2)-v1(2))*(P(2)-v2(2)) + (v2(3)-v1(3))*(P(3)-v2(3)); 

K = (0.25/pi)*(rnot_dot_r1/r1 - rnot_dot_r2/r2)/r1_c_r2_square;

u = K*r1_c_r2_x;
v = K*r1_c_r2_y;
w = K*r1_c_r2_z;

retval = [u,v,w];
end