% Geometry Plot
figure
for i = 1:number_of_spanwise_divisions
    for j = 1:number_of_chordwise_divisions
        plot3(panel_control_points(i,j,1), panel_control_points(i,j,2),0, 'k.');
        hold on
        %plot3(panel_vortex_points(i,j,1), panel_vortex_points(i,j,2), 0, 'b.');
        plot3(-panel_control_points(i,j,1), panel_control_points(i,j,2),0, 'k.');
        %plot3(-panel_vortex_points(i,j,1), panel_vortex_points(i,j,2), 0, 'b.');
        
    end
end
axis equal
%fplot(@(x) tan(sweep*pi/180)*x + rc/4, [0,b/2], 'b') % Quarter Chord Line
fplot(@(x) tan(sweep*pi/180)*x + rc/4 + 0.75*(rc-((rc-tc)/(b/2))*x), [0,b/2], 'r')  % Trailing Edge
fplot(@(x) tan(sweep*pi/180)*x + rc/4 - 0.25*(rc-((rc-tc)/(b/2))*x), [0,b/2], 'r')  % Leading Edge
fplot(@(x) tan(sweep*pi/180)*-x + rc/4 + 0.75*(rc-((rc-tc)/(b/2))*-x), [-b/2,0], 'r')  % Trailing Edge
fplot(@(x) tan(sweep*pi/180)*-x + rc/4 - 0.25*(rc-((rc-tc)/(b/2))*-x), [-b/2,0], 'r')  % Leading Edge
txt = ['Root Chord = ' num2str(rc) 'm'];
text(0.5, -1, txt);
txt = ['Tip Chord = ' num2str(tc) 'm'];
text(0.5, -1.5, txt);
txt = ['Quarter Chord Sweep = ' num2str(sweep) ' degrees'];
text(0.5, -2.0, txt);
txt = ['Airfoil = ' num2str(NACA)];
text(0.5, -2.5, txt);
title('Geometry Plot');

%%
%Delta Cp Plot
figure
colormap Hot;
L1 = 0;
L2 = 0;
D1 = 0;
D2 = 0;
area1 = 0;
area2 = 0;
for i = 1:number_of_spanwise_divisions
    for j = 1:number_of_chordwise_divisions
        subplot(2,1,1);
        hold on
        a = [panel_edges_x(i), panel_edges_x(i+1), panel_edges_x(i+1), panel_edges_x(i)];
        y = [root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i)*2/b), root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b), root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b), root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i)*2/b)];       
        fill(a,y, Del_Cp1((j-1)*number_of_spanwise_divisions+i));
        a = [-panel_edges_x(i), -panel_edges_x(i+1), -panel_edges_x(i+1), -panel_edges_x(i)];
        fill(a,y, Del_Cp1((j-1)*number_of_spanwise_divisions+i));
        
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b);
        areap = (b/(2*number_of_spanwise_divisions))*(trailing_y_edge-leading_y_edge)/2;
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i)*2/b);
        areap = areap + (b/(2*number_of_spanwise_divisions))*(trailing_y_edge-leading_y_edge)/2;
        D1 = D1 + 1*Del_Cp1((j-1)*number_of_spanwise_divisions+i)*areap*sin(panel_slope(i,j,1)-AOA*pi/180);
        L1 = L1 + 2*Del_Cp1((j-1)*number_of_spanwise_divisions+i)*areap*cos(panel_slope(i,j,1)-AOA*pi/180);
        area1 = area1 + 2*areap;

        subplot(2,1,2);
        hold on
        a = [panel_edges_x(i), panel_edges_x(i+1), panel_edges_x(i+1), panel_edges_x(i)];
        y = [root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i)*2/b), root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b), root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b), root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i)*2/b)];       
        fill(a,y, Del_Cp2((j-1)*number_of_spanwise_divisions+i));
        a = [-panel_edges_x(i), -panel_edges_x(i+1), -panel_edges_x(i+1), -panel_edges_x(i)];
        fill(a,y, Del_Cp2((j-1)*number_of_spanwise_divisions+i));
        
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i+1)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i+1)*2/b);
        areap = (b/(2*number_of_spanwise_divisions))*(trailing_y_edge-leading_y_edge)/2;
        leading_y_edge  = root_chord_y(j) + (tip_chord_y(j)-root_chord_y(j))*(panel_edges_x(i)*2/b);
        trailing_y_edge = root_chord_y(j+1) + (tip_chord_y(j+1)-root_chord_y(j+1))*(panel_edges_x(i)*2/b);
        areap = areap + (b/(2*number_of_spanwise_divisions))*(trailing_y_edge-leading_y_edge)/2;
        D2 = D2 + 1*Del_Cp2((j-1)*number_of_spanwise_divisions+i)*areap*sin(panel_slope(i,j,1)-AOA*pi/180);
        L2 = L2 + 2*Del_Cp2((j-1)*number_of_spanwise_divisions+i)*areap*cos(panel_slope(i,j,1)-AOA*pi/180);
        area2 = area2 + 2*areap;
    end
end
Cd1 = D1/area1;
Cl1 = -L1/(area1);
D1 = D1*(0.5*ro*U^2);
Cd2 = D2/area1;
Cl2 = -L2/(area2);
D2 = D2*(0.5*ro*U^2);
subplot(2,1,1);
axis equal
subplot(2,1,2);
axis equal

subplot(2,1,1);
%{
txt = ['Freestream Speed = ' num2str(U) 'm/s'];
text(-0.5*b/2, -1, txt);
txt = ['AOA = ' num2str(AOA) ' degrees'];
text(-0.5*b/2, -1.5, txt);
txt = ['Actuator Vortex Strength = ' num2str(0) 'm^2/s'];
text(-0.5*b/2, -2.0, txt);
txt = ['Total Lift = ' num2str(L1) 'N'];
text(-0.5*b/2, -3.0, txt);
txt = ['Total Drag = ' num2str(D1) 'N'];
text(-0.5*b/2, -3.5, txt);
txt = ['Cl = ' num2str(Cl1)];
text(-0.5*b/2, -4.0, txt);
txt = ['Cd = ' num2str(Cd1)];
text(-0.5*b/2, -4.5, txt);
%}
txt = ['Cl = ' num2str(Cl1)];
text(-0.75*b/2, -1, txt);
txt = ['Cd = ' num2str(Cd1)];
text(+0.5*b/2, -1, txt);
title('Delta Cp Distribution Without Actuator Vortex');
colorbar

subplot(2,1,2);
%{
txt = ['Actuator Vortex Strength = ' num2str(Strength_of_Actuator_Vortex) 'm^2/s'];
text(-0.5*b/2, -1.0, txt);
txt = ['Total Lift = ' num2str(L2) 'N'];
text(-0.5*b/2, -1.5, txt);
txt = ['Total Drag = ' num2str(D2) 'N'];
text(-0.5*b/2, -2.0, txt);
txt = ['Cl = ' num2str(Cl2)];
text(-0.5*b/2, -2.5, txt);
txt = ['Cd = ' num2str(Cd2)];
text(-0.5*b/2, -3.0, txt);
%}
txt = ['Cl = ' num2str(Cl2)];
text(-0.75*b/2, -1, txt);
txt = ['Cd = ' num2str(Cd2)];
text(+0.5*b/2, -1, txt);
txt = ['AVS = ' num2str(Strength_of_Actuator_Vortex)];
text(0, -1, txt);
title('Delta Cp Distribution With Actuator Vortex');
colorbar
%%


% Spanwise Cl Plot
figure
%plot(panel_midpoint_x/b, Spanwise_Cl1, 'k-');
hold on
%plot(panel_midpoint_x/b, Spanwise_Cl2, 'r-');
%plot(-panel_midpoint_x/b, Spanwise_Cl1, 'k-');
%plot(-panel_midpoint_x/b, Spanwise_Cl2, 'r-');

plot(panel_midpoint_x/b, Spanwise_Cl2-Spanwise_Cl1, 'b-');
plot(-panel_midpoint_x/b, Spanwise_Cl2-Spanwise_Cl1, 'b-');
legend('Black: Without Actuator', 'RED: With Actuator');
xlabel('Spanwise Location')
ylabel('Cl')
title('Spanwise Coefficient of Lift')

%%
% Spanwise Moment Plot
figure
plot(panel_midpoint_x/b, Spanwise_Moment1, 'k-');
hold on
plot(panel_midpoint_x/b, Spanwise_Moment2, 'r-');
plot(-panel_midpoint_x/b, Spanwise_Moment1, 'k-');
plot(-panel_midpoint_x/b, Spanwise_Moment2, 'r-');
legend('Black: Without Actuator', 'RED: With Actuator');
xlabel('Spanwise Location')
ylabel('Cm')
title('Spanwise Coefficient of Moment')





