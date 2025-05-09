clear
close all
clc

run("advection_solver.m")
run("advection_solver_fd.m")
x_FD = x;
u_FD = u(:,end);
run("advection_solver_fem.m")
x_FEM = x;
u_FEM = u(:,end);

figure(2)
plot(xx(1,:),uu(1,:),'r-',xx(1,:),u_anal(1,:),'b:');
hold on
h1 = plot(xx',u_anal','b:', 'DisplayName','Analytical Method', LineWidth=2);
h3 = plot(x_FD,u_FD,'k-', 'DisplayName','FD Method', LineWidth=1.5);
h4 = plot(x_FEM, u_FEM,'m-', 'DisplayName', 'FEM Method', LineWidth=0.7);
h2 = plot(xx',uu','r-', 'DisplayName', 'DG Method', LineWidth=1.5);

xlabel('x')
ylabel('u_h(x)')