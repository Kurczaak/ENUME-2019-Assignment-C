clear all
close all
%t=[0,10]
t=linspace(0,10,1000); %h=10/length(t)



%task1

%Graph of the Lobatto IIID solution
lobatto_solution=lobatto_solved_ode(t);
figure()
plot(t, lobatto_solution);
xlabel('t');
ylabel('y');
title('Lobatto IIID solution for the differential equation');
legend('y',"y'");

%Graph of the ode45 MatLab's function solution
ode45_solution=ode45_solved_ode(t);
figure()
plot(t,ode45_solution);
xlabel('t');
ylabel('y');
title('ODE45 solution for the differential equation');
legend('y',"y'");

%Graph of the Explicit Euler method solution
euler_solution=euler_solved_ode(t);
figure()
plot(t,euler_solution);
xlabel('t');
ylabel('y');
title('Explicit Euler solution for the differential equation');
legend('y',"y'");

%Visual comparison of the methods
%Not required in the assignment, just for my own curiosity
figure()
plot(t, euler_solution(1,:),'b');
hold on
plot(t,lobatto_solution(1,:),'g');
hold on
plot(t,ode45_solution(:,1),'r');
xlabel('t');
ylabel('y');
title('Visual comparison of the methods');
legend('Explicit Euler','Lobatto IIID','ODE45');


%task2 & task3
num_of_points=200; %How many values of h
h=logspace(-4,0,num_of_points);
%data preallocation
rms=zeros(1,num_of_points);
rms_euler=zeros(1,num_of_points);
mx=zeros(1,num_of_points);
mx_euler=zeros(1,num_of_points);
%counting errors
for i=1:length(h)
rms(i)=RMS(h(i));
mx(i)=MxError(h(i));
rms_euler(i)=RMS_euler(h(i));
mx_euler(i)=MxError_euler(h(i));
end

%RMS dependency on h
figure()
loglog(h,rms,'or');
hold on
loglog(h,rms_euler,'ob');
xlabel('h');
ylabel('\delta_2');
title('RMS Error dependency on h');
legend('Lobatto IIID', 'Explicit Euler');

%Maximum error dependency on h
figure()
loglog(h,mx,'or');
hold on
loglog(h,mx_euler,'ob');
xlabel('h');
ylabel('\delta_\infty');
title('Maximum Error dependency on h');
legend('Lobatto IIID', 'Explicit Euler');


%----functions-----

%Just for checking corectness of the solution
%at the very beginning of writing the code
%Not used in tasks completion
function y=manually_solved_ode(t)
 y=zeros(1,length(t));
 for i=1:(length(t))
   y(i)=2*exp(-t(i)/3)*sin(t(i));
 end
end

%Solving differential equation using ode45 MatLab function
function y=ode45_solved_ode(t)
    %9y''+6y'+10y=0
    y0=0;
    dy0dt=2;

    
    opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
    [t,y]=ode45( @rhs, t, [y0, dy0dt], opts);
    
    function dydt=rhs(t,y)
        dydt = [y(2); -(6*y(2) + 10*y(1))/9];
    end
end

%Solving differential equation using explicit Euler method
function y=euler_solved_ode(t)
y=zeros(2,length(t));
y(1,1)=0; %y0
y(2,1)=2; %dy0dt
h=10/length(t);   %step
A=[0,1;-10/9,-2/3];

for i=2:length(t)
    y(:,i)=y(:,i-1)+h.*A*y(:,i-1);
end
end

%Solving differential equation using LobattoIIID method
function y=lobatto_solved_ode(t)
F=zeros(6,length(t));
y=zeros(2,length(t));
y(1,1)=0; %y0
y(2,1)=2; %dy0dt
h=10/length(t);   %step

A=[0,1;-10/9,-2/3];
L=[eye(2)-1/6*h*A,zeros(2),1/6*h*A;
   -1/12*h*A, eye(2)-5/12*h*A, zeros(2);
   -1/2*h*A, -1/3*h*A, eye(2)-1/6*h*A];
for i=2:length(t)
    R=[A*y(:,i-1);A*y(:,i-1);A*y(:,i-1)];
    F(:,i-1)=1\L*R;
    f1=F(1:2,i-1);
    f2=F(3:4,i-1);
    f3=F(5:6,i-1);
    y(:,i)=y(:,i-1)+h*(1/6*f1+2/3*f2+1/6*f3);
end
end

%---------Errors calculation----

%Root mean square error
%Ode45 solution as the accurate solution
%Lobatto IIID as the approximated solution
function y=RMS(h)
num_of_pts=10/h;
t=linspace(0,10,num_of_pts);
y_lobatto=lobatto_solved_ode(t).';
y_ode45=ode45_solved_ode(t);
    nominator=norm(y_lobatto(:,1)-y_ode45(:,1),2);
    denominator=norm(y_ode45(:,1),2);
    y=nominator/denominator;
end

%Maximum error
%Ode45 solution as the accurate solution
%Lobatto IIID as the approximated solution
function y=MxError(h)
num_of_pts=10/h;
t=linspace(0,10,num_of_pts);
y_lobatto=lobatto_solved_ode(t).';
y_ode45=ode45_solved_ode(t);
    nominator=norm(y_lobatto(:,1)-y_ode45(:,1),inf);
    denominator=norm(y_ode45(:,1),inf);
    y=nominator/denominator;
end

%----Explicit Euler investigation---

%Root mean square error
%Ode45 solution as the accurate solution
%Explicit Euler method as the approximated solution
function y=RMS_euler(h)
num_of_pts=10/h;
t=linspace(0,10,num_of_pts);
y_euler=euler_solved_ode(t).';
y_ode45=ode45_solved_ode(t);
    nominator=norm(y_euler(:,1)-y_ode45(:,1),2);
    denominator=norm(y_ode45(:,1),2);
    y=nominator/denominator;
end

%Maximum error
%Ode45 solution as the accurate solution
%Explicit Euler method as the approximated solution
function y=MxError_euler(h)
num_of_pts=10/h;
t=linspace(0,10,num_of_pts);
y_euler=euler_solved_ode(t).';
y_ode45=ode45_solved_ode(t);
    nominator=norm(y_euler(:,1)-y_ode45(:,1),inf);
    denominator=norm(y_ode45(:,1),inf);
    y=nominator/denominator;
end
