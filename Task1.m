clear
close
second_order_ode();
manually_solved_ode();
function manually_solved_ode()
t=linspace(0,10,100);
 for i=1:(length(t))
   y(i)=2*exp(-t(i)/3)*sin(t(i));
end
plot(t,y, 'or');
end

function second_order_ode
    %SOLVE 9y''+6y'+10y=0
    %initial conditions: y(0)=0, y'(0)=2
    t= 0:0.01:10;
    
    initial_y = 0;
    initial_dydt = 2;
    
    [t,y]=ode45( @rhs, t, [initial_y, initial_dydt]);
    
    plot(t,y(:,1));
    hold on
    xlabel('t'); ylabel('y');
    
    function dydt=rhs(t,y)
        dydt_1 = y(2);
        dydt_2 = -(6*y(2) + 10*y(1))/9;
        
        dydt=[dydt_1; dydt_2];
    end
end