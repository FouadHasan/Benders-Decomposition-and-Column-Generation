clc
clear all
x=sdpvar(2,1,'full');

Objective=2*x(1)+x(2);
Constraints=[];
Constraints=[Constraints,x(1)<=5,x(2)<=5,x(1)+x(2)<=9,x(1)-x(2)<=5,-x(1)-x(2)<=-2,-3*x(1)-x(2)<=-3];
Constraints=[Constraints, x>=0];

ops=sdpsettings('verbose',0','solver','cplex');
solution=optimize(Constraints, Objective, ops);
value(x)
%% Column Generation
clearvars x
x=sdpvar(2,1,'full');

% initialization
c=[-1 -1];
Objective=c(1)*x(1)+c(2)*x(2);
Constraints=[];
Constraints=[Constraints,x(1)<=5,x(2)<=5,x>=0];
solution=optimize(Constraints, Objective, ops);

z(1)=value(2*x(1)+x(2));

r(1,1)=value(x(1)+x(2));
r(1,2)=value(x(1)-x(2));
r(1,3)=value(-x(1)-x(2));
r(1,4)=value(-3*x(1)-x(2));

X(:,1)=value(x);

clearvars x
x=sdpvar(2,1,'full');

c=[-2 1];
Objective=c(1)*x(1)+c(2)*x(2);
Constraints=[];
Constraints=[Constraints,x(1)<=5,x(2)<=5,x>=0];
solution=optimize(Constraints, Objective, ops);

z(2)=value(2*x(1)+x(2));

r(2,1)=value(x(1)+x(2));
r(2,2)=value(x(1)-x(2));
r(2,3)=value(-x(1)-x(2));
r(2,4)=value(-3*x(1)-x(2));

X(:,2)=value(x);
%%
iter=0;
v=0;sigma=1;

while v<sigma
    iter=iter+1;
    % master problem
    clearvars u
    u=sdpvar(iter+1,1,'full');
    Objective=sum(u.*z');
    Constraints=[];
    Constraints=[Constraints,sum(u.*r(:,1))<=9,sum(u.*r(:,2))<=4,sum(u.*r(:,3))<=-2,sum(u.*r(:,4))<=-3,sum(u)==1,u>=0];
    solution=optimize(Constraints, Objective, ops)
    lamda(1)=dual(Constraints(1));
    lamda(2)=dual(Constraints(2));
    lamda(3)=dual(Constraints(3));
    lamda(4)=dual(Constraints(4));
    sigma=dual(Constraints(5));
    % relaxed subproblem
    clearvars c x
    x=sdpvar(2,1,'full');
    c(1)=2-(lamda(1)*1+lamda(2)*1+lamda(3)*(-1)+lamda(4)*(-3));
    c(2)=1-(lamda(1)*1+lamda(2)*(-1)+lamda(3)*(-1)+lamda(4)*(-1));
    
    Objective=c(1)*x(1)+c(2)*x(2);
    Constraints=[];
    Constraints=[Constraints,x(1)<=5,x(2)<=5,x>=0];
    ops=sdpsettings('verbose',0','solver','cplex');
    solution=optimize(Constraints, Objective, ops);
    
    %convergence checking condiiton
    v=value(c(1)*x(1)+c(2)*x(2));
    
    z(iter+2)=value(2*x(1)+x(2));
    r(iter+2,1)=value(x(1)+x(2));
    r(iter+2,2)=value(x(1)-x(2));
    r(iter+2,3)=value(-x(1)-x(2));
    r(iter+2,4)=value(-3*x(1)-x(2));
    X(:,iter+2)=value(x)
end
opt_x=value(u(1)*X(:,1)+u(2)*X(:,2)+u(3)*X(:,3)+u(4)*X(:,4))




