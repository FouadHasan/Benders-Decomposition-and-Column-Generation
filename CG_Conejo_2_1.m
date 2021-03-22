clc
clear all
x=sdpvar(3,1,'full');

Objective=-4*x(1)-x(2)-6*x(3);
Constraints=[];
Constraints=[Constraints,-x(1)<=-1,x(1)<=2,-x(2)<=-1,x(2)<=2,-x(3)<=-1,x(3)<=2,3*x(1)+2*x(2)+4*x(3)<=17];
Constraints=[Constraints, x>=0];

ops=sdpsettings('verbose',0','solver','cplex');
solution=optimize(Constraints, Objective, ops)


%% Column Generation

clearvars x
x=sdpvar(3,1,'full');

% initialization
c=[-1 -1 -1];
Objective=c(1)*x(1);
Constraints=[];
Constraints=[Constraints,-x(1)<=-1,x(1)<=2];
ops=sdpsettings('verbose',0','solver','cplex');
solution=optimize(Constraints, Objective, ops);

Objective=c(2)*x(2);
Constraints=[];
Constraints=[Constraints,-x(2)<=-1,x(2)<=2];
solution=optimize(Constraints, Objective, ops);

Objective=c(3)*x(3);
Constraints=[];
Constraints=[Constraints,-x(3)<=-1,x(3)<=2];
solution=optimize(Constraints, Objective, ops);

z(1)=value(-4*x(1)-x(2)-6*x(3));
r(1)=value(3*x(1)+2*x(2)+4*x(3));
X(:,1)=value(x);
clearvars c x
c=[1 1 -1];
x=sdpvar(3,1,'full');
Objective=c(1)*x(1);
Constraints=[];
Constraints=[Constraints,-x(1)<=-1,x(1)<=2];
solution=optimize(Constraints, Objective, ops);

Objective=c(2)*x(2);
Constraints=[];
Constraints=[Constraints,-x(2)<=-1,x(2)<=2];
solution=optimize(Constraints, Objective, ops);

Objective=c(3)*x(3);
Constraints=[];
Constraints=[Constraints,-x(3)<=-1,x(3)<=2];
solution=optimize(Constraints, Objective, ops);

z(2)=value(-4*x(1)-x(2)-6*x(3));
r(2)=value(3*x(1)+2*x(2)+4*x(3));
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
    Constraints=[Constraints,sum(u.*r')<=17,sum(u)==1,u>=0];
    solution=optimize(Constraints, Objective, ops);
    lamda=dual(Constraints(1));
    sigma=dual(Constraints(2));
    % relaxed subproblem
    clearvars c x
    x=sdpvar(3,1,'full');
    c(1)=-4-lamda*3;
    c(2)=-1-lamda*2;
    c(3)=-6-lamda*4;
    
    Objective=c(1)*x(1);
    Constraints=[];
    Constraints=[Constraints,-x(1)<=-1,x(1)<=2];
    ops=sdpsettings('verbose',0','solver','cplex');
    solution=optimize(Constraints, Objective, ops);
    
    Objective=c(2)*x(2);
    Constraints=[];
    Constraints=[Constraints,-x(2)<=-1,x(2)<=2];
    solution=optimize(Constraints, Objective, ops);
    
    Objective=c(3)*x(3);
    Constraints=[];
    Constraints=[Constraints,-x(3)<=-1,x(3)<=2];
    solution=optimize(Constraints, Objective, ops);
    
    %convergence checking condiiton
    v=value(c(1)*x(1)+c(2)*x(2)+c(3)*x(3));
    
    z(iter+2)=value(-4*x(1)-x(2)-6*x(3));
    r(iter+2)=value(3*x(1)+2*x(2)+4*x(3));
    X(:,iter+2)=value(x);
end
opt_x=value(u(1)*X(:,1)+u(2)*X(:,2)+u(3)*X(:,3))