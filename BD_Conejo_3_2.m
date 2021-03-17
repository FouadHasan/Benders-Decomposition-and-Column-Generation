%% Benders Decomposition
clc
clear all
%% original Problem
t1=tic;
x=sdpvar(2,1,'full');
y=sdpvar(3,1,'full');

Objective=-2*y(1)-y(2)+y(3)+3*x(1)-3*x(2);

Constraints=[];
Constraints=[Constraints, y(1)+x(1)+x(2)<=3];
Constraints=[Constraints, 2*y(2)+3*x(1)<=12];
Constraints=[Constraints, y(3)-7*x(2)<=-16];
Constraints=[Constraints, -x(1)+x(2)<=2];
Constraints=[Constraints, x(1)>=0,x(2)>=0,y(1)>=0,y(2)>=0,y(3)>=0];

ops=sdpsettings('verbose',0,'solver','cplex');
solution=optimize(Constraints,Objective,ops);
value(x)
value(y)
t1_end=toc(t1)
%% master problem
clear all
t2=tic;
iter=0;
eps=1;
MpCuts=[];
MpConstraints=[];

while eps>0.01
    iter=iter+1;
    x=sdpvar(2,1,'full');
    alpha=sdpvar(1,1,'full');
    
    Objective=3*x(1)-3*x(2)+alpha;
    if iter>1
        MpCuts=[MpCuts, Obj+Lamda(:,1).*(ones(iter-1,1)*x(1)-X(:,1))+Lamda(:,2).*(ones(iter-1,1)*x(2)-X(:,2))<=ones(iter-1,1)*alpha];
    end
    MpConstraints=[MpConstraints,-x(1)+x(2)<=2,alpha>=-100,x(1)>=0,x(2)>=0,MpCuts];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(MpConstraints,Objective,ops);
    mpx=value(x);
    
    z_low=value(3*x(1)-3*x(2)+alpha);
    %% sub problem
    y=sdpvar(3,1,'full');
    v=sdpvar(3,1,'full');
    w=sdpvar(1,1,'full');
    
    Objective=-2*y(1)-y(2)+y(3)+20*(v(1)+v(2)+v(3)+w);
    
    SpConstraints=[];   
    SpConstraints=[SpConstraints, y(1)+x(1)+x(2)+v(1)+v(2)+v(3)-w<=3];
    SpConstraints=[SpConstraints, 2*y(2)+3*x(1)+v(1)+v(2)+v(3)-w<=12]; 
    SpConstraints=[SpConstraints, y(3)-7*x(2)+v(1)+v(2)+v(3)-w<=-16];    
    SpConstraints=[SpConstraints, x(1)-mpx(1)==0,x(2)-mpx(2)==0];
    SpConstraints=[SpConstraints, y(1)>=0,y(2)>=0,y(3)>=0,v(1)>=0,v(2)>=0,v(3)>=0,w>=0];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(SpConstraints,Objective,ops);
    
    z_up=value(-2*y(1)-y(2)+y(3)+3*x(1)-3*x(2));
    eps=z_up-z_low;
    
    X(iter,:)=mpx;
    Obj(iter,1)=value(Objective);
    Lamda(iter,1)=dual(SpConstraints(4));
    Lamda(iter,2)=dual(SpConstraints(5));
    gap(iter)=eps
end
t2_end=toc(t2)
value(x)
value(y)


