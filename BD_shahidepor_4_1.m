%% Benders Decomposition
clc
clear all

%% original Problem
x=sdpvar(1,1,'full');
y=sdpvar(1,1,'full');

Objective=x+y;

Constraints=[];
Constraints=[Constraints, 2*x+y>=3,x>=0,-5<=y<=4];
ops=sdpsettings('verbose',0,'solver','cplex');
solution=optimize(Constraints,Objective,ops);
value(x)
value(y)
%% taking x as complicating variable
% clc
clear all

iter=0;
eps=1;
MpCuts=[];
MpConstraints=[];
while eps>0.01
    iter=iter+1;
    x=sdpvar(1,1,'full');
    alpha=sdpvar(1,1,'full');
    %% master problem
    Objective=x+alpha;
    if iter>1
        MpCuts=[MpCuts, Obj'+Lamda'.*(ones(iter-1,1)*x-X')<=ones(iter-1,1)*alpha];
    end
    MpConstraints=[MpConstraints,x>=0,alpha>=-25,MpCuts];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(MpConstraints,Objective,ops);
    mpx=value(x);
    
    z_low=value(x+alpha);
    %% sub problem
    y=sdpvar(1,1,'full');
    
    Objective=y;
    SpConstraints=[];
    SpConstraints=[SpConstraints, 2*x+y>=3,-5<=y<=4,x-mpx==0];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(SpConstraints,Objective,ops);
    
    z_up=value(x+y);
    eps=z_up-z_low;
%     disp('value(x) value(y) value(alpha) dual(SpConstraints(3)) z_up z_low eps');
    s=[value(x) value(y) value(alpha) dual(SpConstraints(3)) z_up z_low eps]
    
    X(iter)=mpx;
    Obj(iter)=value(Objective);
    Lamda(iter)=dual(SpConstraints(3));
end

value(x)
value(y)


%% taking y as complicating variable
%% Benders Decomposition
clc
clear all

%% original Problem
x=sdpvar(1,1,'full');
y=sdpvar(1,1,'full');

Objective=x+y;

Constraints=[];
Constraints=[Constraints, 2*x+y>=3,x>=0,-5<=y<=4];
ops=sdpsettings('verbose',0,'solver','cplex');
solution=optimize(Constraints,Objective,ops);
value(x)
value(y)
%%
clc
clear all

iter=0;
eps=1;
MpCuts=[];

while eps>0.01
    iter=iter+1
    y=sdpvar(1,1,'full');
    alpha=sdpvar(1,1,'full');
    %% master problem
    Objective=y+alpha;
    MpConstraints=[];
    if iter>1
        MpCuts=[MpCuts, Obj'+Lamda'.*(ones(iter-1,1)*y-X')<=ones(iter-1,1)*alpha];
    end
    MpConstraints=[MpConstraints,-5<=y<=4,alpha>=-100000,MpCuts];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(MpConstraints,Objective,ops);
    mpx=value(y);
    
    z_low=value(y+alpha);
    %% sub problem
    x=sdpvar(1,1,'full');
    
    Objective=x;
    SpConstraints=[];
    SpConstraints=[SpConstraints, 2*x+y>=3,x>=0,y-mpx==0];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(SpConstraints,Objective,ops);
    
    z_up=value(x+y);
    eps=z_up-z_low;
    disp('value(x) value(y) value(alpha) dual(SpConstraints(3)) z_up z_low eps')
    s=[value(x) value(y) value(alpha) dual(SpConstraints(3)) z_up z_low eps]
    
    X(iter)=mpx;
    Obj(iter)=value(Objective);
    Lamda(iter)=dual(SpConstraints(3));
end

value(x)
value(y)
