clc
clear all

%% Benders Decomposition
clear all

%% master problem(minimization)
iter=0;
eps=1;
MpCuts=[];
MpConstraints=[];
while eps>0.01
    iter=iter+1
    x=sdpvar(1,1,'full');
    alpha=sdpvar(1,1,'full');
    
    Objective=-x/4+alpha;
    if iter>1
        MpCuts=[MpCuts, Obj'+Lamda'.*(ones(iter-1,1)*x-X')<=ones(iter-1,1)*alpha];
    end
    MpConstraints=[MpConstraints,0<=x<=16,alpha>=-25,MpCuts];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(MpConstraints,Objective,ops);
    mpx=value(x);
    
    z_low=value(-x/4+alpha);
    %% sub problem (maximization)
    y=sdpvar(1,1,'full');
    Objective=0;
    Objective=-y;
    SpConstraints=[];
    SpConstraints=[SpConstraints, y-x<=5,y-.5*x<=15/2,y+.5*x<=35/2,-y+x<=10,y>=0,x-mpx==0];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(SpConstraints,Objective,ops);
    
    z_up=value(-y-x/4);
    eps=z_up-z_low;
    disp('value(x) value(y) value(alpha) dual(SpConstraints(6)) z_up z_low eps')
    s=[value(x) value(y) value(alpha) dual(SpConstraints(6)) z_up z_low eps]
    
    X(iter)=mpx;
    Obj(iter)=value(Objective);
    Lamda(iter)=dual(SpConstraints(6));
end

value(x)
value(y)


