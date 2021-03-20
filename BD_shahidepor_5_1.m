%% Benders Decomposition, x as complicating variable
clc
clear all
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
    
    Objective=x(1)+3*x(2)+alpha;
    if iter>1
        MpCuts=[MpCuts, Obj+Lamda(:,1).*(ones(iter-1,1)*x(1)-X(:,1))+Lamda(:,2).*(ones(iter-1,1)*x(2)-X(:,2))<=ones(iter-1,1)*alpha];
    end
    MpConstraints=[MpConstraints,alpha>=-100,x(1)>=0,x(2)>=0,MpCuts];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(MpConstraints,Objective,ops);
    mpx=value(x);
    
    z_low=value(x(1)+3*x(2)+alpha);
    %% sub problem
    y=sdpvar(2,1,'full');
    v=sdpvar(2,1,'full');
    w=sdpvar(1,1,'full');
    
    Objective=y(1)+4*y(2)+20*(v(1)+v(2)+w);
    
    SpConstraints=[];
    SpConstraints=[SpConstraints, -2*x(1)-x(2)+y(1)-2*y(2)+v(1)+v(2)-w>=1];
    SpConstraints=[SpConstraints, 2*x(1)+2*x(2)-y(1)+3*y(2)+v(1)+v(2)-w>=1];
    SpConstraints=[SpConstraints, x(1)-mpx(1)==0,x(2)-mpx(2)==0];
%     SpConstraints=[SpConstraints, x-mpx==0];
    SpConstraints=[SpConstraints, y(1)>=0,y(2)>=0,v(1)>=0,v(2)>=0,w>=0];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(SpConstraints,Objective,ops)
    
    z_up=value(x(1)+3*x(2)+y(1)+4*y(2));
    eps=z_up-z_low;
    
    X(iter,:)=mpx;
    Obj(iter,1)=value(Objective);
%     Lamda(iter,1:2)=dual(SpConstraints(3));
    Lamda(iter,1)=dual(SpConstraints(3));
    Lamda(iter,2)=dual(SpConstraints(4));
    gap(iter)=eps
end
t2_end=toc(t2)
value(x)
value(y)

%% Benders Decomposition, y as complicating variable
% clc
clear all
%% master problem
clear all
t2=tic;
iter=0;
eps=1;
MpCuts=[];
MpConstraints=[];

while eps>0.01
    iter=iter+1;
    y=sdpvar(2,1,'full');
    alpha=sdpvar(1,1,'full');
    
    Objective=y(1)+4*y(2)+alpha;
    if iter>1
        MpCuts=[MpCuts, Obj+Lamda(:,1).*(ones(iter-1,1)*y(1)-X(:,1))+Lamda(:,2).*(ones(iter-1,1)*y(2)-X(:,2))<=ones(iter-1,1)*alpha];
    end
    MpConstraints=[MpConstraints,alpha>=-100,y(1)>=0,y(2)>=0,MpCuts];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(MpConstraints,Objective,ops);
    mpx=value(y);
    
    z_low=value(y(1)+4*y(2)+alpha);
    %% sub problem
    x=sdpvar(2,1,'full');
    v=sdpvar(2,1,'full');
    w=sdpvar(1,1,'full');
    
    Objective=x(1)+3*x(2)+20*(v(1)+v(2)+w);
    
    SpConstraints=[];
    SpConstraints=[SpConstraints, -2*x(1)-x(2)+y(1)-2*y(2)+v(1)+v(2)-w>=1];
    SpConstraints=[SpConstraints, 2*x(1)+2*x(2)-y(1)+3*y(2)+v(1)+v(2)-w>=1];
    SpConstraints=[SpConstraints, y(1)-mpx(1)==0,y(2)-mpx(2)==0];
    SpConstraints=[SpConstraints, x(1)>=0,x(2)>=0,v(1)>=0,v(2)>=0,w>=0];
    
    ops=sdpsettings('verbose',0,'solver','cplex');
    solution=optimize(SpConstraints,Objective,ops)
    
    z_up=value(x(1)+3*x(2)+20*(v(1)+v(2)+w)+y(1)+4*y(2));
    eps=z_up-z_low;
    
    X(iter,:)=mpx;
    Obj(iter,1)=value(Objective);
    Lamda(iter,1)=dual(SpConstraints(3));
    Lamda(iter,2)=dual(SpConstraints(4));
    gap(iter)=eps
end
t2_end=toc(t2)
value(x)
value(y)
