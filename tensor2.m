

T=Y.*W;

%R is the cp rank，beta and inf_L are paremeters
inf_L=3;
eps=1e-5;
beta=90.533;
R=2;L=cell(1,inf_L);

Xold=tensor(zeros(T.size));
for j=1:100
    for i=1:inf_L
        P=cp_nmu(T+(1-W).*Xold,R);
        L{i}=P;
    end
    Plambda=zeros(1,inf_L);
    for i=1:inf_L
        Nuclearnorm_P=sum(abs(L{i}.lambda));
        Plambda(i)=Nuclearnorm_P;
    end
    [value,index]=min(Plambda);
    P_need=L{index};
    P_new_lambda=max(P_need.lambda-beta,0);
    P_new_U=P_need.U;
    P_new=ktensor(P_new_lambda,P_new_U);
    Xnew=P_new; Xnew=tensor( Xnew);
Fold=(1/2)*(norm(T-(1-W).*Xold))^2+beta*value;    
Fnew=(1/2)*(norm(T-(1-W).*Xnew))^2+beta*value;
    if(abs(Fnew-Fold)/Fold<eps)
        Xopt=Xnew;
        break;
    end
Xold=Xnew;
end
Xopt=Xnew;

Xneed=(T+(1-W).*Xopt);



 






