%data missing and tensor represent
%A_true=ones(4,7);A_missing=ones(4,7);A_missing(:,1:3)=0;
A_true=Full3_data;A_missing=Raw3_data;
Y=tenzeros(size(A_true,1),fix(size(A_true,2)/2),2);
Y_missing=tenzeros(size(A_missing,1),fix(size(A_missing,2)/2),2);
Y(:,:,1)=A_true(:,1:fix(size(A_true,2)/2));
Y(:,:,2)=A_true(:,fix(size(A_true,2)/2)+1:2*fix(size(A_true,2)/2));
Y_missing(:,:,1)=A_missing(:,1:fix(size(A_missing,2)/2));
Y_missing(:,:,2)=A_missing(:,fix(size(A_missing,2)/2)+1:2*fix(size(A_missing,2)/2));
W=(Y_missing-Y);W=W.data;
W(W==0)=1;W(W<=0)=0;
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
% MSE comparison
A_imputed=zeros(size(A_true));
Xneed=(T+(1-W).*Xopt);

A_imputed(:,1:fix(size(A_true,2)/2))=Xneed(:,:,1);
A_imputed(:,fix(size(A_true,2)/2)+1:2*fix(size(A_true,2)/2))=Xneed(:,:,2);
%A_imputed=max(A_imputed,0);
%A_imputed=abs(A_imputed);
Xr = Full3_data;lXr = log10(Xr+1);
Xi=A_imputed; lXi = log10(Xi+1);
MSE_record= norm(lXr-lXi,'fro');

 






