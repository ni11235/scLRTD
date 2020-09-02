%%% First we need install the package "tensor_toolbox"[2]
%%[2]: B. W. Bader, T. G. Kolda, et al., Matlab tensor toolbox version 3.0-dev, Available online (Oct. 2017)

% data input and tensor represent
% A is the dataset 3
% represent A as a third-order tensor Y
A=importdata("Raw_names.txt");
A=A.data;
Y=tenzeros(size(A,1),fix(size(A,2)/2),2);
Y_missing=tenzeros(size(A,1),fix(size(A,2)/2),2);
Y(:,:,1)=A(:,1:fix(size(A,2)/2));
Y(:,:,2)=A(:,fix(size(A,2)/2)+1:2*fix(size(A,2)/2));
Y_missing(:,:,1)=A(:,1:fix(size(A,2)/2));
Y_missing(:,:,2)=A(:,fix(size(A,2)/2)+1:2*fix(size(A,2)/2));
%% W is the set that data can be observed
W=(Y_missing-Y);W=W.data;
W(W==0)=1;W(W<=0)=0;
%%
T=Y.*W;

% R is the cp rank, beta and inf_L are paremeters
R=2;
beta=8066.0; %If you recalculate, please cancel the comment below
%%%% the bata is calculated by:

%     for i=1:inf_L
%         P=cp_nmu(T+(1-W).*Xold,R);
%         L{i}=P;
%     end
%     Plambda=zeros(1,inf_L);
%     for i=1:inf_L
%         Nuclearnorm_P=sum(abs(L{i}.lambda));
%         Plambda(i)=Nuclearnorm_P;
%     end
%     [value,index]=min(Plambda);
%     P_need=L{index};

%%%% The beta is calculated by following, please refer Lines 208-209 on Page 4   
   % beta=0.1P_need.lambda{R}, 0.01P_need.lambda{R}, 0.001P_need.lambda{R}
   
inf_L=3;eps=1e-5;
L=cell(1,inf_L);

%%%% Low rank decomposition and get the imputed tensor Xneed
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
A_imputed=zeros(size(A_true));
Xneed=(T+(1-W).*Xopt); %%% Xneed is the imputed tensor

A_imputed(:,1:fix(size(A_true,2)/2))=Xneed(:,:,1);
A_imputed(:,fix(size(A_true,2)/2)+1:2*fix(size(A_true,2)/2))=Xneed(:,:,2);

Xi=A_imputed;%%% get the imputed dataset 3
 
