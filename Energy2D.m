function [ K ] = Energy2D( Analyze, d )
%ENERGY2D Summary of this function goes here
%   Detailed explanation goes here
% Written by Yidong Xu.

global C1Index;
global C2Index;
global C3Index;
global FindIJK;
global G;

n=(d+2)*(d+1)/2;
m=size(Analyze,1);

if d<2
    K=sparse(n*m,n*m);
    return;
end

% direction_bary contains powers of x_direction_bary and y_direction_bary.
% direction_bary(:,:,1) is 0th power, direction_bary(:,:,2) is 1st power,
% and so on.
direction_bary=ones(m,6);
% direction_bary(:,:,2)=[Analyze(:,[6,7,8]),Analyze(:,[9,10,11])];
direction_bary(:,:,2)=Analyze(:,[6,7,8,9,10,11]);
direction_bary(:,:,3)=direction_bary(:,:,2).^2;

% x_direction_bary=Analyze(:,[6,7,8]);
% y_direction_bary=Analyze(:,[9,10,11]);
% x_direction_bary1powers=ones(m,3);
% x_direction_bary2powers=ones(m,3);
% x_direction_bary3powers=ones(m,3);
% y_direction_bary1powers=ones(m,3);
% y_direction_bary2powers=ones(m,3);
% y_direction_bary3powers=ones(m,3);

% for i=2:3
%     x_direction_bary1powers(:,i)=x_direction_bary1powers(:,i-1).*Analyze(:,6);
%     x_direction_bary2powers(:,i)=x_direction_bary2powers(:,i-1).*Analyze(:,7);
%     x_direction_bary3powers(:,i)=x_direction_bary3powers(:,i-1).*Analyze(:,8);
%     y_direction_bary1powers(:,i)=y_direction_bary1powers(:,i-1).*Analyze(:,9);
%     y_direction_bary2powers(:,i)=y_direction_bary2powers(:,i-1).*Analyze(:,10);
%     y_direction_bary3powers(:,i)=y_direction_bary3powers(:,i-1).*Analyze(:,11);
% end


d1=d-1;
d2=d-2;
% n1=(d+1)*d/2;
n2=d*d1/2;
DerivativeInfluence=false(n2,n);

for i=1:n
    coef=false(n,1);
    coef(i,1)=true;
    coef=coef(C1Index{d,1}) | coef(C2Index{d,1}) | coef(C3Index{d,1});
    coef=coef(C1Index{d1,1}) | coef(C2Index{d1,1}) | coef(C3Index{d1,1});
    DerivativeInfluence(:,i)=coef;
end



DxDxDerivativeInfluence=zeros(n2,m);
DyDyDerivativeInfluence=zeros(n2,m);
DxDyDerivativeInfluence=zeros(n2,m);

% DxDxDerivativeInfluence=cell(n,1);
% DyDyDerivativeInfluence=cell(n,1);
% DxDyDerivativeInfluence=cell(n,1);
temp=d*d1;

for i=1:n
    ijk1=FindIJK{d,1}(i,:);
    if i>1
        DxDxDerivativeInfluence(:,:,i)=0;
        DxDyDerivativeInfluence(:,:,i)=0;
        DyDyDerivativeInfluence(:,:,i)=0;
    end
    
    FindIndex=find(DerivativeInfluence(:,i));
    for j=1:length(FindIndex)
        pos=FindIndex(j);
        if d==2
            ijk2=[0,0,0];
        else
            ijk2=FindIJK{d2,1}(pos,:);
        end
        
        nu=ijk1(1)-ijk2(1);
        mu=ijk1(2)-ijk2(2);
        kappa=ijk1(3)-ijk2(3);
        
        if max([nu,mu,kappa])==2
            index=find([nu,mu,kappa]==2,1);
            DxDxDerivativeInfluence(pos,:,i)=temp*direction_bary(:,index,3);
            DyDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,index+3,3);
            DxDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,index,2).*direction_bary(:,index+3,2);
%             if index==1
%                 DxDxDerivativeInfluence(pos,:,i)=temp*direction_bary(:,1,3);
%                 DyDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,4,3);
%                 DxDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,1,2).*direction_bary(:,4,2);
%             elseif index==2
%                 DxDxDerivativeInfluence(pos,:,i)=temp*direction_bary(:,2,3);
%                 DyDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,5,3);
%                 DxDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,2,2).*direction_bary(:,5,2);
%             else
%                 DxDxDerivativeInfluence(pos,:,i)=temp*direction_bary(:,3,3);
%                 DyDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,6,3);
%                 DxDyDerivativeInfluence(pos,:,i)=temp*direction_bary(:,3,2).*direction_bary(:,6,2);
%             end
            
            
        else
            index=find([nu,mu,kappa]==1);
            DxDxDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,index(1),2).*direction_bary(:,index(2),2);
            DyDyDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,index(1)+3,2).*direction_bary(:,index(2)+3,2);
            DxDyDerivativeInfluence(pos,:,i)=temp*(direction_bary(:,index(1),2).*direction_bary(:,index(2)+3,2)...
                +direction_bary(:,index(2),2).*direction_bary(:,index(1)+3,2));
            
            
%             if nu==1
%                 if mu==1
%                     DxDxDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,1,2).*direction_bary(:,2,2);
%                     DyDyDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,4,2).*direction_bary(:,5,2);
%                     DxDyDerivativeInfluence(pos,:,i)=temp*(direction_bary(:,1,2).*direction_bary(:,5,2)...
%                         +direction_bary(:,4,2).*direction_bary(:,2,2));
%                 else
%                     DxDxDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,1,2).*direction_bary(:,3,2);
%                     DyDyDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,4,2).*direction_bary(:,6,2);
%                     DxDyDerivativeInfluence(pos,:,i)=temp*(direction_bary(:,1,2).*direction_bary(:,6,2)...
%                         +direction_bary(:,4,2).*direction_bary(:,3,2));
%                 end
%                 
%             else
%                 DxDxDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,2,2).*direction_bary(:,3,2);
%                 DyDyDerivativeInfluence(pos,:,i)=temp*2*direction_bary(:,5,2).*direction_bary(:,6,2);
%                 DxDyDerivativeInfluence(pos,:,i)=temp*(direction_bary(:,2,2).*direction_bary(:,6,2)...
%                     +direction_bary(:,5,2).*direction_bary(:,3,2));
%             end
            
        end
        
%         DxDxDerivativeInfluence{i,1}(pos,:)=temp*2/Factorials(nu)/Factorials(mu)/Factorials(kappa)...
%             *x_direction_bary1powers(:,nu).*x_direction_bary2powers(:,mu).*x_direction_bary3powers(:,kappa);
%         DyDyDerivativeInfluence{i,1}(pos,:)=temp*2/Factorials(nu)/Factorials(mu)/Factorials(kappa)...
%             *y_direction_bary1powers(:,nu).*y_direction_bary2powers(:,mu).*y_direction_bary3powers(:,kappa);
        
    end
    
end



row_index=reshape(1:n2*m,[],m);
column_index=repmat(1:n:1+n*(m-1),n2,1);

for i=2:n
    row_index(:,:,i)=row_index(:,:,i-1);
    column_index(:,:,i)=column_index(:,:,i-1)+1;
end

DxDx=sparse(row_index(:),column_index(:),DxDxDerivativeInfluence(:),n2*m,n*m);
DyDy=sparse(row_index(:),column_index(:),DyDyDerivativeInfluence(:),n2*m,n*m);
DxDy=sparse(row_index(:),column_index(:),DxDyDerivativeInfluence(:),n2*m,n*m);

if d==2
    G_Matrix=1;
else
    G_Matrix=G{d2,d2};
end
GG=kron(spdiags(Analyze(:,5)/((d1+d1)*(d1+d1-1)/2*nchoosek(d2+d2,d2)),0,m,m),G_Matrix);

K=DxDx'*GG*DxDx+DyDy'*GG*DyDy+2*DxDy'*GG*DxDy;



end
