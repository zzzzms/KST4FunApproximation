function [ H ] = Smoothness2D( V, T, d, r )
%SMOOTHNESS2D Summary of this function goes here
%   Detailed explanation goes here
% Written by Yidong Xu.

global Factorials;

n=(d+2)*(d+1)/2;
m=size(T,1);

TR=triangulation(T,V);

E=edges(TR);
Attachments=edgeAttachments(TR,E);

% EdgeCandidates:
% 1st and 2nd components: index of the edge vertex
% 3rd and 4th components: index of the triangle attached to the edge 
EdgeCandidates=zeros(size(E,1),4);
num_of_EdgeCandidates=0;

for i=1:size(Attachments,1)
    if length(Attachments{i,1})>1
        num_of_EdgeCandidates=num_of_EdgeCandidates+1;
        EdgeCandidates(num_of_EdgeCandidates,:)=[E(i,:),Attachments{i,1}];
    end
end

EdgeCandidates(num_of_EdgeCandidates+1:end,:)=[];
% EdgeCandidates=EdgeCandidates(1:num_of_EdgeCandidates,:);

% num_of_H_rows=num_of_EdgeCandidates*(d+1-r/2)*(r+1);

% if num_of_H_rows==0
%     H=[];
% else
%     H=spalloc(num_of_H_rows,n*m,num_of_H_rows*(n+1));
% end

% H=spalloc(num_of_H_rows,n*m,num_of_H_rows*(n+1));



% preallocating will keep getting larger.
preallocating=1000;
num_of_H_nonzero_elements=0;
row_index=zeros(preallocating,1);
% row_index=zeros((num_of_EdgeCandidates*(d+1-r/2)*(r+1))*((r+2)*(r+1)/2),1);
column_index=row_index;
H_elements=row_index;
row=0;

for e=1:num_of_EdgeCandidates
    
    e1=EdgeCandidates(e,1);
    e2=EdgeCandidates(e,2);
    T1=EdgeCandidates(e,3);
    T2=EdgeCandidates(e,4);
    
    if T(T1,1)~=e1 && T(T1,1)~=e2
        permutation1=[1,2,3];
        v1=T(T1,1);
        v2=T(T1,2);
        v3=T(T1,3);
    elseif T(T1,2)~=e1 && T(T1,2)~=e2
        permutation1=[3,1,2];
        v1=T(T1,2);
        v2=T(T1,3);
        v3=T(T1,1);
    else
        permutation1=[2,3,1];
        v1=T(T1,3);
        v2=T(T1,1);
        v3=T(T1,2);
    end
    
    if T(T2,1)~=e1 && T(T2,1)~=e2
        permutation2=[1,2,3];
        v4=T(T2,1);
    elseif T(T2,2)~=e1 && T(T2,2)~=e2
        permutation2=[3,1,2];
        v4=T(T2,2);
    else
        permutation2=[2,3,1];
        v4=T(T2,3);
    end
    
    v4bary=[1,1,1;V([v1,v2,v3],1)';V([v1,v2,v3],2)']\[1;V(v4,:)'];
    
    % Note that v4barypowers(1,1)=b1^0=1, v4barypowers(1,2)=b1^1=b1, and so on.
    v4barypowers=ones(3,r+1);
    if r>=1
        v4barypowers(:,2)=v4bary;
    end
    for j=3:r+1
        v4barypowers(:,j)=v4barypowers(:,j-1).*v4bary;
    end
    
%     b1=v4bary(1);
%     b2=v4bary(2);
%     b3=v4bary(3);
    % Note that b1powers(1)=b1^0=1, b1powers(2)=b1^1=b1, and so on.
    % Similarly for b2powers and b3powers.
%     b1powers=ones(d+1,1);
%     b2powers=ones(d+1,1);
%     b3powers=ones(d+1,1);
%     b1powers(1)=1;
%     b2powers(1)=1;
%     b3powers(1)=1;
%     for j=2:d+1
%         b1powers(j)=b1powers(j-1)*b1;
%         b2powers(j)=b2powers(j-1)*b2;
%         b3powers(j)=b3powers(j-1)*b3;
%     end
    
    for s=0:r
        for j=0:d-s
            k=d-s-j;
            row=row+1;
            for nu=s:-1:0
                for mu=s-nu:-1:0
                    kappa=s-nu-mu;
                    coef=Factorials(s+1)/Factorials(nu+1)/Factorials(mu+1)/Factorials(kappa+1)...
                        *v4barypowers(1,nu+1)*v4barypowers(2,mu+1)*v4barypowers(3,kappa+1);
                    index=[nu,k+mu,j+kappa];
                    index=index(permutation1);
                    pos=(1+d-index(1))*(d-index(1))/2+index(3)+1;
                    
                    if num_of_H_nonzero_elements+1>length(H_elements)
                        H_elements=[H_elements;zeros(1+preallocating,1)];
                        row_index=[row_index;zeros(1+preallocating,1)];
                        column_index=[column_index;zeros(1+preallocating,1)];
                        preallocating=preallocating*2;
                    end
                    row_index(num_of_H_nonzero_elements+1)=row;
                    column_index(num_of_H_nonzero_elements+1)=(T1-1)*n+pos;
                    H_elements(num_of_H_nonzero_elements+1)=coef;
                    num_of_H_nonzero_elements=num_of_H_nonzero_elements+1;
                    
%                     num_of_H_nonzero_elements=num_of_H_nonzero_elements+1;
%                     row_index(num_of_H_nonzero_elements)=row;
%                     column_index(num_of_H_nonzero_elements)=(T1-1)*n+pos;
%                     H_elements(num_of_H_nonzero_elements)=coef;
%                     H(row,(T1-1)*n+pos)=coef;
                end
            end
            
            index=[s,j,k];
            index=index(permutation2);
            pos=(1+d-index(1))*(d-index(1))/2+index(3)+1;
            if num_of_H_nonzero_elements+1>length(H_elements)
                H_elements=[H_elements;zeros(1+preallocating,1)];
                row_index=[row_index;zeros(1+preallocating,1)];
                column_index=[column_index;zeros(1+preallocating,1)];
                preallocating=preallocating*2;
            end
            row_index(num_of_H_nonzero_elements+1)=row;
            column_index(num_of_H_nonzero_elements+1)=(T2-1)*n+pos;
            H_elements(num_of_H_nonzero_elements+1)=-1;
            num_of_H_nonzero_elements=num_of_H_nonzero_elements+1;
            
%             num_of_H_nonzero_elements=num_of_H_nonzero_elements+1;
%             row_index(num_of_H_nonzero_elements)=row;
%             column_index(num_of_H_nonzero_elements)=(T2-1)*n+pos;
%             H_elements(num_of_H_nonzero_elements)=-1;
%             H(row,(T2-1)*n+pos)=-1;
            
        end
        
    end
    
end

% H=H(1:row,:);

row_index(num_of_H_nonzero_elements+1:end)=[];
column_index(num_of_H_nonzero_elements+1:end)=[];
H_elements(num_of_H_nonzero_elements+1:end)=[];
H=sparse(row_index,column_index,H_elements,row,n*m);

% index=1:num_of_H_nonzero_elements;
% row_index=row_index(index);
% column_index=column_index(index);
% H_elements=H_elements(index);
% H=sparse(row_index,column_index,H_elements,row,n*m);




end

