function [  ] = Initialization2D( d )
%INITIALIZATION2D Summary of this function goes here
%   Detailed explanation goes here
% Written by Yidong Xu.

global G;

G=cell(d,d);

for i=1:d
    for j=i:d
        n1=(i+2)*(i+1)/2;
        n2=(j+2)*(j+1)/2;
        G{i,j}=zeros(n1,n2);
        
        index1=0;
        for i1=i:-1:0
            for j1=i-i1:-1:0
                k1=i-i1-j1;
                index1=index1+1;
                
                index2=0;
                for i2=j:-1:0
                    for j2=j-i2:-1:0
                        k2=j-i2-j2;
                        index2=index2+1;
                        G{i,j}(index1,index2)=nchoosek(i1+i2,i1)*nchoosek(j1+j2,j1)*nchoosek(k1+k2,k1);
                    end
                end
            end
        end
        
    end
    
end

for i=1:d
    for j=1:i-1
        G{i,j}=G{j,i}';
    end
end

% n1=(d1+2)*(d1+1)/2;
% n2=(d2+2)*(d2+1)/2;
% 
% G=zeros(n1,n2);
% index1=0;
% 
% for i1=n1:-1:0
%     for j1=n1-i1:-1:0
%         k1=n1-i1-j1;
%         index1=index1+1;
%         
%         index2=0;
%         for i2=n2:-1:0
%             for j2=n2-i2:-1:0
%                 k2=n2-i2-j2;
%                 index2=index2+1;
%                 G(index1,index2)=nchoosek(i1+i2,i1)*nchoosek(j1+j2,j1)*nchoosek(k1+k2,k1);
%             end
%         end
%         
%     end
% end











global C1Index;%logical value
global C2Index;
global C3Index;

C1Index=cell(d,1);
C2Index=cell(d,1);
C3Index=cell(d,1);

for deg=1:d
    C1Index{deg,1}=false((deg+2)*(deg+1)/2,1);
    C1Index{deg,1}(1:(deg+1)*deg/2)=true;
    
    C2Index{deg,1}=false((deg+2)*(deg+1)/2,1);
    C3Index{deg,1}=false((deg+2)*(deg+1)/2,1);
    
    index=0;
    for i=deg:-1:1
        for j=deg-i:-1:0
            index=index+1;
            C2Index{deg,1}(index+deg-i+1)=true;
            C3Index{deg,1}(index+deg-i+2)=true;
        end
    end
    
%     disp([C1Index{deg,1},C2Index{deg,1},C3Index{deg,1}]);
    
end








% n=(d+2)*(d+1)/2;
% C1Index=false(n,d);
% C2Index=false(n,d);
% C3Index=false(n,d);
% 
% for i=1:d
%     C1Index(1:(i+1)*i/2,i)=true;
% end
% 
% for deg=1:d
%     index=0;
%     for i=deg:-1:1
%         for j=deg-i:-1:0
%             index=index+1;
%             C2Index(index+deg-i+1,deg)=true;
%             C3Index(index+deg-i+2,deg)=true;
%         end
%     end
% end









global Factorials;
% Note that Factorials(0)=0!=1 Factorials(n)=(n-1)!
Factorials=ones(d+1,1);

for i=3:d+1
    Factorials(i)=Factorials(i-1)*(i-1);
end











% Given the index, find i, j, and k.
global FindIJK;
FindIJK=cell(d,1);
% FindIJK{1,1}=zeros(1,3);
for deg=1:d
    FindIJK{deg,1}=zeros((deg+2)*(deg+1)/2,3);
    index=0;
    for i=deg:-1:0
        for j=deg-i:-1:0
            k=deg-i-j;
            index=index+1;
            FindIJK{deg,1}(index,:)=[i,j,k];
        end
    end
    
end



end

