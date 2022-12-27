function nI=LaisCrossApproximation(nPhi)
%This is a code written by Dr. Ming-Jun Lai on Nov. 10
%2022 based on some codes jointly developed by Dr. Ming-Jun Lai
% and his Ph.D. student Kenneth Allen in 2020 who 
%received a Ph.D. in Dec. 2021. 
[n1,r]=size(nPhi); 
I_initial = randperm(n1,r);
[Im,km] = maxvol(nPhi,I_initial);
[Igs,kgs] = simple_greedy_maxvol(nPhi,I_initial);
[Ig,kg] = greedy_maxvol(nPhi,I_initial);
%nI=[Im;Igs]; nI=unique(nI);
[a,m]=max([abs(det(nPhi(Im,:)))
 abs(det(nPhi(Igs,:)))
abs(det(nPhi(Ig,:)))]); 
if m==1
nI=Im;
else if m==2
        nI=Igs;
else
    nI=Ig;
end
end
end