% input scores file
MScore=dlmread('networkfiles/emt26network.phs');
Sc=zeros(26,1);
for i=1:length(MScore(:,1))
    Sc(MScore(i,1)+1)=1;
    Sc(MScore(i,2)+1)=-1;
end
tries=100;
% input Clamp1 results
Clamp1=dlmread('Clamp1.dat');
Clamp2d=zeros(26);
com0='/Users/kessler/Work/boolean/TRANSIT/TRANS10';
for i=1:25
    for j=i+1:26
        [st,cout]=system([com0,' "emt" 26 3 400 ',num2str(tries),' ',...
            num2str(i*Clamp1(i,2)),' ',num2str(j*Clamp1(j,2))]);
        Clamp2d(i,j)=str2double(cout);
        disp([i,j,Clamp2d(i,j)])
    end
end
