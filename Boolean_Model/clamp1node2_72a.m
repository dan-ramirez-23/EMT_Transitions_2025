% Clamp 4th node ON and a 2nd node
% input scores file
MScore=dlmread('networkfiles/emt72network.phs');
Sc=zeros(72,1);
for i=1:length(MScore(:,1))
    Sc(MScore(i,1)+1)=1;
    Sc(MScore(i,2)+1)=-1;
end
tries=20;
Vp=find(Sc==1);
SuccessP=zeros(size(Vp));
com0='/Users/kessler/Work/boolean/TRANSIT/TRANS10 "emt" 72 1 8 ';
parfor i=1:length(SuccessP)
    [st,cout]=system([com0,num2str(tries),' 4 ',num2str(-(Vp(i)))]);
    SuccessP(i)=str2double(cout);
end
disp(SuccessP)
Vn=find(Sc==-1);
SuccessN=zeros(size(Vn));
parfor i=1:length(SuccessN)
    [~,cout]=system([com0,num2str(tries),' 4 ',num2str(Vn(i))]);
    SuccessN(i)=str2double(cout);
end
disp(SuccessN)
V0=find(Sc==0);
Success0p=zeros(length(V0),1);
Success0m=zeros(length(V0),1);
parfor i=1:length(V0)
    [~,cout]=system([com0,num2str(tries),' 2 ',num2str(-(V0(i)))]);
    Success0p(i)=str2double(cout);
    [~,cout]=system([com0,num2str(tries),' 2 ',num2str(V0(i))]);
    Success0m(i)=str2double(cout);
end
Success0=[Success0p,Success0m];
disp(Success0);

All=zeros(72,2);
for i=1:length(Vp)
    ndx=Vp(i);
    All(ndx,1)=SuccessP(i);
    All(ndx,2)=-1;
end
for i=1:length(Vn)
    ndx=Vn(i);
    All(ndx,1)=SuccessN(i);
    All(ndx,2)=1;
end
for i=1:length(V0)
    ndx=V0(i);
    if Success0(i,1)>Success0(i,2)
        All(ndx,1)=Success0(i,1);
        All(ndx,2)=-1;
    else
        All(ndx,1)=Success0(i,2);
        All(ndx,2)=1;
    end
end


