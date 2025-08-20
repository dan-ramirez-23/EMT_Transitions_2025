% input scores file
Nodes=72;
filename=['networkfiles/emt',num2str(Nodes),'network.phs']
MScore=dlmread(filename);
Sc=zeros(Nodes,1);
for i=1:length(MScore(:,1))
    Sc(MScore(i,1)+1)=1;
    Sc(MScore(i,2)+1)=-1;
end
tries=10;
Vp=find(Sc==1);
SuccessP=zeros(size(Vp));
Vn=find(Sc==-1);
SuccessN=zeros(size(Vn));
V0=find(Sc==0);
Success0p=zeros(length(V0),1);
Success0m=zeros(length(V0),1);
com0='/Users/kessler/Work/boolean/TRANSIT/TRANS10';
com1=[com0,' "emt" 72 3 32 ']
parfor i=1:length(SuccessP)
    [st,cout]=system([com1,num2str(tries),' ',num2str(-(Vp(i)))]);
    SuccessP(i)=str2double(cout);
end
disp(SuccessP)

parfor i=1:length(SuccessN)
    [~,cout]=system([com1,num2str(tries),' ',num2str(Vn(i))]);
    SuccessN(i)=str2double(cout);
end
disp(SuccessN)

parfor i=1:length(V0)
    [~,coutp]=system([com1,num2str(tries),' ',num2str(-(V0(i)))]);
    Success0p(i)=str2double(coutp);
    [~,coutm]=system([com1,num2str(tries),' ',num2str(V0(i))]);
    Success0m(i)=str2double(coutm);
end
disp([Success0p,Success0m])
Success0=[Success0p,Success0m];

All=zeros(Nodes,2);
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


