% input scores file
MScore=dlmread('networkfiles/emt26network.phs');
Sc=zeros(26,1);
for i=1:length(MScore(:,1))
    Sc(MScore(i,1)+1)=1;
    Sc(MScore(i,2)+1)=-1;
end
tries=100;
Vp=find(Sc==1);
SuccessP=zeros(size(Vp));
com0='/Users/kessler/Work/boolean/TRANSIT/TRANS10 "emt" 26 2 400 ';
for i=1:length(SuccessP)
    [st,cout]=system([com0,num2str(tries),' 2 7 ',num2str(-(Vp(i)))]);
    SuccessP(i)=str2double(cout);
end
disp(SuccessP)
Vn=find(Sc==-1);
SuccessN=zeros(size(Vn));
for i=1:length(SuccessN)
    [~,cout]=system([com0,num2str(tries),' 2 7',num2str(Vn(i))]);
    SuccessN(i)=str2double(cout);
end
disp(SuccessN)
V0=find(Sc==0);
Success0=zeros(length(V0),2);
for i=1:length(V0)
    [~,cout]=system([com0,num2str(tries),' 2 7',num2str(-(V0(i)))]);
    Success0(i,1)=str2double(cout);
    [~,cout]=system([com0,num2str(tries),' 2 7',num2str(V0(i))]);
    Success0(i,2)=str2double(cout);
end
disp(Success0)

All=zeros(26,2);
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


