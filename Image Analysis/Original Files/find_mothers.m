function mask=find_mothers(A)

amax=double(max(max(A)))/65535;
A=imadjust(A,[0.004 amax],[0 1]);

mask=false(size(A));

x=sum(double(A),1);
[bla,x_ch]=findpeaks(x(21:end-20),'MinPeakDistance',20,'MinPeakHeight',500000);
x_ch=x_ch+20;

cutox1=[0.5 0.5 0.5 0.5];
cutox2=[0.7 0.7 0.7 0.7];

y=zeros(size(A,1),1);
yy=zeros(size(A,1),1);
for j=1:length(x_ch)

    cutoff=cutox1(j);
    
    rangex=x_ch(j)+[-10:10];
    for i=1:size(A,1)    
        y(i)=max(double(A(i,rangex)));
    end
    
    % find skeleton of top cell
    yy=find(y>100/amax);
    y1=yy(1);
    [cc,dx]=max(A(y1,rangex));
    c1=cc(1);
    x1=dx(1)+x_ch(j)-11;
    
    y2=y1+1;
    [c2,dx]=max(A(y2,x1+[-1:1]));
    x2=x1+dx(1)-2;
    c=max(c1,c2(1));
    Cx=[c1 c2(1)];
    Tx=[x1 x2];
    
    cutoff=cutox2(j);
    
    while(A(y2,x2)>cutoff*c && y2<size(A,1))
        y2=y2+1;
        [c2,dx]=max(A(y2,x2+[-1:1]));
        x2=x2+dx(1)-2;
        c=max(c,c2(1));
        Cx=[Cx c2(1)];
        Tx=[Tx x2];
    end
    Tx=round(smooth(Tx,5));

    cutoff=cutox1(j);
    
    ic=find(Cx>cutoff*c);
    for i=ic:length(Tx)
        mask(y1+i-1,Tx(i))=true;
    end
    
end

se=strel('disk',2);
mask=imdilate(mask,se);
mask(size(mask,1),:)=false;

