frames=91;
Nx=10;

dir='raw_pix/';
data=NaN(4,Nx*4,frames);

for nx=1:Nx

for i=1:91
    
str='induction05135t';
    
    I1=imread([dir str num2str(i,'%02u') 'xy' num2str(nx,'%02u') 'c1.tif']); % blue
    I2=imread([dir str num2str(i,'%02u') 'xy' num2str(nx,'%02u') 'c2.tif']); % green
    I3=imread([dir str num2str(i,'%02u') 'xy' num2str(nx,'%02u') 'c3.tif']); % red

    I1(:,1:end-1)=I1(:,2:end);

    load(['masks/' num2str(nx) '/mask' num2str(i,'%03u')])
    
    if i==1
        D=regionprops(mask,'Centroid');
        N=length(D);
        v=[D.Centroid];
        p0=v(1:2:length(v)-1);
        b=NaN(N,frames);
        g=NaN(N,frames);
        r=NaN(N,frames);
        s=NaN(N,frames);
        t=NaN(frames,1);
    end
    t(i)=(i-1)*10;
    
    [bi,gi,ri,si,p]=frames_info(mask,I1,I2,I3);
    
    for j=1:length(p)
        for n=1:N
            if abs(p(j)-p0(n))<30
                data(1,(nx-1)*4+n,i)=bi(j);
                data(2,(nx-1)*4+n,i)=gi(j);
                data(3,(nx-1)*4+n,i)=ri(j);
                data(4,(nx-1)*4+n,i)=si(j);
            end
        end
    end
    imwrite(mask,['mask' num2str(i,'%03u') '.jpeg'])
end

end

save('data','data')