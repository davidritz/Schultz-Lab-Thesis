function N=resistances(M,drugx,numbers)

d=length(M);
N=NaN(d,2);

drug=log10(drugx);
drug(1)=drug(2)^2/drug(3);

% figure

for i=1:d
    n=numbers(i);
    if n=='PAO1'
        name='PAO1';
    elseif n=='PA14'
        name='PA14';
    else
        name=['IPC' num2str(n)];
    end

    %drug concentration where delay is above a doubling time
    
    v=M(i,:,1);
    thresh=1/M(i,1,2);
    
    [a,b]=max(v);
    if(a<thresh) b=length(v); end
    ii=[b-5:b];
    ii=ii(ii>0);

    f1=fittype( 'exp1' );
    pd = fit( drug(ii)' , v(ii)' , f1 ,'Normalize','on');
    
    drug_range=linspace(min(drug(ii)),1+max(drug(ii)),1000);
    vx = feval( pd, drug_range );

    ix = find( vx > thresh , 1);
    x = drug_range(ix);
    
    if(isempty(x))
        N(i,1)=(10^max(drug))*1.2;
    else
        N(i,1)=10^x;
    end
        
    
    % subplot(8,ceil(d/4),i)
    % plot(drugx,v,'ko')
    % hold on
    % plot(drugx(ii),v(ii),'ro')
    % vx(vx>v(ii(end)))=NaN;
    % plot(10.^drug_range,vx,'b')
    % plot(10.^[min(drug-1) max(drug+1)],thresh*[1 1],'--k')
    % axis([drugx(1) drugx(end)+100 -2.0 21])
    % set(gca,'XScale','log')
    % if(i==1) 
    %     set(gca,'box','off','Xtick',[1 10 100 1000],'Ytick',[0 20])
    %     ylabel('Delay (hours)')
    %     xlabel('(\mug/ml')
    % else
    %     axis off
    % end
    % 
    % title(name)

    % drug concentration where growth is reduced by half
    
    v=M(i,:,2);
    thresh=M(i,1,2)/2;
    
    i1=find(v==0,1);
    if isempty(i1)
        i1=length(v);
    end
    
    i1=max(i1,3);
    ii=[i1-5:i1];
    ii=ii(ii>0);
    
    f2=fittype( 'power2' );
    pg = fit( drug(ii)'+4 , v(ii)' , f2 );
    
    drug_range=linspace(min(drug(ii)),1+max(drug(ii)),1000); 
    vy = feval( pg, drug_range+4 );
    
    iy = find( vy < thresh ,1);
    y = drug_range(iy);
    
    if(isempty(y))
        N(i,2)=(10^max(drug))*1.2;%10^max(drug);
    else
        N(i,2)=10^y;
    end
    
    % subplot(8,ceil(d/4),4*ceil(d/4)+i)
    % plot(drugx,v,'ko')
    % hold on
    % plot(drugx(ii),v(ii),'ro')
    % vy(vy<0)=NaN;
    % plot(10.^drug_range,vy,'b')
    % plot(10.^[min(drug-1) max(drug+1)],thresh*[1 1],'--k')
    % axis([drug(1) drugx(end)*1.1 -0.1 1])
    % set(gca,'XScale','log')
    % if(i==1) 
    %     set(gca,'box','off','Xtick',[1 10 100 1000],'Ytick',[0 1])
    %     ylabel('Growth (doub./h)')
    %     xlabel('(\mug/ml')
    % else
    %     axis off
    % end
    % 
    % title(name)

end

N(N>(drugx(end)*1.2)) = drugx(end)*1.2;