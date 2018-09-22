function REP = MOPSOCD(params,MultiObj)

    % Parameters
    Np      = params.Np;
    Nr      = params.Nr;
    maxgen  = params.maxgen;
    W       = params.W;
    C1      = params.C1;
    C2      = params.C2;
    %ngrid   = params.ngrid;
    maxvel  = params.maxvel;
    %u_mut   = params.u_mut;
    fun     = MultiObj.fun;
    nVar    = MultiObj.nVar;
    var_min = MultiObj.var_min(:);
    var_max = MultiObj.var_max(:);


    
     
    d=nVar;   %�����ռ�ά����δ֪�������� 

    xmin=var_min;
    xmax=var_max;
    xrange=xmax-xmin;


    n=Np;  %��ʼ��Ⱥ�������Ŀ 
    MaxDT=maxgen;    %���������� 

    rp=[];  %�ⲿ�浵������������ÿ�����������ķ�֧��⣬Ҳ���㷨���Ľ�� 
    rpmax=Nr;  %������ģ 
    c1=C1;c2=C2;
    w=W;   %����Ȩ�� 
    % pbx n*d   ����---�������Ž� 
    % gbx n*d   ����---ȫ�����Ž⣨ÿ�����ӵ�ȫ�����Žⶼ��ͬ�� 
    % pbf n*k   ����----��������ֵ 
    % gbf n*k   ����----ȫ������ֵ 

    %������� n ������,��ʼ�����ӵ�λ�ú��ٶ� 
    %x=xmin+xrange.*rand(n,d);

    x=[];
    for i=1:n
        x(i,:)=xmin'+xrange'.*rand(1,d);
    end

    v=zeros(n,d); 

    maxvel   = (var_max-var_min).*maxvel./100;

    tmp = fun(x);
    %�������ӵ���Ӧ�� 
    k=size(tmp,2);  %Ŀ�꺯������
    x(:,2*d+1:2*d+k)=tmp; 
    x(:,d+1:2*d)=v; 
    
    %����弫ֵ��ȫ�ּ�ֵ 
    pbx=x(:,1:d);   %��������λ��Ϊ��ʼ���ӱ��� 
    pbf=x(:,2*d+1:2*d+k); 
    gbx=x(:,1:d);   %ȫ������λ��Ϊ��ʼ���ӱ��� 
    gbf=x(:,2*d+1:2*d+k); 
    xk=x(:,1:d); 
    v=x(:,(d+1):2*d); 

    rp = x;
    
    
    % Plotting and verbose
    if(k==2)
        h_fig = figure(1);
        h_par = plot(x(:,2*d+1),x(:,2*d+2),'or'); hold on;
        h_rep = plot(rp(:,2*d+1),rp(:,2*d+2),'ok'); hold on;
        try
            set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)');
            axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                  min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
            grid on; xlabel('f1'); ylabel('f2');
        end
        drawnow;
    end
    if(k==3)
        h_fig = figure(1);
        h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
        h_rep = plot3(REP.pos_fit(:,1),REP.pos_fit(:,2),REP.pos_fit(:,3),'ok'); hold on;
        try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)','ztick',REP.hypercube_limits(:,3)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
        end
        grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
        drawnow;
        axis square;
    end
    display(['Generation #0 - Repository size: ' num2str(size(rp,1))]);
    
    
    
    
    

    %��ʼѭ�� 
    for t=1:MaxDT 
        %��������Ⱥ��ʽ������λ�á��ٶȸ��� 
        gen = t;
        for i=1:n           
            %v(i,:)=0.3.*randn(1,d)+c1.*rand.*(pbx(i,:)-xk(i,:))+c2.*rand.*(gbx(i,:)-xk(i,:)); 
            v(i,:)=w.*v(i,:)+c1.*rand.*(pbx(i,:)-xk(i,:))+c2.*rand.*(gbx(i,:)-xk(i,:)); 
            for j=1:d 
                if v(i,j)>maxvel(j) 
                    v(i,j)=maxvel(j); 
                elseif v(i,j)<-maxvel(j) 
                    v(i,j)=-maxvel(j); 
                end 
            end 
            xk(i,:)=xk(i,:)+v(i,:); 
            for j=1:d 
                if xk(i,j)>xmax(j) 
                    xk(i,j)=xmax(j);
                elseif xk(i,j)<xmin(j) 
                    xk(i,j)=xmin(j);
                end 
            end 
        end 


        for i=1:n   
            mu=0.1;             % Mutation Rate
            pm=(1-(t-1)/(MaxDT-1))^(1/mu);
            if rand<pm
                xk(i,:) = Mutate(xk(i,:),pm,xmin,xmax);
            end    
        end
        
        
   
        
        
        x(:,1:d)=xk; 
        x(:,(d+1):2*d)=v; 
        %�������ӵ���Ӧ�� 
        x(:,2*d+1:2*d+k)=fun(xk); 

        %���֧��⼯ np���������򷨣� 
        q=x;np=[]; 
        while isempty(q)==0 
            xsb=q(1,:); 
            q(1,:)=[]; 
            flag=1;    %�� flag==1���� xsb �Ƿ�֧��⣬�Ϳ��Էŵ���֧��⼯���� 
            [column,row]=size(q); 
            for i=1:column 
                dom_less=0; 
                dom_equal=0; 
                dom_more=0; 
                for l=1:k 
                    if xsb(2*d+l)<q(i,2*d+l) 
                        dom_less=dom_less+1; 
                    elseif xsb(2*d+l)==q(i,2*d+l) 
                        dom_equal=dom_equal+1; 
                    else 
                        dom_more=dom_more+1; 
                    end 
                end 

                %�� xsb ֧�� q(i)���� q(i)���ڵ��б��Ϊ 0�����������ɾ����������� 
                if dom_more==0&dom_equal~=k   
                    q(i,:)=0;  
                %�� q(i)֧�� xsb����ô xsb ���Ƿ�֧��⣬���ܱ������֧��⼯
                elseif dom_less==0&dom_equal~=k  
                    flag=0;break 
                end 
            end 
            if flag==1     %�� xsb �Ƿ�֧��⣬��������֧��⼯ np �� 
                np=[np;xsb]; 
            end 
            %�� q �б��Ϊ 0 ����ɾ����ʣ�²��� xsb ֧������ӣ�����������ļ��֮�� 
            q(~any(q,2),:)=[];  
        end 

        %���¸��弫ֵ������ǰλ��֧������弫ֵλ�ã������Ϊ��ǰλ�ã� 
        for i=1:n 
            dom_less=0; 
            dom_equal=0; 
            dom_more=0; 
            for l=1:k 
                if x(i,2*d+l)<pbf(i,:) 
                    dom_less=dom_less+1; 
                elseif x(i,2*d+l)==pbf(i,:) 
                    dom_equal=dom_equal+1; 
                else 
                    dom_more=dom_more+1; 
                end 
            end 
            if dom_more==0&dom_equal~=k 
                pbf(i,:)=x(i,2*d+1:2*d+k); 
                pbx(i,:)=x(i,1:d); 
            end 
        end 
        %�����ⲿ�� rp(����֧�伯��֧���ϵ�����ⲿ����
        [column,row]=size(rp); 
        if column==0 
            rp=np; 
        else 
            [column2,row2]=size(np); 
            for i=1:column2 
                flag=1;   %�� flag==1���� np(i,:)�ǾͿ��Էŵ��ⲿ������ 
                [column1,row1]=size(rp); 
                for j=1:column1 
                    dom_less=0; 
                    dom_equal=0; 
                    dom_more=0; 
                    for l=1:k 
                        if np(i,2*d+l)<rp(j,2*d+l) 
                            dom_less=dom_less+1; 
                        elseif np(i,2*d+l)==rp(j,2*d+l) 
                            dom_equal=dom_equal+1; 
                        else 
                            dom_more=dom_more+1; 
                        end 
                    end 
                %����֧�伯�е����� np(i,:)֧���ⲿ���е����� rp(j,:)�� 
                %�� rp(j,:)���ڵ��б��Ϊ 0�����������ɾ����������� 
                    if dom_more==0&dom_equal~=k 
                        rp(j,:)=0;  
                %�� rp(j,:)֧�� np(i,:)�����ʾ np(i,:)�� rp(j,:)֧�䣬 
                %��ô np(i,:)��һ�����Ƿ�֧��⣬�Ͳ��ܱ������֧��⼯���� 
                    elseif dom_less==0&dom_equal~=k  
                        flag=0;break 
                    end 
                end 
                if flag==1     %�� flag==1���� np(i,:)�ǾͿ��Էŵ��ⲿ������ 
                                        rp=[rp;np(i,:)]; 
                end 
                rp(~any(rp,2),:)=[]; 
            end 
        end 

        np = rp;
        [column,row]=size(rp); 
        for i=1:column
            for j=i:column
                if i~=j && sum(np(i,1:d)-np(j,1:d) == 0) == d
                    rp(i,:)=0; 
                    break;
                end
            end
        end
        rp(~any(rp,2),:)=[]; 
        %����ӵ��������Ƶ�����ģ 
        %�������ⲿ�������ӵ�ӵ�����룬��ӵ��������н�������ɾ������Ķ�����壩 
        [column,row]=size(rp); 
        indexrp=[1:column]'; 
        rp=[rp indexrp];   %���� 
        rp(:,2*d+k+2)=zeros(column,1);   %��������ӵ������ 
        for i=1:k 
            rp=sortrows(rp,2*d+i);   %�����ӽ��жԵ� i ��Ŀ�꺯��ֵ������������ 
            if column>2 
                for j=2:column-1 
                    rp(j,2*d+k+2)= rp(j,2*d+k+2)+(rp(j+1,2*d+i)-rp(j-1,2*d+i))/(rp(column,2*d+i)-rp(1,2*d+i)); 
                end 
            end 
            rp(1,2*d+k+2)=rp(1,2*d+k+2)+inf; 
            rp(column,2*d+k+2)=rp(column,2*d+k+2)+inf; 
        end 
        rp=sortrows(rp,-(2*d+k+2));     %�ⲿ������ӵ�����뽵������ 
        if column>rpmax 
            rp=rp(1:rpmax,:); 
        end 

        %����ȫ�ּ�ֵ����ʱ�ⲿ���Ѹ���ӵ�����뽵�����򣬴��ⲿ������ǰ�������ѡ��) 
        for i=1:n 
        %����һ������������ⲿ��ǰ 10%����������ѡ 
            randomnumber=ceil(rand*column/10);  
            gbx(i,:)=rp(randomnumber,1:d);     %ȫ������λ�� 
            gbf(i,:)=rp(randomnumber,2*d+1:2*d+k); 
        end 
        rp(:,2*d+k+1:2*d+k+2)=[]; 
        
        
         % Plotting and verbose
        if(k==2)
            figure(h_fig); delete(h_par); delete(h_rep);
            h_par = plot(x(:,2*d+1),x(:,2*d+2),'or'); hold on;
            h_rep = plot(rp(:,2*d+1),rp(:,2*d+2),'ok'); hold on;
            try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
            end
            if(isfield(MultiObj,'truePF'))
                try delete(h_pf); end
                h_pf = plot(MultiObj.truePF(:,1),MultiObj.truePF(:,2),'.','color',0.8.*ones(1,3)); hold on;
            end
            grid on; xlabel('f1'); ylabel('f2');
            drawnow;
            axis square;
        end
        if(k==3)% �Ժ��ٸ�
            figure(h_fig); delete(h_par); delete(h_rep); 
            h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
            h_rep = plot3(REP.pos_fit(:,1),REP.pos_fit(:,2),REP.pos_fit(:,3),'ok'); hold on;
            try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)','ztick',REP.hypercube_limits(:,3)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2)) ...
                      min(REP.hypercube_limits(:,3)) max(REP.hypercube_limits(:,3))]);
            end
            if(isfield(MultiObj,'truePF'))
                try delete(h_pf); end
                h_pf = plot3(MultiObj.truePF(:,1),MultiObj.truePF(:,2),MultiObj.truePF(:,3),'.','color',0.8.*ones(1,3)); hold on;
            end
            grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
            drawnow;
            axis square;
        end
        display(['Generation #' num2str(gen) ' - Repository size: ' num2str(size(rp,1))]);
        
        
        
        
        
        
        
    end 

    REP.pos = rp(:,1:d);
    REP.pos_fit = rp(:,2*d+1:2*d+k);
    hold off;
%     hold on;
%     rpcd=rp;
%     plot(rpcd(:,2*d+1),rpcd(:,2*d+2),'^');
end

function xnew=Mutate(x,pm,VarMin,VarMax)
    nVar=numel(x);
    j=randi([1 nVar]);
    dx=pm*(VarMax(j)-VarMin(j));
    lb=x(j)-dx;
    ub=x(j)+dx;
    lb = max(lb,VarMin(j));
    lb = min(lb,VarMax(j));
    ub = max(ub,VarMin(j));
    ub = min(ub,VarMax(j));    
    xnew=x;
    xnew(j)=unifrnd(lb,ub);
end
