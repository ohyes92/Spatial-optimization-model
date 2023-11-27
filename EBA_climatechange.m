%% 생태적 리질리언스 계획은 통합적인 기후변화 완화 및 적응 전략에 시너지 효과를 가지는가?
% Case Study : Suwon
% Decision varialbe : urban park
% Fiteness value : ecological resilience, reduction of urban, carbon stock
% Constratint : Cost.


%% Read data (%100:urban, 200:agriculture, 300:forest, 400: grass, 500:wetland, 600:bare, 700:water, 800:New habitat)
carbon = dlmread('E:\논문투고\EBA\Fitness_V/carbonv.txt');
cost = dlmread('E:\논문투고\EBA\Fitness_V/costv.txt');
land = dlmread('E:\논문투고\EBA\Fitness_V/landv.txt');
lst = dlmread('E:\논문투고\EBA\Fitness_V/lstv.txt');

%% remove NaN value
for i = 1: numel(carbon)
    if carbon(i) == -9999
       carbon(i) = 0;
    end
    if cost(i) == -9999
       cost(i) = 0;
    end
    if land(i) == -9999
       land(i) = 0;
    end
    if lst(i) == -9999
       lst(i) = 0;
    end
end


%% Constraint(비용설정)
% 비용 제한에 따른 시나리오 설정
sc_k= 0.2; %설치개수 (0.2, 0.4, 0.6) 

%% optimization
% Initialization

MaxIt=100;
nPop=1000;
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
empty_individual.Best=[];
pop=repmat(empty_individual,nPop,1);

nh = find(land==600); % 새로운 녹지 설치할 수 있는 공간 위치
Num_nh = numel(nh); % 새로운 녹지 설치할 수 있는 공간 수
domain = land; % Case study, 도메인 설정

for i=1:nPop
    pop(i).Position=domain;
    num_applied_nh = round(Num_nh*sc_k); %설치개수 (0.2, 0.4, 0.6) 
    location_nh=datasample(nh,num_applied_nh); %새로운 서식지 설치할 공간 찾기
    pop(i).Position(location_nh)=800; %새로운 서식지 설치하기
    
    % Edge Density
    ED_data = pop(i).Position;
    ED_data(ED_data == 300 | ED_data == 800 | ED_data == 900) = 10;
    ED_data(ED_data ~= 10) = 0;
    forest = ED_data;

    EDB_data = bwlabel(ED_data,4); % make habitat patches
    m=max(EDB_data,[],'all'); m; %max 값 찾기

    for edb = 1:m
            [x,y]=find(EDB_data==edb);
            k=[x y];
            p=numel(find(x == min(x,[],'all')));
            q=numel(find(y == min(y,[],'all')));
            length_B=zeros(1,m);
            length_B(edb)=(p*2)+(q*2);
    end
    
    ED = 1/(sum(length_B)/numel(domain)); % Boundary density
    
    % Connectivity
    for con = 1:m
        [x,y]=find(EDB_data==con);
        k=round(numel(x)/2);
        coor=zeros(m,2);
        coor(con,1) = x(k);
        coor(con,2) = y(k);
    end
    DC=sum(pdist(coor))/numel(domain);
     
    % Carbon stock
    CS = sum(carbon(pop(i).Position==800));
    
    % Land Surface Temperature (ha = m2*0.0001)
    % analy_lst = lst; % LST도메인 설정

 
    for lk=1:m
        habitat_area= numel(find(EDB_data==lk))*0.0001;
        buffer_area = 26.8*log(habitat_area)+94.172; %buffer area = m단위, habitat_area= ha단위
        cooling_extent = 0.152*log(habitat_area)+0.156;
        
        % urban_area
        [x,y]=find(EDB_data==lk);
        min_x=x(1); max_x=x(end); min_y=y(1); max_y=y(end);
        xv=(min_x:max_x); yv=(min_y:max_y);
        edge_ax=repmat(min_x,numel(yv),1); %윗면x
        edge_ay=(min_y:max_y); edge_ay=edge_ay.';%윗면y
        edge_up=[edge_ax,edge_ay];
        
        edge_bx=repmat(max_x,numel(yv),1); %밑면x
        edge_by=(min_y:max_y); edge_by=edge_by.';%밑면y
        edge_donw=[edge_bx,edge_by];
        
        edge_cx=(min_x+1:max_x-1); edge_cx=edge_cx.';%왼쪽면x
        edge_cy=repmat(min_y,numel(xv)-2,1); %왼쪽면y
        edge_left=[edge_cx,edge_cy];
        
        edge_dx=(min_x+1:max_x-1); edge_dx=edge_dx.';%오른쪽면x
        edge_dy=repmat(max_y,numel(xv)-2,1); %오른쪽면y
        edge_right=[edge_dx,edge_dy];
        
        buffer_d=fix(buffer_area/30); %버퍼 총 몇칸까지 가는지 확인하기
        count_b=zeros(1,buffer_d);
        for buf=1:buffer_d
            edge_up=[edge_ax-buf,edge_ay];
            edge_down=[edge_bx+buf,edge_ay];
            edge_left=[edge_cx,edge_cy-buf];
            edge_right=[edge_dx,edge_dy+buf];
            num_up = numel(find(land(edge_up)==100));
            num_down = numel(find(land(edge_down)==100));
            num_left = numel(find(land(edge_left)==100));
            num_right = numel(find(land(edge_right)==100));
            count_b(buf)=num_up+num_down+num_left+num_right;           
        end
        count_buffer=sum(count_b);
        if count_buffer < 1
            urban_area=0.001;
        else
            urban_area=count_buffer;
        end
        Cooling=zeros(1,lk);
        Cooling(lk) = (2.52* urban_area * cooling_extent)/(1/buffer_area);
    end
    CE=sum(Cooling); %cooling effect

    

    z=[ED DC CS CE]';
    
    pop(i).Cost=z;
end

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);


%% NSGA-II Main Loop

for it=1:MaxIt
    %% Crossover
     nCrossover=nPop*0.6;
     popc=repmat(empty_individual,nCrossover,1);
     for s=1:nCrossover/2
         random_num=datasample(1:nPop,2,'Replace',false);
         p=random_num(1); q=random_num(2);
         c1=pop(p).Position;
         c2=pop(q).Position;
         
         c1_nh=datasample(find(c1==800),round(numel(find(c1==800))/2),'Replace',false);
         c2_nh=datasample(find(c2==800),round(numel(find(c2==800))/2),'Replace',false);
         
         c1(c1_nh)=600; c1(c2_nh)=800;
         c2(c2_nh)=600; c2(c1_nh)=800; 
         
         popc(s+(s-1)).Position=c1;
         popc(s+s).Position=c2;
     end
        
                         
    %% Mutation
     nMutation=nPop*0.4;
     popm=repmat(empty_individual,nMutation,1);
     for k = 1:nMutation
         t=datasample(1:nPop,1);
         m1=pop(t).Position;
         nMu_remove=datasample(1:numel(find(m1==800)),round(numel(find(m1==800))*0.001));
         nMu_new=datasample(1:numel(find(m1==600)),round(numel(find(m1==800))*0.001));
         m1(nMu_remove)=600; m1(nMu_new)=800;
         popm(k).Position=m1;
     end

    %% Crossover에 대한 cost 분석
    % Edge Density
    for i =1:nCrossover
        ED_data=popc(i).Position;
        ED_data(ED_data == 300 | ED_data == 800 | ED_data == 900) = 10;
        ED_data(ED_data ~= 10) = 0;
        forest = ED_data;
        EDB_data = bwlabel(ED_data,4); % make habitat patches
        m=max(EDB_data,[],'all'); m; %max 값 찾기
    
        for edb = 1:m
                [x,y]=find(EDB_data==edb);
                k=[x y];
                p=numel(find(x == min(x,[],'all')));
                q=numel(find(y == min(y,[],'all')));
                length_B=zeros(1,m);
                length_B(edb)=(p*2)+(q*2);
        end  
        ED = 1/(sum(length_B)/numel(domain)); % Boundary density
       
        % Connectivity
        for con = 1:m
            [x,y]=find(EDB_data==con);
            k=round(numel(x)/2);
            coor=zeros(m,2);
            coor(con,1) = x(k);
            coor(con,2) = y(k);
        end
        DC=sum(pdist(coor))/numel(domain);

        % Carbon stock
        CS = sum(carbon(pop(i).Position==800));       
   
        % Land Surface Temperature (ha = m2*0.0001)
        % analy_lst = lst; % LST도메인 설정
        
        for lk=1:m
            habitat_area= numel(find(EDB_data==lk))*0.0001;
            buffer_area = 26.8*log(habitat_area)+94.172; %buffer area = m단위, habitat_area= ha단위
            cooling_extent = 0.152*log(habitat_area)+0.156;
        
        % urban_area
            [x,y]=find(EDB_data==lk);
            min_x=x(1); max_x=x(end); min_y=y(1); max_y=y(end);
            xv=(min_x:max_x); yv=(min_y:max_y);
            edge_ax=repmat(min_x,numel(yv),1); %윗면x
            edge_ay=(min_y:max_y); edge_ay=edge_ay.';%윗면y
            edge_up=[edge_ax,edge_ay];
        
            edge_bx=repmat(max_x,numel(yv),1); %밑면x
            edge_by=(min_y:max_y); edge_by=edge_by.';%밑면y
            edge_donw=[edge_bx,edge_by];
        
            edge_cx=(min_x+1:max_x-1); edge_cx=edge_cx.';%왼쪽면x
            edge_cy=repmat(min_y,numel(xv)-2,1); %왼쪽면y
            edge_left=[edge_cx,edge_cy];
        
            edge_dx=(min_x+1:max_x-1); edge_dx=edge_dx.';%오른쪽면x
            edge_dy=repmat(max_y,numel(xv)-2,1); %오른쪽면y
            edge_right=[edge_dx,edge_dy];
        
            buffer_d=fix(buffer_area/30); %버퍼 총 몇칸까지 가는지 확인하기
            count_b=zeros(1,buffer_d);
            for buf=1:buffer_d
                edge_up=[edge_ax-buf,edge_ay];
                edge_down=[edge_bx+buf,edge_ay];
                edge_left=[edge_cx,edge_cy-buf];
                edge_right=[edge_dx,edge_dy+buf];
                num_up = numel(find(land(edge_up)==100));
                num_down = numel(find(land(edge_down)==100));
                num_left = numel(find(land(edge_left)==100));
                num_right = numel(find(land(edge_right)==100));
                count_b(buf)=num_up+num_down+num_left+num_right;           
            end
            count_buffer=sum(count_b);
            if count_buffer < 1
                urban_area=0.001;
            else
                urban_area=count_buffer;
            end
            Cooling=zeros(1,lk);
            Cooling(lk) = (2.52* urban_area * cooling_extent)/(1/buffer_area);
        end
        CE=sum(Cooling); %cooling effect 
        
        z=[ED DC CS CE]';
        popc(i).Cost=z;
    end
    
    %% Mutation에 대한 cost 분석
    % Edge Density
    for i =1:nMutation
        ED_data=popm(i).Position;
        ED_data(ED_data == 300 | ED_data == 800 | ED_data == 900) = 10;
        ED_data(ED_data ~= 10) = 0;
        forest = ED_data;
        EDB_data = bwlabel(ED_data,4); % make habitat patches
        m=max(EDB_data,[],'all'); m; %max 값 찾기
    
        for edb = 1:m
                [x,y]=find(EDB_data==edb);
                k=[x y];
                p=numel(find(x == min(x,[],'all')));
                q=numel(find(y == min(y,[],'all')));
                length_B=zeros(1,m);
                length_B(edb)=(p*2)+(q*2);
        end
    
        ED = 1/(sum(length_B)/numel(domain)); % Boundary density
       
        
        % Connectivity
        for con = 1:m
            [x,y]=find(EDB_data==con);
            k=round(numel(x)/2);
            coor=zeros(m,2);
            coor(con,1) = x(k);
            coor(con,2) = y(k);
        end
        DC=sum(pdist(coor))/numel(domain);
 
        % Carbon stock
        CS = sum(carbon(pop(i).Position==800));       
   
        % Land Surface Temperature (ha = m2*0.0001)
        % analy_lst = lst; % LST도메인 설정
        tic
        for lk=1:m
            habitat_area= numel(find(EDB_data==lk))*0.0001;
            buffer_area = 26.8*log(habitat_area)+94.172; %buffer area = m단위, habitat_area= ha단위
            cooling_extent = 0.152*log(habitat_area)+0.156;
        
        % urban_area
            [x,y]=find(EDB_data==lk);
            min_x=x(1); max_x=x(end); min_y=y(1); max_y=y(end);
            xv=(min_x:max_x); yv=(min_y:max_y);
            edge_ax=repmat(min_x,numel(yv),1); %윗면x
            edge_ay=(min_y:max_y); edge_ay=edge_ay.';%윗면y
            edge_up=[edge_ax,edge_ay];
        
            edge_bx=repmat(max_x,numel(yv),1); %밑면x
            edge_by=(min_y:max_y); edge_by=edge_by.';%밑면y
            edge_donw=[edge_bx,edge_by];
        
            edge_cx=(min_x+1:max_x-1); edge_cx=edge_cx.';%왼쪽면x
            edge_cy=repmat(min_y,numel(xv)-2,1); %왼쪽면y
            edge_left=[edge_cx,edge_cy];
        
            edge_dx=(min_x+1:max_x-1); edge_dx=edge_dx.';%오른쪽면x
            edge_dy=repmat(max_y,numel(xv)-2,1); %오른쪽면y
            edge_right=[edge_dx,edge_dy];
        
            buffer_d=fix(buffer_area/30); %버퍼 총 몇칸까지 가는지 확인하기
            count_b=zeros(1,buffer_d);
            for buf=1:buffer_d
                edge_up=[edge_ax-buf,edge_ay];
                edge_down=[edge_bx+buf,edge_ay];
                edge_left=[edge_cx,edge_cy-buf];
                edge_right=[edge_dx,edge_dy+buf];
                num_up = numel(find(land(edge_up)==100));
                num_down = numel(find(land(edge_down)==100));
                num_left = numel(find(land(edge_left)==100));
                num_right = numel(find(land(edge_right)==100));
                count_b(buf)=num_up+num_down+num_left+num_right;           
            end
            count_buffer=sum(count_b);
            if count_buffer < 1
                urban_area=0.001;
            else
                urban_area=count_buffer;
            end
            Cooling=zeros(1,lk);
            Cooling(lk) = (2.52* urban_area * cooling_extent)/(1/buffer_area);
        end
        CE=sum(Cooling); %cooling effect 
        tok
        z=[ED DC CS CE]'; 
        popm(i).Cost=z;
    end    
        
    
     %% Merge
     pop=[pop
          popc
          popm]; %#ok
     
     % Non-Dominated Sorting
     [pop, F]=NonDominatedSorting(pop);

     % Calculate Crowding Distance
     pop=CalcCrowdingDistance(pop,F);

     % Sort Population
     pop=SortPopulation(pop);
    
     % Truncate
     pop=pop(1:nPop);
    
     % Non-Dominated Sorting
     [pop, F]=NonDominatedSorting(pop);

     % Calculate Crowding Distance
     pop=CalcCrowdingDistance(pop,F);

     % Sort Population
     [pop, F]=SortPopulation(pop);
    
     % Store F1
     F1=pop(F{1});
    
     % Show Iteration Information
     disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
     tiledlayout(2,3)

     v=genvarname('solution',who);
     eval([v ' = table(F1.Cost)']); 
     
     fitness=[F1.Cost];
     iter(it)=sum(fitness(1,:)+fitness(2,:)+fitness(3,:)+fitness(4,:));
 
end

