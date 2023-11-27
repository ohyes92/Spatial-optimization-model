%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA120
% Project Title: Non-dominated Sorting Genetic Algorithm II (NSGA-II)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

%% Problem Definition

circuit_second=dlmread('C:\Users\mrsol\OneDrive\Desktop\connect_150/4class_connect1501.txt');
lst_second=dlmread('C:\Users\mrsol\OneDrive\Desktop\lst_150/exact_lst1501.txt');
unusedland_second=dlmread('C:\Users\mrsol\OneDrive\Desktop\empty lot_150/unusedland_1501.txt');
weight=dlmread('C:\Users\mrsol\OneDrive\Desktop\connect_150/weight11.txt');
seoul_domain=lst_second;
for i = 1: numel(seoul_domain)
    if seoul_domain(i)>0
        seoul_domain(i)=1;
    end
end
q=find(unusedland_second==1);
seoul_domain(q)=2;
A=seoul_domain;



%% Initialization

MaxIt=50;
nPop=1000;
nCol=274;
nRow=202;
nVar=nCol*nRow;
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
pop=repmat(empty_individual,nPop,1);

NV = 30;  % number of new vegetation
k=find(unusedland_second==1);
for i=1:nPop
    pop(i).Position=unusedland_second;
    p=datasample(k,NV,'Replace',false);
    pop(i).Position(p)=4;
end
 for i = 1:nPop
     z1 = -sum(circuit_second(pop(i).Position==4).*weight(pop(i).Position==4));
     z2 = -sum(lst_second(pop(i).Position==4));
     pop(i).Cost=[z1 z2];
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
     nCr=8;
     popc=repmat(empty_individual,nCrossover,1);
     for s=1:nCrossover
         q=datasample(1:100,1);
         c1=pop(q).Position;
         ck1=find(c1==1);
         cp=datasample(ck1,nCr,'Replace',false);
         c1(cp)=4;
         m=find(c1==4);
         n=datasample(m,nCr,'Replace',false);
         c1(n)=1;
         popc(s).Position=c1;
     end
        
            
    %% Crossover(Cost)
     for i = 1:nCrossover
         z1 = -sum(circuit_second(pop(i).Position==4).*weight(pop(i).Position==4));
         z2 = -sum(lst_second(popc(i).Position==4));
         popc(i).Cost=[z1 z2];
     end
              
    %% Mutation
     nMutation=nPop*0.4;
     nMu=10;
     popm=repmat(empty_individual,nMutation,1);
     for i = 1:nMutation
         t=datasample(1:100,1);
         m1=pop(t).Position;
         mf=find(m1==4);
         mn=datasample(mf,nMu,'Replace',false);
         m1(mn)=1;
         mk=find(m1==1);
         mkn=datasample(mk,nMu,'Replace',false);
         m1(mkn)=4;
         popm(i).Position=m1;
     end
         
    %% Mutation (Cost)
     for i = 1:nMutation
         z1 = -sum(circuit_second(pop(i).Position==4).*weight(pop(i).Position==4));
         z2 = -sum(lst_second(popm(i).Position==4));
         popm(i).Cost=[z1 z2];
     end
         
     
     
     % Merge
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
    
     % Plot F1 Costs
     figure(1);
     PlotCosts(F1);
     pause(0.01);
    
end

%% Results

