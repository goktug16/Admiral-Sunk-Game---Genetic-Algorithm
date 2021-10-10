clear all;
% fitness1 : Setting the ball throws to make the least angle change of the ball after each throw
% fitness2 : After each shot, we want the artillery pattern to be selected, in which the ships and bridges with the highest value of the next cannon shot are hit.
% We want to minimize the  fitness 1 function
% We want to maximize the  fitness 2 function

% possible moves [0 8] 
% 0 stand in same position
% 1-8 directions

hareket=[0 0; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1 ; 1 0; 1 1];

% cost between consecutive moves:
ac_t=[2:9];
ac=(ac_t-2)*45; ac=[0 ac]; % 0 0 45 90 135 180 225 270 315

G=200; % number of generations 
P=100; % population size 
cross=2; % crossover 1: single point, crossover 2: from two points
bas=[1 1]; % starting point

prompt1 = 'Enter the population value';
prompt2 = 'Enter the generation value';
prompt3 = 'Enter the crossover type : 1 or 2';
prompt4 = 'Starting point : ';

P = input(prompt1);
G = input(prompt2);
cross = input(prompt3);
bas= input(prompt4);

Mx=9;My=Mx; % matrix size
M=zeros(Mx,My);% M matrix
u=Mx*Mx-1; % len of each pattern
mu=0.01;  % mutation rate 
BK= P - P/4; % the best number of individuals directly copied to the next generation
B=zeros(P,u); % patterns 
total_points = 0;

max_f1 = 180*(u-1);% The maximum value of the fitness1 function determined to be used in normalization
max_f2 = 345; % The maximum value of the fitness2 function determined to be used in normalization

bb = zeros(1,G); %holds the  value of the best individual of each generation
bb_yol = zeros(G,u); % holds the firing pattern of the best individual of each generation

ob=zeros(1,G); % holds the average value of the each generation
j=1;
index_value = 0;
tic

original_map = zeros(Mx,My);

original_map(1,1) = 15;  % point values of enemy ship coordinates and assignment of points to enemy ship coordinates
original_map(1,2) = 15;
original_map(1,3) = 15;

original_map(1,5) = -10; % point values to be lost if the friendly ship is hit and point assignments to the coordinates where the friendly ship is located
original_map(1,6) = -10;
original_map(1,7) = -10;

original_map(3,6) = 15;  %  point values of enemy ship coordinates and assignment of points to enemy ship coordinates
original_map(3,5) = 15;


original_map(5,6) = 25;  % strategic bridge coordinates and the points to be earned by the accurate hit
original_map(6,6) = 25;
original_map(7,6) = 25;

original_map(9,4) = 15;  % point values of enemy ship coordinates and assignment of points to enemy ship coordinates
original_map(9,5) = 15;
original_map(9,6) = 15;

original_map(9,7) = -10; % point values to be lost if the friendly ship is hit and point assignments to the coordinates where the friendly ship is located
original_map(9,8) = -10;

original_map(5,1) = 25;  % strategic bridge coordinates and the points to be earned by the accurate hit
original_map(6,1) = 25;
original_map(7,1) = 25;

original_map(5,9) = 25;  % strategic bridge coordinates and the points to be earned by the accurate hit
original_map(6,9) = 25;
original_map(7,9) = 25;

original_map(8,1) = -10; % point values to be lost if the friendly ship is hit and point assignments to the coordinates where the friendly ship is located
original_map(9,1) = -10;

B=round(8*rand(P,u)); % fill with numbers between 0 and 8

for i=1:G
    f1 = zeros(1,P);
    f2 = f1(1:P);
    
    for j=1:P
        index_value = 0;
        M = original_map;
        birey=B(j,:); % movement
    
        kxy = bas;
        M(kxy(1),kxy(2))=1;  % at start 
        
        for k = 1:u % bireylerin hareketlerinin hesaplanmasi
            p_kxy =kxy +hareket(birey(k)+1,:);
             if p_kxy(1)<1 || p_kxy(2)<1 || p_kxy(1)>Mx  || p_kxy(2)>Mx % stay in place if it goes out
             else % didn't go out, move
                    kxy = p_kxy;
                    index_value = index_value + M(kxy(1),kxy(2));
             end
              M(kxy(1),kxy(2))=0; % mark the destination on the map
        end
        f2(j) = index_value;
        index_value = 0;
         %Calculation of deflection angle for the next ball toss
        birey_arti=birey+1;
        birey_2=birey_arti(2:end);
        birey_1=birey_arti(1:end-1);
        a1=ac(birey_1); 
        a2=ac(birey_2);
        fark=abs(a1-a2);
        fark(fark>180)=360-fark(fark>180);
        f1(j)=sum(fark);
    end
    % normalization for fitness functions 
    n_f1 = f1/max_f1;
    n_f2 = f2/max_f2;
    n_f2 = 1 - n_f2;
    
    w =n_f1 + n_f2; % sum of fitness values
    n_w = w/sum(w); % Mins must have a high chance of being selected in n_w
    yollar = B(:,:);
    n_w =1-n_w; % highlighting the most wanted cases that we want to be selected from the functions by inverting
    n_w=n_w/sum(n_w);
    % roulette wheel
    % using the normalized fitness functions
    [sorted,inds]=sort(n_w);
    rn_w(inds)=1:P;
    rn_w=double(rn_w)/double(sum(rn_w));
    [val best_ind]=max(rn_w);
     %best_ind
    bb(i)=w(best_ind); % the best value is chosen
    bb_yol(i,:)=B(best_ind,:); % the firing order of the one with the best value is chosen
    ob(i)=mean(w);
    secilenler = randsample(P,P,true,rn_w);% According to the values of the indices, indices with a high fitness function result are selected so that the chance of selection is higher.
    
    % breeding new individuals single-point crossover
    YB=zeros(P,u); % new 
    for j=1:P/2
        b1=B(secilenler(j),:);
        b2=B(secilenler(j+(P/2)-1),:);
        if cross==1 % single-point crossover
            kesme=round((u-3)*rand(1,1))+2; % Number between 2 - (u-1)
            YB(j,:) = [b1(1:kesme) b2(kesme+1:end)];
            YB(j+(P/2),:)=[b2(1:kesme) b1(kesme+1:end)];
        else
            kesme=round((u-3)*rand(1,2))+2; % 2 - (u-1) 2 points
            kesme=sort(kesme); % sort from smallest to largest
            YB(j,:)      =[b1(1:kesme(1)) b2(kesme(1)+1:kesme(2)) b1(kesme(2)+1:end) ];
            YB(j+(P/2),:)=[b2(1:kesme(1)) b1(kesme(1)+1:kesme(2)) b2(kesme(2)+1:end) ];
        end
    end
    if BK>0 % Copy the best BK value in B to YB
        YB(inds(BK+1:end),:)=B(inds(BK+1:end),:); 
    end
    % apply mutation
    d_ind=rand(P,u)<mu; % cells to change
    yy=round(8*rand(P,u)); % what will they change
    YB(d_ind)=yy(d_ind);
    B=YB; % new generation ready 
end

plot(bb,'r'); % draw the best
hold on;
plot(ob); % draw the avg

xlabel('Jenerasyon Sayisi');
ylabel('Fitness Degeri');
hold on 
plot(ob,'-g'); % draw the avg

[v bb_best] = max(bb);
birey = bb_yol(bb_best,:);
rota =zeros(2,u);
M=original_map;
kxy = bas;
M(kxy(1),kxy(1))=1; % at start 

for k=1:u  % obtaining the most successful firing order
    p_kxy=kxy+hareket(birey(k)+1,:);
    if  p_kxy(1)<1 || p_kxy(2)<1 || p_kxy(1)>Mx  || p_kxy(2)>Mx 
    else 
        kxy=p_kxy;
    end
        total_points = total_points + M(kxy(1),kxy(2));
        M(kxy(1),kxy(2))= 0;
        rota(:,k)=kxy;
end

M=original_map;
kxy = bas;
M(kxy(1),kxy(1))=1; 

for k=1:u  
    p_kxy=kxy+hareket(birey(k)+1,:);
    if  p_kxy(1)<1 || p_kxy(2)<1 || p_kxy(1)>Mx  || p_kxy(2)>Mx 
    else  
        kxy=p_kxy;
    end
        M(kxy(1),kxy(2))= k+1;
end


for i=1:Mx   % changing values to understand colors in charts
    for j=1:My
        if original_map(i,j) < 0
            original_map(i,j) = 300;
        elseif original_map(i,j) == 15
            original_map(i,j) = 100;
        end
        if original_map(i,j) == 25
            original_map(i,j) = 150;
        end
    end
end

success_rate = (total_points / max_f2) * 100;
success_rate
bb(end) % 
figure; image(original_map); % 
figure;image((64/u)*M);
figure; scatter(rota(1,:),rota(2,:));
hold on


axis([0 Mx+1 0 Mx+1]);
