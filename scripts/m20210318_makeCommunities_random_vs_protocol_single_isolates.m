%% init

clc
close all
clearvars

fig1 = figure(1);
set(fig1,'Position',get(fig1,'Position').*[0 0 0.7 0.4])

fig2 = figure(2);
set(fig2,'Position',get(fig1,'Position'))

%% parameters

Ns = 8*8; % species pool size
ns = 8; % species per community
Nc = 64; % number of communities

%% analytical

% probability of a species appearing in k communities
p_ind = @(k) nchoosek(Nc,k)*(ns/Ns)^k*(1-ns/Ns)^(Nc-k);

% subplot(3,2,1)
k = (0:1:15);
p_k_ind = arrayfun(p_ind,k);
% plot(k,p_k_ind)
% %xlabel('# of appearances')
% ylabel('probability')
% title('individual species')
% niceAxis()
% set(gca,'XTick',[])

% probability of a pair appearing in k communities
p_pairs = @(k) nchoosek(Nc,k)*(ns*(ns-1)/Ns/(Ns-1))^k*(1-ns*(ns-1)/Ns/(Ns-1))^(Nc-k);

% subplot(3,2,2)
p_k_pairs = arrayfun(p_pairs,k);
% plot(k,p_k_pairs)
% %xlabel('# of appearances')
% %ylabel('probability')
% title('species pairs')
% niceAxis()
% set(gca,'XTick',[])
% set(gca,'YTick',[])

%% true random sampling

% make communities
C = zeros(Ns,Nc);
for i=1:Nc
    comm_i = randsample(Ns,ns);
    C(comm_i,i) = 1;
end

figure(1)
makePlots(C)

%% protocols

% single plates
plate_singles = zeros(sqrt(Nc),sqrt(Nc),Ns);

for i=1:Ns
    row = ceil(i/sqrt(Nc));
    col = i - sqrt(Nc)*(row-1);
    plate_singles(row,col,i) = 1;
end

% use protocols to generate communities
[C,id] = protocol(plate_singles);

figure(2)
makePlots(C)

% plot plate
figsc = figure(3);
set(figsc,'Position',get(figsc,'Position').*[0 0 0.7 0.7])
imagesc(id)
set(gca,'XAxisLocation','top')
set(gca,'XTick',(1:1:sqrt(Nc)),'XTickLabel',(3:1:10))
set(gca,'YTick',(1:1:sqrt(Nc)),'YTickLabel',{'A','B','C','D','E','F','G','H'})
set(gca,'TickLength',[0 0])
hold all
for i=2:sqrt(Nc)
    plot([i-.5 i-.5],[.5 sqrt(Nc)+.5],'-k')
    plot([.5 sqrt(Nc)+.5],[i-.5 i-.5],'-k')
end
C_unique = unique(C','rows')';

% are there multi-sampled (repeated) species?
max_num_samples = max(max(C));

%% functions

function plate_out = rotate_left(plate)
plate_out = rot90(plate);
end
function plate_out = rotate_right(plate)
plate_out = rot90(plate,3);
end
function plate_out = shuffle_rows(plate,transformation_table)
plate_out = zeros(size(plate));
for i=1:size(plate,1)
    plate_out(transformation_table(i,2),:,:) = plate(transformation_table(i,1),:,:);
end
end

function makePlots(C)

Nc = size(C,2);
Ns = size(C,1);
ns = unique(sum(C));

% analytical: probability of a species appearing in k communities
p_ind = @(k) nchoosek(Nc,k)*(ns/Ns)^k*(1-ns/Ns)^(Nc-k);

% analytical: probability of a pair appearing in k communities
p_pairs = @(k) nchoosek(Nc,k)*(ns*(ns-1)/Ns/(Ns-1))^k*(1-ns*(ns-1)/Ns/(Ns-1))^(Nc-k);

% analytical plots
k = (0:1:Nc);
p_k_ind = arrayfun(p_ind,k);
p_k_pairs = arrayfun(p_pairs,k);

% individual species
[counts,edges] = histcounts(sum(C,2),linspace(0,Nc,Nc+1)-0.5);
counts = counts/sum(counts);
center = 0.5*(edges(1:end-1)+edges(2:end));

subplot(1,2,1)
bar(center,counts,0.7,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','BarWidth',0.8);
hold all
plot(k,p_k_ind)
xlabel('# of appearances')
ylabel('probability')
niceAxis()
title('singles')

% species pairs
pairs = zeros(Ns,Ns);
for i=1:(Ns-1)
    for j=(i+1):Ns
        for n=1:Nc
            pairs(i,j) = pairs(i,j) + double(isequal(C([i,j],n),[1;1]));
        end
    end
end
pairs = pairs(triu(true(size(pairs)),1));

[counts,edges] = histcounts(pairs,linspace(0,Nc,Nc+1)-0.5);
counts = counts/sum(counts);
center = 0.5*(edges(1:end-1)+edges(2:end));

subplot(1,2,2)
bar(center,counts,0.7,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','BarWidth',0.8);
hold all
plot(k,p_k_pairs)
xlabel('# of appearances')
%ylabel('probability')
niceAxis()
set(gca,'YTick',[])
title('pairs')

end

function niceAxis()

%{
h = findobj(gca,'Type','line');
x = get(h,'Xdata');
y = get(h,'Ydata');

dx = max(x) - min(x);
dy = max(y) - min(y);

set(gca,'XLim',[min(x) max(x)] + dx*0.15*[-1 1])
set(gca,'YLim',[min(y) max(y)] + dy*0.15*[-1 1])
set(gca,'TickLength',[0 0])
%}

set(gca,'XLim',[-1 15+1])
set(gca,'YLim',[-0.05 1.05])
set(gca,'TickLength',[0 0])

end

function [C,id] = protocol_1(plate_singles)

% transformation_table
transformation_table = [1,2,3,4,5,6,7,8;3,4,5,6,7,8,1,2]';
transformation_table = [1,2,3,4,5,6,7,8;7,8,1,2,3,4,5,6]';

% pairs
plate_pairs_1 = plate_singles + rotate_left(rotate_left(plate_singles));
plate_pairs_2 = shuffle_rows(plate_pairs_1,transformation_table);

% quadruplets
plate_quad_1 = plate_pairs_1 + plate_pairs_2;

% octuplets
plate_oct_1 = plate_quad_1 + rotate_left(plate_quad_1);

% fix structure
Ns = size(plate_singles,3);
Nc = size(plate_singles,1)*size(plate_singles,2);
C = zeros(Ns,Nc);
for i=1:Nc
    row = ceil(i/sqrt(Nc));
    col = i - sqrt(Nc)*(row-1);
    C(:,i) = squeeze(plate_oct_1(row,col,:));
end

% community IDs
C_unique = unique(C','rows')';
id = zeros(1,Nc);
for i=1:Nc
    for j=1:size(C_unique,2)
        if isequal(C(:,i),C_unique(:,j))
            id(i) = j;
        end
    end
end

id_plate = zeros(sqrt(Nc),sqrt(Nc));
for i=1:Nc
    row = ceil(i/sqrt(Nc));
    col = i - sqrt(Nc)*(row-1);
    id_plate(row,col) = id(i);
end
id = id_plate;

end

function [C,id] = protocol_2(plate_singles)

% transformation_table
transformation_table = [1,2,3,4,5,6,7,8;7,8,1,2,3,4,5,6]';

% pairs
plate_pairs_1 = plate_singles + rotate_left(plate_singles);
plate_pairs_2 = shuffle_rows(plate_pairs_1,transformation_table);

% quadruplets
plate_quad_1 = plate_pairs_1 + plate_pairs_2;
plate_quad_2 = shuffle_rows(plate_quad_1,transformation_table);

% octuplets
plate_oct_1 = plate_quad_1 + rotate_left(plate_quad_2);

% fix structure
Ns = size(plate_singles,3);
Nc = size(plate_singles,1)*size(plate_singles,2);
C = zeros(Ns,Nc);
for i=1:Nc
    row = ceil(i/sqrt(Nc));
    col = i - sqrt(Nc)*(row-1);
    C(:,i) = squeeze(plate_oct_1(row,col,:));
end

% community IDs
C_unique = unique(C','rows')';
id = zeros(1,Nc);
for i=1:Nc
    for j=1:size(C_unique,2)
        if isequal(C(:,i),C_unique(:,j))
            id(i) = j;
        end
    end
end

id_plate = zeros(sqrt(Nc),sqrt(Nc));
for i=1:Nc
    row = ceil(i/sqrt(Nc));
    col = i - sqrt(Nc)*(row-1);
    id_plate(row,col) = id(i);
end
id = id_plate;

end

function [C,id] = protocol(plate_singles)
[C,id] = protocol_2(plate_singles);
end
