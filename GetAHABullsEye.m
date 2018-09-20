function GetAHABullsEye(name, alg_label)
if nargin == 1
    alg_label = 'LGE';
end
base_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
file_path = cat(2, base_dir, 'rocForReal.csv');
df = xlsread(file_path);
Names = GetSubjectName(base_dir);
RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);
%name = 'OH_SANG_MOO';
idx = find(strcmp(name,Names));
%alg_label = 'Mean5SD';
switch alg_label
    case 'LGE'
        segPerc = df((16*idx - 15):(16*idx), 1);
    case 'Mean5SD'
        segPerc = df((16*idx - 15):(16*idx), 2);
    case 'Otsu'
        segPerc = df((16*idx - 15):(16*idx), 3);
    case 'Kmeans'
        segPerc = df((16*idx - 15):(16*idx), 4);
    case 'GMM'
        segPerc = df((16*idx - 15):(16*idx), 5);
end

% Create the AHA 17-segment bullseye
figure();
c = createBullseye([0 1 1 0; 1 2 4 45; 2 3 6 0; 3 4 6 0]);
set(c,'Color','w','LineWidth',1)

% Filling the bullseye, vector by vector
fillBullseye(segPerc(1:4),1,2,-225, 135);
fillBullseye(segPerc(5:10),2,3,0,360);
fillBullseye(segPerc(11:16),3,4,0,360);
uistack(c,'top');
set(gca,'Color',[50/255 50/255 50/255]);
caxis([0, 40]);
colormap parula;
clr = colorbar;
clr.Location = 'east';
clr.Color = [1 1 1];
clr.TickLabels = {'0%', '40%'};
clr.Ticks = [0, 40];
clr.Limits = [0, 40];


end