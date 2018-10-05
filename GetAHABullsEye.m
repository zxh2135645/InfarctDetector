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
fillBullseye(flip(segPerc(1:6)),3,4,60,420);
fillBullseye(flip(segPerc(7:12)),2,3,60,420);
fillBullseye(flip(segPerc(13:16)),1,2,45,405);

uistack(c,'top');
set(gca,'Color',[50/255 50/255 50/255]);
caxis([0, 100]);
colormap parula;
clr = colorbar;
clr.Location = 'east';
clr.Color = [1 1 1];
clr.TickLabels = {'0%', '100%'};
clr.Ticks = [0, 100];
clr.Limits = [0, 100];

color_spec = zeros(16, 3);
for i = 1:16
   if segPerc(i) > 85
       color_spec(i,:) = [0 0 0];
   else
       color_spec(i,:) = [1 1 1];
   end
end

h1 = text(0, 3.5, num2str(round(segPerc(6))), 'HorizontalAlignment', 'center', 'Color', color_spec(6,:), 'FontSize', 12);
h7 = text(0, 2.5, num2str(round(segPerc(12))), 'HorizontalAlignment', 'center', 'Color', color_spec(12,:), 'FontSize', 12);
h13 = text(0, 1.5, num2str(round(segPerc(16))), 'HorizontalAlignment', 'center', 'Color', color_spec(16,:), 'FontSize', 12);
h15 = text(0, -1.5, num2str(round(segPerc(14))), 'HorizontalAlignment', 'center', 'Color', color_spec(14,:), 'FontSize', 12);
h10 = text(0, -2.5, num2str(round(segPerc(9))), 'HorizontalAlignment', 'center', 'Color', color_spec(9,:), 'FontSize', 12);
h4 = text(0, -3.5, num2str(round(segPerc(3))), 'HorizontalAlignment', 'center', 'Color', color_spec(3,:), 'FontSize', 12);
h16 = text(1.5, 0, num2str(round(segPerc(13))), 'HorizontalAlignment', 'center', 'Color', color_spec(13,:), 'FontSize', 12);
h14 = text(-1.5, 0, num2str(round(segPerc(15))), 'HorizontalAlignment', 'center', 'Color', color_spec(15,:), 'FontSize', 12);
h12 = text(2.1, 1.4, num2str(round(segPerc(7))), 'HorizontalAlignment', 'center', 'Color', color_spec(7,:), 'FontSize', 12);
h6 = text(2.9, 2.0, num2str(round(segPerc(1))), 'HorizontalAlignment', 'center', 'Color', color_spec(1,:), 'FontSize', 12);
h8 = text(-2.1, 1.4, num2str(round(segPerc(11))), 'HorizontalAlignment', 'center', 'Color', color_spec(11,:), 'FontSize', 12);
h2 = text(-2.9, 2.0, num2str(round(segPerc(5))), 'HorizontalAlignment', 'center', 'Color', color_spec(5,:), 'FontSize', 12);
h11 = text(2.1, -1.4, num2str(round(segPerc(8))), 'HorizontalAlignment', 'center', 'Color', color_spec(8,:), 'FontSize', 12);
h5 = text(2.9, -2.0, num2str(round(segPerc(2))), 'HorizontalAlignment', 'center', 'Color', color_spec(2,:), 'FontSize', 12);
h9 = text(-2.1, -1.4, num2str(round(segPerc(10))), 'HorizontalAlignment', 'center', 'Color', color_spec(10,:), 'FontSize', 12);
h3 = text(-2.9, -2.0, num2str(round(segPerc(4))), 'HorizontalAlignment', 'center', 'Color', color_spec(4,:), 'FontSize', 12);
%set(h,'Rotation',90);
end