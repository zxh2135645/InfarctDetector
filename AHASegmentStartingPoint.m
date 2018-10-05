clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Desktop\contour_exporting_Guan\';
out_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference'};
anatomy = anatomy_label{1};
label = sequence_label{1};

%% Get the angle
outputFileName = [out_dir, label, '_coords.csv'];
[num,txt,raw] = xlsread(outputFileName);

x = num(1,1);
y = num(1,2);
x_centroid = num(1,3);
y_centroid = num(1,4);

atan2(y - y_centroid, x - x_centroid) * 180 / pi
angle(x - x_centroid + 1i*(y - y_centroid)) * 180 / pi

atan2(y_centroid-y, x_centroid-x) * 180 / pi
angle(x_centroid-x + 1i*(y_centroid-y)) * 180 / pi

atan2(x - x_centroid, y - y_centroid) * 180 / pi % These are what we wanted
angle(y - y_centroid + 1i*(x - x_centroid)) * 180 / pi % These are what we wanted 

atan2(x_centroid-x, y_centroid-y) * 180 / pi
angle(y_centroid-y + 1i*(x_centroid-x)) * 180 / pi

%% Show if basal and apical is correct
name = 'YOO_ON';
[LGE_Myo3D, LGE_dicom_idx] = ReadMatFile3D(name, sequence_label{1}, anatomy);
figure();
for i = 1:size(LGE_Myo3D, 3)
    subplot(3,3,i)
    imagesc(LGE_Myo3D(:,:,i))
    axis equal
end

figure();
LGE_Myo3D = LGE_Myo3D(:,:,end:-1:1);
for i = 1:size(LGE_Myo3D, 3)
    subplot(3,3,i)
    imagesc(LGE_Myo3D(:,:,i))
    axis equal
end

%% Bullseye
name = 'CHOI_DAE_SUK';
GetAHABullsEye(name);

for alg_idx = 1:length(alg)
    alg_label = alg{alg_idx};
    GetAHABullsEye(name, alg_label)
end

%% Try different thresholds
[num,txt,raw] = xlsread(cat(2, out_dir, 'rocForReal.csv'));
thresh_array = 0:1:90;
LGE_labels = zeros(size(num,1), length(thresh_array));

for i = 1:length(thresh_array)
   LGE_labels(:, i) = num(:,1) > thresh_array(i); 
end

auc_t15sd = zeros(length(thresh_array), 1);
auc_t1otsu = zeros(length(thresh_array), 1);
auc_t1kmeans = zeros(length(thresh_array), 1);
auc_t1gmm = zeros(length(thresh_array), 1);
for i = 1:length(thresh_array)
    [X, Y, T, AUC] = perfcurve(LGE_labels(:, i), num(:,2), 1);
    auc_t15sd(i) = AUC;
    [X, Y, T, AUC] = perfcurve(LGE_labels(:, i), num(:,3), 1);
    auc_t1otsu(i) = AUC;
    [X, Y, T, AUC] = perfcurve(LGE_labels(:, i), num(:,4), 1);
    auc_t1kmeans(i) = AUC;
    [X, Y, T, AUC] = perfcurve(LGE_labels(:, i), num(:,5), 1);
    auc_t1gmm(i) = AUC;
end

figure();
plot(thresh_array, auc_t15sd)
hold on;
plot(thresh_array, auc_t1otsu)
plot(thresh_array, auc_t1kmeans)
plot(thresh_array, auc_t1gmm)
grid on;
legend({'Mean+5SD', 'Otsu', 'Kmeans', 'GMM'}, 'Location', 'SouthEast');
xlabel('Thresholding');
ylabel('AUC');

auc_t15sd(16)
auc_t1otsu(16)
auc_t1kmeans(16)
auc_t1gmm(16)

figure();
plot(X,Y)
xlabel('False Positive Rate');
ylabel('True Positive Rate');

%% Algorithms
[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,2), 1);
AUC
figure();
plot(X,Y)
hold on;

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,3), 1);
AUC
plot(X,Y)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,4), 1);
AUC
plot(X,Y)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,5), 1);
AUC
plot(X,Y)

xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;
legend({'Mean+5SD', 'Otsu', 'Kmeans', 'GMM'}, 'Location', 'SouthEast');

%% 5SD
[X, Y, T, AUC] = perfcurve(LGE_labels(:, 2), num(:,2), 1);
AUC
figure();
plot(X,Y, 'LineWidth', 1.5)
hold on;

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 11), num(:,2), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,2), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 31), num(:,2), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 51), num(:,2), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;
legend({'1%','10%', '15%', '30%', '50%'}, 'Location', 'SouthEast');
title('T1 Mean+5SD');

%% Otsu
[X, Y, T, AUC] = perfcurve(LGE_labels(:, 2), num(:,3), 1);
AUC
figure();
plot(X,Y, 'LineWidth', 1.5)
hold on;

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 11), num(:,3), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,3), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 31), num(:,3), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 51), num(:,3), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;
legend({'1%','10%', '15%', '30%', '50%'}, 'Location', 'SouthEast');
title('T1 Otsu');

%% Kmeans
[X, Y, T, AUC] = perfcurve(LGE_labels(:, 2), num(:,4), 1);
AUC
figure();
plot(X,Y, 'LineWidth', 1.5)
hold on;

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 11), num(:,4), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,4), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 31), num(:,4), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 51), num(:,4), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;
legend({'1%','10%', '15%', '30%', '50%'}, 'Location', 'SouthEast');
title('T1 Kmeans');

%% GMM
[X, Y, T, AUC] = perfcurve(LGE_labels(:, 2), num(:,5), 1);
AUC
figure();
plot(X,Y, 'LineWidth', 1.5)
hold on;

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 11), num(:,5), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 16), num(:,5), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 31), num(:,5), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

[X, Y, T, AUC] = perfcurve(LGE_labels(:, 51), num(:,5), 1);
AUC
plot(X,Y,'LineWidth', 1.5)

xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;
legend({'1%','10%', '15%', '30%', '50%'}, 'Location', 'SouthEast');
title('T1 GMM');
%% ROC of two continuous variables
thresh_array_complete = 0:1:100;
LGE_labels_complete = zeros(size(num, 1), length(thresh_array_complete));
for i = 1:length(thresh_array)
   LGE_labels_complete(:, i) = num(:,1) > thresh_array(i); 
end

TPR = zeros(1,length(thresh_array_complete));
FPR = zeros(1, length(thresh_array_complete));

for i = 1:length(thresh_array_complete)
    T1_labels = num(:,2) > thresh_array_complete(i);
    TP = sum(and(LGE_labels_complete(:,i), T1_labels));
    TN = sum(and(~LGE_labels_complete(:,i), ~T1_labels));
    FP = sum(and(~LGE_labels_complete(:,i), T1_labels));
    FN = sum(and(LGE_labels_complete(:,i), ~T1_labels));
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
    TPR(i) = Sensitivity;
    FPR(i) = 1 - Specificity;
end


TPR(isnan(TPR)) = 0;
TPR_mod = [1 TPR];
FPR_mod = [1 FPR];
figure();
plot(FPR_mod, TPR_mod)
xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;
xlim([0 1]);
ylim([0 1]);
hold on;

TPR = zeros(1,length(thresh_array_complete));
FPR = zeros(1, length(thresh_array_complete));

for i = 1:length(thresh_array_complete)
    T1_labels = num(:,3) > thresh_array_complete(i);
    TP = sum(and(LGE_labels_complete(:,i), T1_labels));
    TN = sum(and(~LGE_labels_complete(:,i), ~T1_labels));
    FP = sum(and(~LGE_labels_complete(:,i), T1_labels));
    FN = sum(and(LGE_labels_complete(:,i), ~T1_labels));
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
    TPR(i) = Sensitivity;
    FPR(i) = 1 - Specificity;
end
TPR(isnan(TPR)) = 0;
TPR_mod = [1 TPR];
FPR_mod = [1 FPR];

plot(FPR_mod, TPR_mod)
TPR = zeros(1,length(thresh_array_complete));
FPR = zeros(1, length(thresh_array_complete));

for i = 1:length(thresh_array_complete)
    T1_labels = num(:,4) > thresh_array_complete(i);
    TP = sum(and(LGE_labels_complete(:,i), T1_labels));
    TN = sum(and(~LGE_labels_complete(:,i), ~T1_labels));
    FP = sum(and(~LGE_labels_complete(:,i), T1_labels));
    FN = sum(and(LGE_labels_complete(:,i), ~T1_labels));
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
    TPR(i) = Sensitivity;
    FPR(i) = 1 - Specificity;
end
TPR(isnan(TPR)) = 0;
TPR_mod = [1 TPR];
FPR_mod = [1 FPR];

plot(FPR_mod, TPR_mod)

TPR = zeros(1,length(thresh_array_complete));
FPR = zeros(1, length(thresh_array_complete));

for i = 1:length(thresh_array_complete)
    T1_labels = num(:,5) > thresh_array_complete(i);
    TP = sum(and(LGE_labels_complete(:,i), T1_labels));
    TN = sum(and(~LGE_labels_complete(:,i), ~T1_labels));
    FP = sum(and(~LGE_labels_complete(:,i), T1_labels));
    FN = sum(and(LGE_labels_complete(:,i), ~T1_labels));
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
    TPR(i) = Sensitivity;
    FPR(i) = 1 - Specificity;
end
TPR(isnan(TPR)) = 0;
TPR_mod = [1 TPR];
FPR_mod = [1 FPR];

plot(FPR_mod, TPR_mod)
legend({'Mean+5SD', 'Otsu', 'Kmeans', 'GMM'}, 'Location', 'SouthEast');