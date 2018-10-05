clear all;
close all;

camera = webcam;
net = alexnet;

%%
while true
    im = snapshot(camera);
    image(im);
    im = imresize(im, [227 227]);
    label = classify(net,im);
    title(char(label));
    drawnow
end