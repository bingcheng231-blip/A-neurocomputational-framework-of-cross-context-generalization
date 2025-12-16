function [outputArg1] = s001_showRDM(aps, it4)
%S001_SHOWRDM 此处显示有关此函数的摘要
%   此处显示详细说明
    map_corr=[];
    for i=1:size(aps,1)
        for j=1:size(aps,2)
            map_corr(i,j,1:3)=it4(max(min(round(255*aps(i,j))+1,256),1),1:3);
        end
    end
    f = figure;
    f.Position(3:4) = [500 500];
    imshow(map_corr, 'border', 'tight', 'InitialMagnification','fit');
    outputArg1 = 'Done';
end

