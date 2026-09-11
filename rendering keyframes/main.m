format compact
clear
clc
clf reset

% ----------

fps = 30;

%1. frame
%2. launch TA
%3. eject v
%4. azimuth
%5. elevation
%6. center x
%7. center y
%8. width x/2
%9. 1440 resolution div factor
%10. time

keyframe_raw = [
1, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1, 3600*24*135
5, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
100, 0.25, 3e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
125, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
400, 0.25, 2.35e3, 90+360, 0, 0, 0, 180, 1.5, 3600*24*135
401, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
600, 0.25, 2.35e3, 90, 90, 0, 0, 180, 1.5, 3600*24*135
625, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
800, 0.25+1, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
801, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*135
830, 0.25, 2.35e3, 90, 0, 0, 0, 180, 1.5, 3600*24*1e3
];
keyframe_raw = single(keyframe_raw);

for n=1:height(keyframe_raw)-1
    if abs(keyframe_raw(n,1) - keyframe_raw(n+1,1)) <= 2
        keyframe_raw(n+1,1) = keyframe_raw(n+1,1)+2;
    end
end

keyframes_duplicate = keyframe_raw;
keyframes_tmp = keyframe_raw;
for n=1:2
    keyframes_duplicate = [keyframes_duplicate; [keyframes_tmp(:,1)+n , keyframes_tmp(:,2:end) ] ];
end
[~,ind_sort] = sort(keyframes_duplicate(:,1));
keyframes_duplicate = keyframes_duplicate(ind_sort,:);

frames = keyframes_duplicate(1,1):keyframes_duplicate(end,1);

keyframes_interp(:,1) = frames;

for n=2:width(keyframes_duplicate)
    keyframes_interp(:,n) = interp1(keyframes_duplicate(:,1),keyframes_duplicate(:,n),keyframes_interp(:,1),"pchip");
    keyframes_interp(:,n) = round(keyframes_interp(:,n),5);
end
keyframes_interp(:,2) = rem(keyframes_interp(:,2),1);
keyframes_interp(:,4) = rem(keyframes_interp(:,4),360);

start_frame = 1;

for ind_t = start_frame:height(keyframes_interp)

    pass_matrix = [];
    for m=2:width(keyframes_interp)
        pass_matrix(m-1) = keyframes_interp(ind_t,m);
    end
    build_fractal_image(pass_matrix);

    frame = [];
    frame_new = uint8([]);
    frame = getframe(gcf);
    for n=1:3
        frame_new(:,:,n) = uint8(imresize(squeeze(frame.cdata(:,:,n)), [nan, 1920*2 ]));
    end

    fprintf("rendered frame %i.\n",ind_t)

    imwrite(frame_new, "frame_"+string(ind_t)+".png");

end

for n=1:100
    pause(1);
    sound(sin(2*pi*400*(0:1/14400:0.3)), 14400);
end




