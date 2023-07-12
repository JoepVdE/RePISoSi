figure(2)

    set(gcf,'Color','w')
   s1 = 't = ';
   s3 = ' s';


% v = VideoWriter('Vground5.mp4','MPEG-4');
% v.FrameRate = 1;
% open(v);
% 
% for k = 1:size(VGroundVectorhistory,2)
%     if tArrayhistory(k) > 40
%         break
%     end
% 
% 
%    titlestring = strcat(s1,num2str(tArrayhistory(k),3),s3);
% 
%    plot(1e6.*(VGroundVectorhistory(:,k))./(5.*2.*pi.*0.12));
% 
%    title(titlestring)
% 
%    legend('E [\muV/m]',Location='southwest')
%    grid on
%    %    ylim([-4000 8000])
%    % xlim([1 450])
% ylabel('E [\muV/m]')
%    xlabel('line element nr [-]')
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% 
% end
% 
% close(v);




index1 = round(length(VGroundVector)*1/10);
index2 = round(length(VGroundVector)*9/10);


for k = 1:size(VGroundVectorhistory,2)-1

Ecoil(k) = 1e6*(mean(VGroundVectorhistory(1:round(length(VGroundVector)*1/5),k)) - mean(VGroundVectorhistory(round(length(VGroundVector)*4/5):end,k)))/(2*pi*0.12);
end
figure(4)
plot(tArrayhistory(1:length(Ecoil)),Ecoil)
xlabel('time [s]')
ylabel('E [\muV/m]')
grid on
figure(5)
plot(tArrayhistory(1:length(centerBhistory)),centerBhistory)
xlabel('time [s]')
ylabel('B [T]')
grid on
figure(6)
plot(tArrayhistory(1:length(temperatureRingBottomhistory)),temperatureRingBottomhistory)
hold on
plot(tArrayhistory(1:length(temperatureRingBottomhistory)),temperatureRingTophistory)
grid on
xlabel('time [s]')
ylabel('T [K]')
legend('T_{bottomlead}','T_{toplead}')


% N = 31    ;
% 
% for i = 1:N
%     figure(1)  
%   %  imshow(processo(:,:,1,i))
% 
% plot(IArrayhistory(:,N))
% 
%       F(i) = getframe(gcf) ;
%       drawnow
%     end
%   % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo2.avi');
%   writerObj.FrameRate = 1;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);