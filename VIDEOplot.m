figure(2)

    set(gcf,'Color','w')
   s1 = 't = ';
   s3 = ' s';


v = VideoWriter('Iarray_Isc_400mm15turntry2.mp4','MPEG-4');
v.FrameRate = 2;
open(v);

for k = 1:size(IArrayhistory,2)
    % if tArrayhistory(k) > 11
    %     break
    % end


   titlestring = strcat(s1,num2str(tArrayhistory(k),3),s3);

   plot(IArrayhistory(:,k))
   hold on
   plot(IcArrayhistory(:,k))
   hold off
   title(titlestring)
   
   legend('I in line element','I_c of line element',Location='southwest')
   grid on
      ylim([-4000 8000])
   xlim([1 length(IcArray)])
ylabel('Current [A]')
   xlabel('Line element number [-]')
   frame = getframe(gcf);
   writeVideo(v,frame);

end

close(v);















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