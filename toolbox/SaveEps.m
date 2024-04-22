function SaveEps(dir,name)

mymkdir(dir);
saveas(gcf,[dir,filesep,name,'.eps'],'epsc'); %vector 常用；
%saveas(gcf,[dir,filesep,name,'.emf'],'meta');
saveas(gcf,[dir,filesep,name,'.bmp'],'bmp'); %bitmap
saveas(gcf,[dir,filesep,name,'.fig'],'fig');
saveas(gcf,[dir,filesep,name,'.tif'],'tiff'); %常用 bitmap  不会因为编译改变
saveas(gcf,[dir,filesep,name,'.pdf'],'pdf'); %vector

end