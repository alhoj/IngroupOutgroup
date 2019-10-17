function cropAndAddAlpha(inputDir) 
%function cropAndAddAlpha(inputDir) 
%----------------------------------
%Takes input images with a black background and crops them to fit the
%contents and makes the background transparent.
%Input directory (inputDir variable) should be a directory that contains
%only images you want to crop and convert the background to transparent.
%Non-image files will be skipped but all the jpg and png files will be
%overwritten with png-files (jpg files will be converted to png files and
%the originals are removed in the process).
%
%16.11.2012 Juha Lahnakoski
%juha.lahnakoski@aalto.fi

curDir=pwd;
cd(inputDir);
d=dir;
d=d(find(~cellfun(@sum,{d(:).isdir})));

for k=1:length(d)
    %If the files are in png format we can use them as they are
    if strcmp('png',d(k).name(end-2:end))
        %Read the input file
        img=imread(d(k).name);
        %Find the black voxels
        [x,y]=find(max(img,[],3)>10);
        %Write the file with alpha channel
        imwrite(img(min(x):max(x),min(y):max(y),:),d(k).name,'Alpha',double(max(img(min(x):max(x),min(y):max(y),:),[],3)>10));
        
        %If the files are jpegs we can just change the name and resave the
        %png-files on top of the old ones. (Other file formats could be
        %added)
    elseif strcmp('jpg',d(k).name(end-2:end))
        %Rename the jpeg to png and read
        oldname=d(k).name;
        d(k).name=sprintf('%spng',d(k).name(1:end-3));
        movefile(oldname,d(k).name)
        img=imread(d(k).name);
        %Find the black voxels
        [x,y]=find(max(img,[],3)>10);
        %Write the file with alpha channel
        imwrite(img(min(x):max(x),min(y):max(y),:),d(k).name,'Alpha',double(max(img(min(x):max(x),min(y):max(y),:),[],3)>10));
    end;
        
    
end;

cd(curDir);

end