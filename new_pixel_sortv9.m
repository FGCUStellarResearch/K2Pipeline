format long
clear all
close all

%folder = input('Assuming you are Carly on Colossus. Is this correct? (Y/N) ','s');
%if (folder == 'n') || (folder == 'N')
%    folder = input('Please specify correct folder name/path: ','s');
%else
   folder = '/home/derek/K2_PHunt/Data';
   folder = '/media/derek/TOSHIBA/K2_data/K2_456/test';
   folder = '/media/derek/TOSHIBA/K2_data/Dhital/test';
   %folder = '/media/derek/TOSHIBA/temp';
   %folder = '/media/derek/TOSHIBA/temp/temp2';
   %%folder = '/media/derek/TOSHIBA/K2_data/PH';
   %folder = '/media/derek/TOSHIBA/Clayton';
   folder = '/media/derek/TOSHIBA/K2_data/Stello/C8/raw';
   folder = '/media/derek/TOSHIBA/K2_data/Stello/C5';
   folder = '/media/derek/TOSHIBA/K2_data/PH/2017';
   folder = '/media/derek/TOSHIBA/K2_data/WideBinaries';
   
   %folder = '/media/derek/TOSHIBA/K2_data/Trappist';
   folder = '/media/derek/TOSHIBA/K2_data/PH/2017/C10';
%end

    % uncommented for loop allows multiple files to be run at once:
  %  file = input('full file name:  ','s');
  %  file = strcat(folder,'/',file); %

kplrfiles = dir(strcat(folder,'/ktwo*'));
   for fileNum = 1:(length(kplrfiles))
     %kplrfiles = dir(folder)
   %for fileNum = 3:(length(kplrfiles))
         
         file = kplrfiles(fileNum).name;
         file = strcat(folder,'/',file);
     %   disp('getting data')
     clear time flux bkg xcen ycen
     clear mask
        [time,apdim,n,dim,series,data,cnum] = get_k2_data(file);
        
        series(isnan(series))= 0;
      %  meanseries = mean(series(:))
      %  disp('calling mean_image_gen')
        sat_flag = 0;
        [mean_image] = mean_image_gen(series,dim,n,file);
        tmax = max(mean_image(:));
       % if max(mean_image(:))>1.8e5
       %     sat_flag = 1;
       % end
        
        %Need to get smarter about how we handle saturated frames, to deal
        %with the case where there's a saturated and a nonsaturated star in
        %the same frame. I think the code handles it OK, the problem is
        %identifying that a saturated or near-saturated cluster is really only ONE target and not
        %multiples...
        
      %  disp('back from mean_image_gen')
        filename = strrep(file,'_lpd-targ.fits','');
        filename = strrep(filename,'ktwo','outputs/pipeout_ktwo');
        filename = strcat(filename,'_','mean_image');
        dlmwrite(filename, mean_image);
        
        %Characterize and remove background
        bkg = fit_background2(apdim,series);
        
      %  disp('Background removed')
        
        
        for count = 1:apdim
            series(:,:,count) = series(:,:,count) - bkg(count);
            temp = series(:,:,count);
            if min(temp(:))<0
                series(:,:,count) = series(:,:,count)-min(temp(:));
            end
        end
    
        %sat_flag = 1;
        if sat_flag>0
        
        %First start to handle saturated frames
        test_image = zeros(dim(1),dim(2));
        test_image(mean_image>3*mean(bkg))=1;
        cc = bwconncomp(test_image);
        numPixels = cellfun(@numel,cc.PixelIdxList);
        [biggest,idx] = max(numPixels);
        mask{1} = cc.PixelIdxList{idx};
        num_mask = 1;
        num = num_mask;
        for count = 1:apdim
            temp = series(:,:,count);
            temp = temp(:);
            test_image(~cc.PixelIdxList{idx}) = 0;
            temp = temp.*test_image(:);
            flux(1,count) = sum(temp(:));
        end
        
     
        
        %adjust mask-finding routine to be water-level type....but start at
        %level where there's only one? or just a fractin of the peak value?
        % and LOWER water level until there are >1 or until 3-bkg level
        % is reached
        
        %%OK probably should update variable names above
        %%next extend to non-saturated case
        
        else
           % temp_image = mean_image;
            temp_image = zeros(dim(1),dim(2));
            cutoff = max(5*mean(bkg),2e-3*max(mean_image(:)));
            temp_image(mean_image>cutoff)=1;
            temp_image = mean_image;
            cutoff = min(3*mean(bkg),max(mean_image(:))-1);
            temp_image(mean_image<cutoff)=0; %undeleted this one
            if max(mean_image(:))>1.8e5
                temp_image = imgaussfilt(temp_image,1);
            else
                temp_image = imgaussfilt(temp_image);
            end
            bw = imregionalmax(temp_image);
            [out,num] = bwlabel(bw);
            
            %now we've determined the number num of maxima in the image
            %remove the maxes that are too faint...
            
           % stop
           
           %%still need to fix what happens when num = 1, and also do we
           %%still need the saturated kluge?
            
           
           if num>1
           tmp_num = num;
           for i=1:num
               tmp_out = out;
               tmp_out(tmp_out~=i) = 0;
               tmp_out = tmp_out/i;
               max_test = mean_image.*tmp_out;
               if max(max_test(:)) < 50 && tmp_num > 1 %is this how we want to set num?
                   tmp_num = tmp_num-1;
                   out(out==i) = 0;
               end
           end
           num = tmp_num;
           end
           
           %% at this point the number of targets in the frame, num,  has been determined
           
           id_vals = unique(out);
           if id_vals(1) == 0
               id_vals(1) = [];
           end
           clear peak_val
           for i = 1:numel(id_vals)
               peak_val(i) = mean_image(out == id_vals(i)); %%these are the values at the peaks
           end
            
           %% at this point the values of the peaks for each target have been determined
           
           %% adding a section to handle situation where there is only one mask
           
           if num == 1
               %divide range between peak and mean background level into
               %100 pieces
               min_val = ceil(mean(bkg(:)));
               %min_val = 1000;
               max_val = floor(peak_val(1))-1;
               range_val = 0.01*(max_val-min_val);
               %for each value need to calculate the relevant statistic...
               for calc_val = 1:100
                    t_image = mean_image;
                    t_image(t_image<(min_val+calc_val*range_val)) = 0; %zero all values below the reference value
                    cct = bwconncomp(t_image,4);
                    mask{1} = cct.PixelIdxList{1};
                    clear tflux
                    for count = 1:apdim
                        temp = series(:,:,count);
                        temp = temp(mask{1});
                        tflux(count) = sum(temp(:));
                    end
                    %tstat(calc_val) = sum(abs(diff(diff(diff(tflux)))));
                    tstat(calc_val) = sum(abs(diff(tflux)));
                    tstat(calc_val) = tstat(calc_val)/mean(tflux);
                    %plot(tflux/mean(tflux))
                    %calc_val, std(tflux)/mean(tflux)
                    %hold on
               end
               
               [aa bb] = min(tstat);
               t_image = mean_image;
               t_image(t_image<(min_val+bb*range_val)) = 0;
               
               %bb=1; %this approach doesn't work for saturated targets!
               %%look at brightest 9 pixel group mean
                sort_pix = sort(mean_image(:),'descend');
                sort_pix = mean(sort_pix(1:9));
                if sort_pix > 1.7e5
                    bb = 1;
                end
               %t_image(t_image<(min_val+bb*range_val)) = 0;
               %t_image(t_image<450) = 0;
               cct = bwconncomp(t_image);
               %mask{1} = cct.PixelIdxList{2};
               mask{1} = cct.PixelIdxList{1};
               %numel(mask{1})
           else
           
           temp = find(out);
           fin_ref_value = zeros(1,num);

 for ii=1:num
        ploc = temp(ii); % location of this peak
        parr = setdiff(temp,ploc); %all peaks not this one
        ref_value = ceil(peak_val(ii))-1; %start at value at this peak
        %fin_ref_value(ii) = 1;
         %check if current peak is inside
         for j = ref_value:-1:0 %Good for step 1, but will have to replace with binary search tree 
            t_image = mean_image;
            t_image(t_image<j) = 0; %zero all values below the reference value
            cct = bwconncomp(t_image);
            for i=1:cct.NumObjects
                if ismember(ploc,cct.PixelIdxList{i}) & any(ismember(parr,cct.PixelIdxList{i})) 
                    fin_ref_value(ii) = j+1; %overlapping masks...increase ref value and solution is found
                    
                    %cct.PixelIdxList{ii} %Need to figure out how to label
                    %the elements in the mask!
                end
            end
            if fin_ref_value(ii) > 0 %if this has been set, then we are done!
                t_image = mean_image;
                t_image(t_image<fin_ref_value(ii)) = 0;
                cct = bwconncomp(t_image,4);
                for k = 1:cct.NumObjects
                   % disp(ploc)
                    if ismember(ploc,cct.PixelIdxList{k})
                      %  cct.PixelIdxList{1};
                      %  cct.PixelIdxList{2};
                        mask{ii} = cct.PixelIdxList{k};
                    end
                end
                break
            end
         end
 end
 num = numel(mask);
 
 %%I think in here is where we need to check mask sizes for the multi-mask
 %%case, to see if shrinking the mask improves the SNR of the dewhitened
 %%time series. Should be able to steal code from single-mask case!
 
 %for ii=1:num
 %     %divide range between peak and mean background level into
 %              %100 pieces
 %              min_val = ceil(fin_ref_value(ii));
 %              %min_val = 1000;
 %              max_val = floor(peak_val(ii))-1;
 %              range_val = 0.01*(max_val-min_val);
 %              %for each value need to calculate the relevant statistic...
 %              for calc_val = 1:100
 %                   t_image = mean_image;
 %                   t_image(t_image<(min_val+calc_val*range_val)) = 0; %zero all values below the reference value
 %                   cct = bwconncomp(t_image);
 %                   mask{ii} = cct.PixelIdxList{ii};
 %                   for count = 1:apdim
 %                       temp = series(:,:,count);
 %                       temp = temp(mask{ii});
 %                       tflux(count) = sum(temp(:));
 %                   end
 %                   tstat(calc_val) = sum(abs(diff(diff(diff(tflux)))));
 %                   tstat(calc_val) = tstat(calc_val)/mean(tflux);
 %              end
 %              
 %              [aa bb] = min(tstat);
 %              t_image = mean_image;
 %              t_image(t_image<(min_val+bb*range_val)) = 0;
 %              cct = bwconncomp(t_image);
 %              mask{ii} = cct.PixelIdxList{ii};
 %end
     
 
           end
  %% at this point the masks have been determined

            end
            

          %  end
            

  
  
  %Calculate fluxes and centroids
  
             for i=1:num
                 %prep for centroid calculations
                 %first convert to x-y positions
                 [xx yy] = ind2sub([dim(1), dim(2)],mask{i});
                for count = 1:apdim
                    temp = series(:,:,count);
                    temp = temp(mask{i});
                    flux(i,count) = sum(temp(:));
                    %calculate 
                    xflux = 0.0;
                    yflux = 0.0;
                    for j = 1:numel(mask{i})
                        xflux = xflux + xx(j)*series(xx(j),yy(j),count);
                        yflux = yflux + yy(j)*series(xx(j),yy(j),count);
                    end
                    xcen(i,count) = xflux/flux(i,count);
                    ycen(i,count) = yflux/flux(i,count);
                end
                xc(i) = mean(xcen(i,:));
                yc(i) = mean(ycen(i,:));
             end
            
                     
             
  %%temporarily try centroids which are simply the weighted sum of the entire background-subtracted frame
   for i=1:num
                 %prep for centroid calculations
                 %first convert to x-y positions
                 [xx yy] = ind2sub([dim(1), dim(2)],[1:1:dim(1)*dim(2)]);
               
                for count = 1:apdim
                    temp = series(:,:,count);
                   % temp = temp(mask{i});
                    tflux(i,count) = sum(temp(:));
                    %calculate 
                    xflux = 0.0;
                    yflux = 0.0;
                    for j = 1:(dim(1)*dim(2))
                        xflux = xflux + xx(j)*series(xx(j),yy(j),count);
                        yflux = yflux + yy(j)*series(xx(j),yy(j),count);
                    end
                    xcen(i,count) = xflux/tflux(i,count);
                    ycen(i,count) = yflux/tflux(i,count);
                end
                xcen(i,:) = xcen(i,:)-mean(xcen(i,:))+xc(i);
                ycen(i,:) = ycen(i,:)-mean(ycen(i,:))+yc(i);
   end
  
   %reorder based on mean flux!
  % reorganize mask order so 1 is brightest
  clear nmask
  clear nmean
  nmask = mask;
  nflux = flux;
  nxcen = xcen;
  nycen = ycen;
  if num > 1
      for i=1:num
          nmean(i) = mean(flux(i,:));
      end
      [ai bi] = sort(nmean,'descend');
      for i=1:num
          mask{i} = nmask{bi(i)};
          flux(i,:) = nflux(bi(i),:);
          xcen(i,:) = nxcen(bi(i),:);
          ycen(i,:) = nycen(bi(i),:);
      end
  end
   
  %put raw lightcurves out to file
            out_to_file2(time,flux,bkg,xcen,ycen,file,num);    
            
            
            % export image of mask
            mask_image = zeros(dim(1),dim(2));
            for i=1:num
                mask_image(mask{i}) = i;
            end
            hmask = imagesc(mask_image);
            colorbar;
            
            filename = strrep(file,'_lpd-targ.fits','');
            filename = strrep(filename,'ktwo','outputs/pipeout_ktwo');
            filename = strcat(filename,'_','mask.jpg');
            saveas(hmask,filename);
            
   end
   
