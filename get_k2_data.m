function [time,apdim,n,dim,series,data,cnum] = get_k2_data(file)

    %file = input('EPIC number: ','s');
    
    fitsdata = fitsread(file,'binarytable');
    time = fitsdata{1}; % data time stamp
    data = fitsdata{5}; % calibrated pixel flux
    qual = fitsdata{10}; % quality flag
    cnum = fitsdata{3}; % cadence number, for dealing w/ campaign 0's safe mode
 
    % remove those points associated with thruster firings - Campaign 9
    % only!!!
   % test = csvread('/media/derek/TOSHIBA/Clayton/c9_thruster_firing.csv');
   % flag1 = test(:,2);
   % flag2 = test(:,3);
   % cadence = test(:,1);
   % rem_cadence = cadence(flag1>0 | flag2>0);
   % [C ia ib] = intersect(cnum,rem_cadence);
   % rem_rows = [ia]';
   
    %time(ia) = [];
    %data(ia) = [];
    %qual(ia) = [];
    %cnum(ia) = [];
 
    
%whos('rem_cadence')
%whos('qual')
    
    [apdim,n] = size(data); % apdim -> number of images
    dim = fitsinfo(file);
    dim = dim.Image.Size; % dim(1) x dim(2) -> size of total image


    % remove images where the quality is flagged non-zero or is bad in general
    % Remember: Changed to remove all data from before safe mode! C#0 sp.
    % check fits file for campaign 0
    %ind = find(cnum==89347); % get the index of the first good cadence
    %rem_rows = [1:ind]; % indicies of rows of 'bad' data
    rem_rows = [];
    

    for i = apdim:-1:1 %was ind
        if (qual(i) > 0)
            rem_rows = [rem_rows,i];
            % check for any images with quality tags
       % else if (max(data(i,:)) < 1000)
                rem_rows = [rem_rows,i];
                % check that the image max flux has expected magnitude
       %     end
        end
    end


    %[i,j] = find(data < 0);
    
    data(rem_rows,:)=[];
    %data(:,j)=[];
    time(rem_rows,:)=[];
    cnum(rem_rows,:)=[];
    [apdim,n] = size(data);    
    
    
    % because the matrix had 'bad' data removed
    % apdim will change, but n shouldn't. If it does, that's bad
    
    temp = dim(1,1);
    dim(1,1) = dim(1,2);
    dim(1,2) = temp;
    
    
    series = zeros(dim(1),dim(2),apdim);
    for i = 1:apdim
        % create rectangular images out of data
        series(:,:,i) = reshape(data(i,:),dim(1),dim(2));
    end
    
end