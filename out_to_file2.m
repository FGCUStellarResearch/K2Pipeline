function out_to_file2(time,flux,bkg,xcentroid,ycentroid,file,targetNum)
    %whos('time')
    %whos('flux')
    %whos('bkg')
    %whos('xcentroid')
    %whos('ycentroid')
    
   for i=1:targetNum
        output_data = [time'; flux(i,:) ;bkg; xcentroid(i,:); ycentroid(i,:);];
     % places all output data into same matrix to print
    
    
        filename = strrep(file,'_lpd-targ.fits','');
        filename = strrep(filename,'ktwo','outputs/pipeout_ktwo');
        filename = strcat(filename,'_target', num2str(i));
        % creates output filename based on input filename, preserving target ID


        ID = fopen(filename,'w');
    
        fprintf(ID,'%12.8f %f %f %f %f\n',output_data);
        fclose(ID);
        % will overwrite any existing target output file
    end
end
