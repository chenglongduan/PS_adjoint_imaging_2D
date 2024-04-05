function rec_dat=extract_shot_gathers(data_mat, geom_tab, outsrc_filename, outrec_dir, outdat_dir, tshift, dt)
% convert to source gathers with corresponding geometry files from offset-aligned files


ntr_data = size(data_mat,2);
nt = size(data_mat,1);

% data time shift
if tshift~=0
    itshift = round(abs(tshift)/dt);
    data_tmp = data_mat;
    data_mat = zeros(nt,ntr_data);
    if(tshift>0)
        data_mat(itshift+1:nt,:) = data_tmp(1:nt-itshift,:);
    else
        data_mat(1:nt-itshift,:) = data_tmp(itshift+1:nt,:);
    end
end

% source
src_uniq = unique(geom_tab(:,6));
nsrc = length(src_uniq);
if ~isempty(outsrc_filename)
    source_file(src_uniq, zeros(nsrc,1), outsrc_filename);
end

% receiver, data
rec_dat = cell(nsrc,3);
for isrc=1:nsrc
    index = geom_tab(:,6)==src_uniq(isrc);
    rec_dat{isrc,1} = sum(index);
    rec_dat{isrc,2} = geom_tab(index,7);
    rec_dat{isrc,3} = data_mat(:,index);

    if ~isempty(outrec_dir)
        receiver_file(isrc, rec_dat{isrc,2}, zeros(rec_dat{isrc,1},1), outrec_dir);
    end

    if ~isempty(outdat_dir)
        dat_file = strcat(outdat_dir,'obs.bin.shot',num2str(isrc));
        fp=fopen(dat_file,'wb');fwrite(fp,rec_dat{isrc,3},'float32');fclose(fp);
    end
end

if (sum([rec_dat{:,1}]) ~= ntr_data)
    error('Err: check trace number in the source gather!');
end


end
