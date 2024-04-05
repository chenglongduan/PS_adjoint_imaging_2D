function source_file(sxpos, szpos, outfilename)
% sxpos, szpos are vectors


% convert to column vector
if (size(sxpos,1)==1)
    sxpos = sxpos';
end
if (size(szpos,1)==1)
    szpos = szpos';
end

source_xz = [sxpos, szpos];

nsrc = size(source_xz,1);

fileID = fopen(outfilename,'w');

fprintf(fileID,'%d\n',nsrc);
fprintf(fileID,'%.2f %.2f\n',source_xz');

fclose(fileID);

end