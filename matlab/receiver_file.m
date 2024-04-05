function receiver_file(isrc, rxpos, rzpos, outfile_directory)
% rxpos, rzpos are vectors


outfilename = strcat(outfile_directory,'rec_shot',num2str(isrc),'.txt');

% convert to column vector
if (size(rxpos,1)==1)
    rxpos = rxpos';
end
if (size(rzpos,1)==1)
    rzpos = rzpos';
end

receiver_xz = [rxpos, rzpos];

fileID = fopen(outfilename,'w');

fprintf(fileID,'%.2f %.2f\n',receiver_xz');

fclose(fileID);

end
