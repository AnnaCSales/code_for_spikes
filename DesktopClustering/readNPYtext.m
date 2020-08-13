

filename='text.npy';
fid = fopen(filename);
contents = fread(fid);
contents(contents==0)=[];
contents=contents';
message_all=string(char(contents));
message_parts=strsplit(message_all, ',');
message_parts(1:5)=[];
