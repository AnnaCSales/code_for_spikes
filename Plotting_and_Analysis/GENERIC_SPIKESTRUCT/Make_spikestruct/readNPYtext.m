function [messages] = readNPYtext(filename)
%A very rough and ready function to read in an NPY array containing
%characters...as returned from openEphys recordings in binary
%format.
% Expects to be passed a string containing a filename (with path)
% Returns a string array, where each array element is the text
% content from one message. 
  
    fid = fopen(filename);
    contents = fread(fid);
    contents(contents==0)=[];
    contents=contents';
    message_all=string(char(contents));
    messages=strsplit(message_all, ',');
    messages(1:4)=[];
    fclose(fid);
end

