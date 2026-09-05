iVers = 1;
iVers = 2;
iVers = 3;

if iVers == 1
  oni = load(oni_data);
elseif iVers == 2  
  oni = readmatrix(oni_data, 'Delimiter', ' ', 'CommentStyle', '%');
  disp(size(oni));
elseif iVers == 3
  %Open the file
  fid = fopen(oni_data, 'r');
  
  % Read data, automatically skipping lines that start with '%'
  % '%f' will read both integers and real numbers into doubles
  % 'NumHeaderLines', 7 explicitly skips the first 7 lines
  data = textscan(fid, repmat('%f', 1, 13), 'Delimiter', ' ', 'CommentStyle', '%', 'NumHeaderLines', 7);
  
  % Close the file
  fclose(fid);
  
  % Convert the resulting cell array into a clean 26x13 matrix
  oni = cell2mat(data);
  oni = oni(1:27,:);
end  
