[filename1 pathname1]=uigetfile('*.fig','Select a figure file to extract x-y data');
if isequal(filename1,0) | isequal(pathname1,0); return; end %User pressed Cancel
[filename2, pathname2] = uiputfile('*.xls', 'Provide an Excel file name to save the x-y data');
if isequal(filename2,0) | isequal(pathname2,0)
    return
end %User pressed Cancel
    
if isempty(strfind(filename2,'.xls'));
    filename2=[filename2 '.xls'];
end
s=hgload(fullfile(pathname1,filename1));   % open figure and get handle to figure
delete(legend);  %legends add extra data to output. Remove legends

%get line handles
h = findobj(s,'Type','line'); %line is the type of your figure file.
if isempty(h)
line{1}='Your figure file does not contain x-y data.';
    line{2}='Did you accidentally load a wrong file?';
line{3}='Your can view your figures using the File menu.';
line{4}='Make sure you do not load a wrong file.';
uiwait(msgbox(line, 'No x-y data in your file','warn'));
close(s)
    return
end
x=get(h,'xdata');
y=get(h,'ydata');
close(s) %close figure with legends removed
s=hgload(fullfile(pathname1,filename1));  %show original figure with legends