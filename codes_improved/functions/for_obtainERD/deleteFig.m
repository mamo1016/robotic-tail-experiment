function deleteFig()
handles = findall(groot,'Type','figure');
if ~isempty(handles)
    for i = 1:length(handles)
        try
            delete(handles);
            str = 'deleted';
        catch
            str = 'There are no processed windows';
        end
    end
    % disp(str)
end