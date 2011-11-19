function S = my_parseXML(filename)
% PARSEXML converts an XML file to a MATLAB structure. This XML parser will
% only store elements, attributes, text, comments, and CDATA in a struct.
% Note that indexing the XML object starts at 0, not 1.
% input: filename = the name of the XML file to be parsed.
% output: S = the contents of the XML file are converted into a
% struct.
%
%Note Java Objects in MATLAB do not share the same memory limits as
%traditional MATLAB variables. Therefore, my_parseXML may choke on files
%larger than 64MB, the assumed default limit.

try
    xdoc = xmlread(filename);
catch err
    disp(err.message)
    error('Failed to read XML file %s.',filename);
end

%Initialize S
root_node = xdoc.getDocumentElement; %Locates the node at the root of the XML tag tree
my_temp_name = regexprep(root_node.getNodeName.toString.toCharArray','[-:.]','_'); %Struct names cannot contain [-:.]
my_tree = ['S.' my_temp_name]; %The location in S where data is currently being created
if root_node.hasAttributes
    my_attributes = root_node.getAttributes;
    my_length = my_attributes.getLength;
    for i = 0:(my_length-1)
        my_temp_string = my_attributes.item(i).toString.toCharArray';
        my_temp_name = regexp(my_temp_string,'.*(?==)','match');
        my_temp_name = regexprep(my_temp_name,'[-:.]','_');
        my_temp_attribute = regexp(my_temp_string,'(?<=").*(?=")','match');
        eval([my_tree '.attd8a.' my_temp_name{1} '=''' my_temp_attribute{1} ''';'])
    end
end

if root_node.hasChildNodes
    current_node = root_node.getFirstChild;
    while ~isempty(current_node)
        switch current_node.getNodeType
            case 1
                %NodeType 1 = element
                S  = parseElement(current_node,S,my_tree);
            case {3,4,8}
                %NodeType 3 = text node
                %NodeType 4 = CDATA
                %NodeType 8 = comment node
                S = parseMiscellaneous(current_node,S,my_tree);
            otherwise
        end
        current_node = current_node.getNextSibling;
    end
end
end

% ----- Subfunction PARSEATTRIBUTES -----
function S = parseAttributes(node,S,tree)
my_attributes = node.getAttributes;
my_length = my_attributes.getLength;
for i = 0:(my_length-1)
    my_temp_string = my_attributes.item(i).toString.toCharArray';
    my_temp_name = regexp(my_temp_string,'.*(?==)','match');
    my_temp_name = regexprep(my_temp_name,'[-:.]','_');
    my_temp_attribute = regexp(my_temp_string,'(?<=").*(?=")','match');
    eval([tree '.attd8a.' my_temp_name{1} '=''' my_temp_attribute{1} ''';']);
end
end

% ----- Subfunction PARSEELEMENT -----
function S = parseElement(node,S,tree)
%WARNING: This is a recursive function and the author of this code does not
%know if this could lead to an infinite loop or crashing MATLAB.
my_temp_name = regexprep(node.getNodeName.toString.toCharArray','[-:.]','_'); %Struct names cannot contain [-:.]
if eval(['isfield(' tree ',''' my_temp_name ''')'])
    my_temp_tree = [tree '.' my_temp_name '(end+1)']; %The location in S where data is currently being created
    eval([my_temp_tree '.txtd8a={};'])
    tree = [tree '.' my_temp_name '(end)']; %The location in S where data is currently being created
else
    tree = [tree '.' my_temp_name '(1)']; %The location in S where data is currently being created
    eval([tree '.txtd8a={};'])
end
if node.hasAttributes
    S = parseAttributes(node,S,tree);
end
if node.hasChildNodes
    current_node = node.getFirstChild;
    while ~isempty(current_node)
        switch current_node.getNodeType
            case 1
                %NodeType 1 = element
                S  = parseElement(current_node,S,tree);
            case {3,4,8}
                %NodeType 3 = text node
                %NodeType 4 = CDATA
                %NodeType 8 = comment node
                S = parseMiscellaneous(current_node,S,tree);
            otherwise
        end
        current_node = current_node.getNextSibling;
    end
else
    my_temp_string = node.getTextContent.toString.toCharArray';
    flag = checkString(my_temp_string);
    if flag
        eval([tree '.txtd8a{1}=my_temp_string;'])
    end
end
end

% ----- Subfunction PARSEMISCELLANEOUS -----
function S = parseMiscellaneous(node,S,tree)
my_temp_string = node.getData.toString.toCharArray';
flag = checkString(my_temp_string);
if flag
    if eval(['isfield(' tree ',''txtd8a'')'])
        eval([tree '.txtd8a{end+1}=my_temp_string;'])
    else
        eval([tree '.txtd8a{1}=my_temp_string;'])
    end
end
end

% ----- Subfunction CHECKSTRING -----
function [flag] = checkString(my_str)
flag = 1;
if strcmp(my_str,sprintf('\n'))
    flag = 0;
end
if isempty(regexp(my_str,'\S','once'))
    flag = 0;
end
end
