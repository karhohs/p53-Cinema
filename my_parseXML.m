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
my_temp_name = regexprep(root_node.getNodeName,'[-:.]','_'); %Struct names cannot contain [-:.]
my_tree = ['S.' my_temp_name]; %The location in S where data is currently being created
if root_node.hasAttributes
    my_attributes = root_node.getAttributes;
    my_length = my_attributes.getLength;
    for i = 0:(my_length-1)
        my_temp_string = my_attributes.item(i).toString.toCharArray';
        my_temp_name = regexp(my_temp_string,'.*(?==)','match');
        my_temp_name = regexprep(my_temp_name,'[-:.]','_');
        my_temp_attribute = regexp(my_temp_string,'(?<=").*(?=")','match');
        eval([my_tree '.attd8a.' my_temp_name '=' my_temp_attribute]);
    end
end

if root_node.hasChildNodes
    current_node = root_node.getFirstChild;
    element_names = cell(1,2);
    element_names{1,1} = ''; %A junk string to initialize the cell
    element_names{1,2} = 0;
    while ~isempty(current_node)
        switch current_node.getNodeType
            case 1
                %NodeType 1 = element
                my_temp_string = current_node.getNodeName.toString.toCharArray';
                [ind,element_names] = countElement(my_temp_string,element_names);
                S  = parseElement(current_node,S,ind,my_tree);
            case [3,4,8]
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
    eval([tree '.attd8a.' my_temp_name '=''' my_temp_attribute ''';']);
end
end

% ----- Subfunction COUNTELEMENT -----
function [ind,names] = countElement(str,names)

end

% ----- Subfunction PARSEELEMENT -----
function S = parseElement(node,S,ind,tree)
%WARNING: This is a recursive function and the author of this code does not
%know if this could lead to an infinite loop or crashing MATLAB.
my_temp_name = regexprep(node.getNodeName,'[-:.]','_'); %Struct names cannot contain [-:.]
tree = [tree '(' ind ').' my_temp_name]; %The location in S where data is currently being created
if node.hasAttributes
    S = parseAttributes(node,S,tree);
end
if node.hasChildNodes
    current_node = root_node.getFirstChild;
    
    while ~isempty(current_node)
        switch node.getNodeType
            case 1
                %NodeType 1 = element
                my_temp_string = current_node.getNodeName.toString.toCharArray';
                [ind,element_names] = countElement(my_temp_string,element_names);
                S  = parseElement(current_node,S,ind,my_tree);
            case [3,4,8]
                %NodeType 3 = text node
                %NodeType 4 = CDATA
                %NodeType 8 = comment node
                S = parseMiscellaneous(current_node,S);
            otherwise
        end
        current_node = current_node.getNextSibling;
    end
else
    my_temp_string = node.getTextContent.toString.toCharArrary';
    eval([tree '=''' my_temp_string ''';'])
end
end

% ----- Subfunction PARSEMISCELLANEOUS -----
function S = parseMiscellaneous(node,S,tree)
try
    my_temp_string = node.getData.toString.toCharArray';
    eval([tree '.txtd8a{end+1}=''' my_temp_string ''';'])
catch %#ok<CTCH>
    my_temp_string = node.getData.toString.toCharArray';
    eval([tree '.txtd8a{1}=''' my_temp_string ''';'])
end
end