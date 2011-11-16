function theStruct = my_parseXML(filename)
% PARSEXML converts an XML file to a MATLAB structure. Note that for
% whatever reason the method .getChildNodes() does not work in MATLAB.
% input: filename = the name of the XML file to be parsed.
% output: theStruct = the contents of the XML file are converted into a
% struct.

try
   xdoc = xmlread(filename);
catch err
   disp(err.message)
   error('Failed to read XML file %s.',filename);
end


%Initialize theStruct
current_node = xdoc.getDocumentElement; %Locates the node at the root of the XML tag tree
my_location = ['theStruct.' current_node.getNodeName]; %The location in theStruct where data is currently being created
.hasChildNodes
.hasAttributes


.getFirstChild
.getNextSibling
.getTextContent %to access text node from parent
.getData %to access the contents of a text node directly
try
   theStruct = parseNode(xdoc);
catch err
   disp(err.message)
   error('Unable to parse XML file %s.',filename);
end


% ----- Subfunction PARSECHILDNODES -----
function children = parseNode(my_node)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Subfunction MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Subfunction PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end