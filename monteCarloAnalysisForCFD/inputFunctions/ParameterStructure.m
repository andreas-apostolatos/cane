% Class definition that makes sure that the problem structure objects get
% passed by reference

classdef ParameterStructure < handle
  % MEMBER VARIABLES ----------------------------------------------------
  properties
    fldMsh
    homDBC
    inhomDBC
    valuesInhomDBC
    parameters
  end
  
  % MEMBER FUNCTIONS ------------------------------------------------------
  methods
      % Constructor
    function this = ParameterStructure(varargin)
        % this = ParameterStructure(varargin)
        % The constructor can be called without arguments (all fields 
        % remain empty) or with a problemStructure, in which case all
        % fields of the problemStructure that are defined in
        % ParameterStructure will be copied.
        if numel(varargin) == 1
            % Copy all relevant fields of the argument
            fieldNames = fields(varargin{1});
            for k=1:length(fieldNames)
                if isprop(this, fieldNames{k})
                    this.(fieldNames{k}) = varargin{1}.(fieldNames{k});
                end
            end
        end
    end
    
    % Set fields
    function [] = setField(this, fieldNames, fieldValues)
        % [] = setField(fieldNames, fieldValues)
        % Set the value(s) of specific field(s). Substructures up to any 
        % depth can be set as well. Since multiple fields can be set at
        % once, the arguments must be cells or cell arrays.
        % Example1: ParameterStructure.setField({'fldMsh.nodes'},{nodes})
        % Example2: .setField({'homDBC','inhomDBC'},{homDBC,inhomDBC})
        for k=1:numel(fieldNames)
            field = textscan(fieldNames{k},'%s','Delimiter','.');
            if ~isempty(field{1})
                this = setfield(this, field{1}{:}, fieldValues{k});
            else
                this.(fieldNames{k}) = fieldValues{k};
            end
        end
    end
    
  end
end