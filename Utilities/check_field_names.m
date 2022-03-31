function check_field_names(struct1, fieldnames_cellarray)
% Throws an error if there are any field names that are not in the struct.
% If there are invalid fields, the error message will give the indices of
% the invalid fields

valid_fieldnames = isfield(struct1, fieldnames_cellarray);
invalid_fields = find(valid_fieldnames == 0);

if ~isempty(invalid_fields)
    error('The following field names are not valid: %d', invalid_fields);
end

end

