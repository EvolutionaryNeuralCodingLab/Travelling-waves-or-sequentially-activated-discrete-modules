function row_num = get_row_num_by_dt_dx(results_table,dt,dx,dt_col_name,dx_col_name)
%GET_ROW_NUM_BY_DT_DX returns the row number corresponding to dt, dx in simulation results table.
%   If dt_col_name,dx_col_name are optional with defaults "deltaT",
%   "distance"

if nargin==3
    dt_col_name = "deltaT";
    dx_col_name = "distance";
elseif nargin~=5
    error("Number of inputs should either be 3 (results_table, dt, dx), or 5 (results_table,dt,dx,dt_col_name,dx_col_name)")
end

row_num = find(results_table.(dt_col_name) == dt & results_table.(dx_col_name)==dx);

end

