from utils import solve_base_solution, solve_first_type_equity_solution, solve_second_type_equity_solution
import gurobipy as gp
from gurobipy import GRB



# for i in range(3, 12):
#     files_name = "instance" + ("0" + str(i) if i<10 else str(i))
solve_second_type_equity_solution("instance07")

