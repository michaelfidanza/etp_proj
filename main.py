from utils_linearized import solve_instance, solve_instance_diff_formulation




for i in range(1, 10):
    files_name = "instance" + ("0" + str(i) if i<10 else str(i))
    if i!= 6:
        solve_instance(files_name, equity_measure=0, additional_restriction=0)

