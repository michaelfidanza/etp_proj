import numpy as np
import gurobipy as gp
from gurobipy import GRB


def read_instances(files_name : str) -> tuple:
    """Method used to initialize all the variables for the selected instance, starting from files .slo, .exm and .stu

    Args:
        files_name {string}: Files name for reading instances file (.slo, stu, .exm)

    Returns:
        slo {integer} : number of time slots available in the examination period
        exm {dict} : dictionary containing number of students enrolled in each exam
        stu {dict} : dictionary containing list of exams for each student
        n {dict} : dictionary that, for each pair of exams, contains number of enrolled students (used to find conflicting exams)
    """

    ##initialize file names
    slo_file_path = "instances/" + files_name + ".slo"
    exm_file_path = "instances/" + files_name + ".exm"
    stu_file_path = "instances/" + files_name + ".stu"

    exm = dict()
    stu = dict()

    #read number of time slots
    with open(slo_file_path, 'r') as f:
        slo = int(f.readline())

    #read number of enrolled students for each exam
    with open(exm_file_path, 'r') as f:
        for line in f.readlines():
            if line != '\n':
                exm[int(line.split(" ")[0])] = int(line.split(" ")[1])

    #read student enrollments in each exam
    with open(stu_file_path, 'r') as f:
        for line in f.readlines():
            if line != '\n':
                str_key = int(line.split(" ")[0].replace("s",""))
                str_value = int(line.split(" ")[1])
                if str_key in stu.keys():
                    stu[str_key].append(str_value)
                else:
                    stu[str_key] = [str_value]
        
    n = dict()
    #initialize n with all zeros
    for e1 in exm.keys():
        #for e2 in range(e1+1,len(exm.keys())+1):
        for e2 in exm.keys():    
            n[e1,e2] = 0

    #populate n with number of conflicts for each pair of exams
    #start from students to find conflicting exams
    for s in stu.keys():
        #if student is enrolled in more than 1 exam
        if len(stu[s]) > 1:
            for i in range(len(stu[s])):
                for j in range(len(stu[s])):
                    if i != j:
                        #increase counter for all conflicting exams for that student
                        n[stu[s][i], stu[s][j]] += 1
    return(slo, exm, stu, n)




def solve_base_solution(files_name : str):
    ##populate vars from files
    slo, exm, stu, n = read_instances(files_name)


    T = list(range(1, slo+1))
    E = list(exm.keys())
    S = list(stu.keys())


    #prepare variables for gurobi
    x_dict = dict()

    #initialize variable x[e,t]
    for e in E:
        for t in T:
            x_dict[e, t] = 0

    # Create Gurobi model
    model = gp.Model("Exam_Timetabling")

    # Create gurobi binary decision variable x[e, t]
    x = model.addVars(x_dict, vtype=GRB.BINARY, name="x")

    #set objective function
    model.setObjective(
        gp.quicksum((2**(5 - abs(t1 - t2)) )* n[e1, e2] * x[e1, t1] * x[e2,t2] / len(S) 
                    for e1 in E for e2 in range(e1+1, len(E)+1) for t1 in T for t2 in T if abs(t1-t2)<=5 and abs(t1-t2)>=1),
        GRB.MINIMIZE
    )

    #Add constraints

    #each exam can be scheduled only once in examination period
    for e in E:
        model.addConstr(gp.quicksum(x[e, t] for t in T) == 1, name=f"{e}_scheduled_once")

    #can't have conflicting exams scheduled in same time slot
    for e1 in E:
        for e2 in range(e1+1, len(E)+1):
            if n[e1, e2] > 0:
                for t in T:
                    model.addConstr(x[e1, t] * x[e2, t] == 0, name=f"no_conflicting_exams_{e1}_{e2}_{t}")

    #enable solution pool, model won't stop at first optimal solution
    model.setParam(GRB.Param.PoolSearchMode, 2)  
    model.setParam(GRB.param.TimeLimit, 1200)
    model.setParam(GRB.param.Presolve, 2)

    #Optimize the model
    model.optimize()

    #Display results
    if model.SolCount > 0:
        if check_feasibility(x, n):
            #save solution found
            model.write(f"solutions/{files_name}_base.sol")
            if model.status == GRB.OPTIMAL:
                print("Found optimal solution")
            for e, t in x.keys():
                if x[e, t].x > 0.5:
                    print(f"Exam {e} is scheduled in time-slot {t}")
        else:
            print("Solution not feasible")
    else:
        print("No solution found.")

    model.dispose()




def solve_first_type_equity_solution(files_name : str):
    ##populate vars from files
    slo, exm, stu, n = read_instances(files_name)


    T = list(range(1, slo+1))
    E = list(exm.keys())
    S = list(stu.keys())

    #compute total number of conflicting exmams
    num_conflicting_exm = np.sum(1 for e1 in E for e2 in range(e1+1, len(E)+1) if n[e1, e2] > 0)

    #prepare variables for gurobi
    x_dict = dict()

    #initialize variable x[e,t]
    for e in E:
        for t in T:
            x_dict[e, t] = 0

    # Create Gurobi model
    model = gp.Model("Exam_Timetabling")

    # Create gurobi binary decision variable x[e, t]
    x = model.addVars(x_dict, vtype=GRB.BINARY, name="x")

    #set objective function
    model.setObjective(
        gp.quicksum((abs(t1-t2) * x[e1, t1] * x[e2, t2]) / num_conflicting_exm
                    for e1 in E for e2 in range(e1+1, len(E)+1) for t1 in T for t2 in T if n[e1, e2] > 0),
        GRB.MAXIMIZE
    )

    #Add constraints

    #each exam can be scheduled only once in examination period
    for e in E:
        model.addConstr(gp.quicksum(x[e, t] for t in T) == 1, name=f"{e}_scheduled_once")

    #can't have conflicting exams scheduled in same time slot
    for e1 in E:
        for e2 in range(e1+1, len(E)+1):
            if n[e1, e2] > 0:
                for t in T:
                    model.addConstr(x[e1, t] * x[e2, t] == 0, name=f"no_conflicting_exams_{e1}_{e2}_{t}")

    #enable solution pool, model won't stop at first optimal solution
    model.setParam(GRB.Param.PoolSearchMode, 2)  
    model.setParam(GRB.param.TimeLimit, 600)
    # model.setParam(GRB.param.Presolve, 2)

    #Optimize the model
    model.optimize()



    #Display results
    if model.SolCount > 0:
        if check_feasibility(x, n):
            #save solution found
            model.write(f"solutions/{files_name}_first_type_equity.sol")
            if model.status == GRB.OPTIMAL:
                print("Found optimal solution")
            for e, t in x.keys():
                if x[e, t].x > 0.5:
                    print(f"Exam {e} is scheduled in time-slot {t}")
        else:
            print("Solution not feasible")
    else:
        print("No solution found.")
    model.dispose()




def solve_second_type_equity_solution(files_name : str):
    ##populate vars from files
    slo, exm, stu, n = read_instances(files_name)


    T = list(range(1, slo+1))
    E = list(exm.keys())
    S = list(stu.keys())


    #prepare variables for gurobi
    x_dict = dict()
    z_dict = dict()
    y_dict = dict()

    #initialize variable x[e,t]
    for e in E:
        for t in T:
            x_dict[e, t] = 0

    #initialize variable z[s, t]
    for s in S:
        for t in T:
            z_dict[s, t] = 0

    #initialize variable y[s]
    for s in S:
        y_dict[s] = 0

    # Create Gurobi model
    model = gp.Model("Exam_Timetabling")

    #Create gurobi binary decision variable x[e, t]
    x = model.addVars(x_dict, vtype=GRB.BINARY, name="x")

    #binary variable z[s, t] = 1 if conflicting exams for student s scheduled in t and t+1
    z = model.addVars(z_dict, vtype=GRB.BINARY, name="z")

    #binary variable y[s] = 1 if student s has at least 2 conflicting exams scheduled back to back
    y = model.addVars(y_dict, vtype=GRB.BINARY, name="y")

    #minimize number of students with back to back conflicting exams 
    model.setObjective(
        gp.quicksum(y[s] for s in S),
        GRB.MINIMIZE
    )

    #Add constraints

    #each exam can be scheduled only once in examination period
    for e in E:
        model.addConstr(gp.quicksum(x[e, t] for t in T) == 1, name=f"{e}_scheduled_once")

    #can't have conflicting exams scheduled in same time slot
    for e1 in E:
        for e2 in range(e1+1, len(E)+1):
            if n[e1, e2] > 0:
                for t in T:
                    model.addConstr(x[e1, t] * x[e2, t] == 0, name=f"no_conflicting_exams_{e1}_{e2}_{t}")

    #link z to x
    for s in S:
        for t in range(1, len(T)):
            model.addConstr(
                #z[s, t] can have value 1 at most, otherwise it would against previous constraing(no conflicts in same time slot)
                z[s, t] == np.sum(x[e1, t] * x[e2, t+1] for e1 in stu[s] for e2 in stu[s] if n[e1, e2] > 0) 
                , name=f"back_to_back_{s}_{e1}_{e2}_{t}"
                )
            
    #link y to z
    for s in S:
        model.addConstr(
            #y[s] is forced by minimization objective function to 0 when the sum is 0, and it is pushed to 1 when when the sum is > 0
            #because it is divided by a number which is greater than the maximum value obtainable by the sum
            y[s] >= (gp.quicksum(z[s, t] for t in T) / float(len(stu[s])))
            , name=f"{s}_at_least_one_back_to_back"
        )


    #enable solution pool, model won't stop at first optimal solution
    model.setParam(GRB.Param.PoolSearchMode, 2)  
    model.setParam(GRB.param.TimeLimit, 600)
    # model.setParam(GRB.param.Presolve, 2)

    #Optimize the model
    model.optimize()




    #Display results
    if model.SolCount > 0:
        if check_feasibility(x, n):
            #save solution found
            model.write(f"solutions/{files_name}_second_type.sol")
            if model.status == GRB.OPTIMAL:
                print("Found optimal solution")
            for e, t in x.keys():
                if x[e, t].x > 0.5:
                    print(f"Exam {e} is scheduled in time-slot {t}")
        else:
            print("Solution not feasible")
    else:
        print("No solution found.")
    model.dispose()








def check_feasibility(solution, n : dict) -> bool:
    """Method used to check the feasibility of a solution

    Args:
        solution {gurobi solution}: Gurobi solution to be verified
        n {dict}: Dictionary containing number of students enrolled in pairs of conflicting exams
    Returns:
        feasibile {bool} : True if the solution respects the constraints
    """
    dict_schedules = dict()
    
    #1st constraint: all exams must be scheduled exactly once
    i = 1
    scheduled = 0
    for e, t in solution.keys():
        if i != e:
            #when we start checking a new exam, the sum of scheduling for prev exam must be 1
            if scheduled != 1:
                return False
            
            #check next exam
            i += 1
            scheduled = 0
        #if exam is scheduled, increase counter
        if solution[e, t].x > 0.5:
            scheduled += 1
            #build dictionary containing scheduling solution
            dict_schedules[e] = t

    #2nd constraint: can't have conflicting exam in same time slot
    for e1 in dict_schedules.keys():
        for e2 in dict_schedules.keys():
            if (e1, e2) in n.keys() and n[e1, e2] > 0:
                if dict_schedules[e1] == dict_schedules[e2]:
                    return False

    
            
    return True


