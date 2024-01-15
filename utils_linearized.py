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




def solve_instance(files_name : str, equity_measure : int = 0, additional_restriction : int = 0):
    """_summary_

    Args:
        files_name (str): name of the files containing the instance to solve
        equity_measure (int): integer variable containing the type of equity measure (objective function) to apply
                                        0: base equity measure
                                        1: equity measure based on average distance between conflicting exams
                                        2: equity measure based on number of students having back to back conflicting exams
        additional_restriction (int): integer variable containing the type of restriction to apply
                                        0: base solution
                                        1: at most 3 consecutive time slots with conflicting exams
                                        2: if two consecutive slots have conflicting exams, then no conflicts for 3 slots
                                        3: bonus profit each time no conflicting exams are scheduler for 6 slots
                                        4: at most 3 conflicting pairs of exams can be scheduled in same time slot
    """

##############################################################################################
####################################   INITIALIZATION   ######################################
##############################################################################################
    

    ##populate vars from files
    slo, exm, stu, n = read_instances(files_name)
    
    T = list(range(1, slo+1))
    E = list(exm.keys())
    S = list(stu.keys())

    #prepare variables for gurobi
    x_dict = dict()
    z_dict = dict()

    #initialize variable x[e,t]
    for e in E:
        for t in T:
            x_dict[e, t] = 0
    
    #linearization: new binary decision var z[e1,e2,t1, t2] = 1 if x[e1, t1]=1 and x[e2,t2]=1
    for e1 in E:
        for e2 in range(e1+1, len(E)+1):
            for t1 in T:
                for t2 in T:
                    if abs(t1-t2) <= 5 and abs(t1-t2)>=1 and n[e1,e2]>0:
                        z_dict[e1,t1,e2,t2] = 0

    # Create Gurobi model
    model = gp.Model("Exam_Timetabling")

    # Create gurobi binary decision variable x[e, t]
    x = model.addVars(x_dict, vtype=GRB.BINARY, name="x")

    #linearization
    z = model.addVars(z_dict, vtype=GRB.BINARY, name="z")

##############################################################################################
####################################    OBJ FUNCTION   #######################################
##############################################################################################

    match equity_measure:
        case 0:
            model.setObjective(
                gp.quicksum((2**(5 - abs(t1-t2)) )* n[e1, e2] *  z[e1,t1,e2,t2] / len(S) 
                            for e1 in E for e2 in range(e1+1,len(E)+1) for t1 in T for t2 in T if abs(t1-t2)<=5 and abs(t1-t2)>=1 and n[e1,e2]>0),
                GRB.MINIMIZE
            )
        case 1:
            return
        case 2:
            return

##############################################################################################
#####################################   CONSTRAINTS   ########################################
##############################################################################################

    #each exam can be scheduled only once in examination period
    for e in E:
        model.addConstr(gp.quicksum(x[e, t] for t in T) == 1, name=f"{e}_scheduled_once")

    #can't schedule conflicting exams in same time slot
    match additional_restriction:
        case 0,1,2,3:
            for e1 in E:
                for e2 in range(e1+1, len(E)+1):
                    if n[e1, e2] > 0:
                        for t in T:
                            model.addConstr(x[e1,t] + x[e2,t] <= 1)
        case 4:
            return

    #link z to x (linearize x[e1,t1] * x[e2,t2])
    for e1 in E:
        for e2 in range(e1+1, len(E)+1):
            for t1 in T:
                for t2 in T:
                    if abs(t1-t2)<=5 and abs(t1-t2)>=1 and n[e1, e2] > 0:
                        model.addConstr(z[e1,t1,e2,t2]>=x[e1,t1] + x[e2,t2] - 1 ,
                                        name=f"z_{e1}_{t1}_{e2}_{t2}_gt")
                        model.addConstr(z[e1,t1,e2,t2]<=x[e1,t1] ,
                                        name=f"z_{e1}_{t1}_{e2}_{t2}_lt_x_{e1}_{t1}")
                        model.addConstr(z[e1,t1,e2,t2]<=x[e2,t2] ,
                                        name=f"z_{e1}_{t1}_{e2}_{t2}_lt_x_{e2}_{t2}")



##############################################################################################
#####################################   GUROBI PARS   ########################################
##############################################################################################

    #enable solution pool, model won't stop at first optimal solution
    #model.setParam(GRB.Param.PoolSearchMode, 2)  
    model.setParam(GRB.param.TimeLimit, 1800)
    #model.setParam(GRB.param.Presolve, 2)
    # model.setParam(GRB.Param.Method, 0)  # 0 = primal simplex
    # model.setParam(GRB.Param.Cuts, 0)    # 0 = disable cuts
    # model.setParam(GRB.Param.MIPFocus, 3)  # 3 = focus on finding feasible solutions
    # model.setParam(GRB.Param.Heuristics, 0)  # 0 = disable heuristics

    #Optimize the model
    model.optimize()



    #Display results
    if model.SolCount > 0:
        for e, t in x.keys():
                if x[e, t].x > 0.5:
                    print(f"Exam {e} is scheduled in time-slot {t}")

        if check_feasibility(x, n):
            #save solution found
            model.write(f"solutions/{files_name}_base_restriction_{additional_restriction}.sol")
            if model.status == GRB.OPTIMAL:
                print("Found optimal solution")
            
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
            if n[e1, e2] > 0:
                if dict_schedules[e1] == dict_schedules[e2]:
                    return False

    
            
    return True


