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
    """Method used to solve the instance of an ETP problem contained in 3 files with name {files_name}

    Args:
        files_name (str): name of the files containing the instance to solve
        equity_measure (int): integer variable containing the type of equity measure (objective function) to apply
                                        0: base equity measure
                                        1: equity measure based on average distance between conflicting exams
                                        2: equity measure based on number of students having back to back conflicting exams
                                        3: equity measure based on minimum distance between any 2 conflicting exams
                                        4: equity measure based on average distance between conflicting exams using MAD as penalty
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
    
    #force base problem in case of additional restrictions selected
    if additional_restriction > 0: equity_measure=0

    ##populate vars from files
    slo, exm, stu, n = read_instances(files_name)
    
    T = list(range(1, slo+1))
    E = list(exm.keys())
    S = list(stu.keys())

    #prepare variables for gurobi
    x_dict = dict()
    y_dict = dict()

    #initialize variable x[e,t]
    for e in E:
        for t in T:
            x_dict[e, t] = 0
    
    #linearization: new binary decision var y[e1,e2,t1, t2] = 1 if x[e1, t1]=1 and x[e2,t2]=1
    for e1 in E:
        #commented to make it work for equity measure =2 
        #for e2 in range(e1+1, len(E)+1):
        for e2 in E:
            for t1 in T:
                for t2 in T:
                    #linearize x[e1,t1] * x[e2,t2] using y[e1,t1,e2,t2]
                    if n[e1,e2]>0:
                        y_dict[e1,t1,e2,t2] = 0



    # Create Gurobi model
    model = gp.Model("Exam_Timetabling")

    # Create gurobi binary decision variable x[e, t]
    x = model.addVars(x_dict, vtype=GRB.BINARY, name="x")

    #linearization
    y = model.addVars(y_dict, vtype=GRB.BINARY, name="y")


    #compute total number of conflicts between exams
    num_conflicting_exm = np.sum(1 for e1 in E for e2 in range(e1+1, len(E)+1) if n[e1, e2] > 0)


##############################################################################################
####################################    OBJ FUNCTION   #######################################
##############################################################################################

    match equity_measure:
        case 0:
            if additional_restriction != 3:
                model.setObjective(
                    gp.quicksum((2**(5 - abs(t1-t2)) )* n[e1, e2] *  y[e1,t1,e2,t2] / len(S) 
                                for e1 in E for e2 in range(e1+1,len(E)+1) for t1 in T for t2 in T if abs(t1-t2)<=5 and n[e1,e2]>0),
                    GRB.MINIMIZE
                )
            else:
                #introduce bonus variable b[t] =1 if 6 timeslots with no conflicting exams scheduled
                b_dict = dict()
                for t in range(1,len(T)-4):
                    b_dict[t]=0
                b = model.addVars(b_dict, vtype=GRB.BINARY, name="b")

                model.setObjective(
                    gp.quicksum((2**(5 - abs(t1-t2)) )* n[e1, e2] *  y[e1,t1,e2,t2] / len(S) 
                                for e1 in E for e2 in range(e1+1,len(E)+1) for t1 in T for t2 in T if abs(t1-t2)<=5 and n[e1,e2]>0) 
                                #subtract all the bonus (1 every 6 timeslots withoutconflicting exams, even if overlapping)
                                - gp.quicksum(b[t] for t in range(1, len(T)-4)),
                    GRB.MINIMIZE
                )
        case 1:
            #maximize average distance between any two conflicting exams
            model.setObjective(
                gp.quicksum((abs(t1-t2) * y[e1,t1,e2,t2]) / num_conflicting_exm
                            for e1 in E for e2 in range(e1+1, len(E)+1) for t1 in T for t2 in T if n[e1, e2] > 0),
                GRB.MAXIMIZE
            )
        case 2:
            #initialize variable z[s]
            z_dict = dict()
            for s in S:
                z_dict[s] = 0

            #binary variable z[s] = 1 if student s has at least 2 conflicting exams scheduled back to back
            z = model.addVars(z_dict, vtype=GRB.BINARY, name="z")   
            
            #minimize number of students with back to back conflicting exams 
            model.setObjective(
                gp.quicksum(z[s] for s in S),
                GRB.MINIMIZE
            )
        case 3:
            w = model.addVar(0, vtype=GRB.INTEGER, name="w")
            #Maximize the minimum distance between any 2 conflicting exams
            model.setObjective(
                w,
                GRB.MAXIMIZE
            )
        case 4:
            #d[e1,e2] integer variable containing absolute difference of the distance of e1-e2 
            #from the average distance of any two conflicting exams
            d_dict = dict()
            for e1 in E:
                for e2 in range(e1+1, len(E)+1):
                    if n[e1,e2]>0:
                        d_dict[e1,e2] = 0

            d = model.addVars(d_dict, vtype=GRB.INTEGER, name="d")

            #we want to minimize d => distance of two conflicting exams should be 
            #as close as possible to average distance between any 2 conflicting exams 
            # model.setObjective(gp.quicksum(d[e1, e2] / num_conflicting_exm for e1 in E for e2 in range(e1+1, len(E)+1) if n[e1,e2]>0)
            #     ,GRB.MINIMIZE
            # )

            #objective function taking into consideration both average value and MAD trying to maximize average while minimizing MAD 
            #with certain weights assigned to both measures
            model.setObjective(
                 gp.quicksum((abs(t1-t2) * y[e1,t1,e2,t2]) / num_conflicting_exm
                            for e1 in E for e2 in range(e1+1, len(E)+1) for t1 in T for t2 in T if n[e1, e2] > 0) - 
                            1.5 *(gp.quicksum(d[e1, e2] / num_conflicting_exm for e1 in E for e2 in range(e1+1, len(E)+1) if n[e1,e2]>0)),
                GRB.MAXIMIZE
            )

##############################################################################################
#####################################   CONSTRAINTS   ########################################
##############################################################################################

    #each exam can be scheduled only once in examination period
    for e in E:
        model.addConstr(gp.quicksum(x[e, t] for t in T) == 1, name=f"{e}_scheduled_once")

    #can't schedule conflicting exams in same time slot (valid for every case except for additional restriction 4)
    match additional_restriction:
        case 0|1|2|3:
            for e1 in E:
                for e2 in range(e1+1, len(E)+1):
                    if n[e1, e2] > 0:
                        for t in T:
                            model.addConstr(x[e1,t] + x[e2,t] <= 1)
        case 4:
            for t in T:
                #at most 3 conflicting pairs of exams in same time slot
                model.addConstr(gp.quicksum(y[e1,t,e2,t] for e1 in E for e2 in range(e1+1, len(E)+1) if n[e1, e2]>0) <= 3, name=f"max_3_conflicting_exams_{t}")


    #link y to x (linearize x[e1,t1] * x[e2,t2])
    for e1 in E:
        for e2 in E:
            if n[e1, e2]>0:
                for t1 in T:
                    for t2 in T:
                        model.addConstr(y[e1,t1,e2,t2]>=x[e1,t1] + x[e2,t2] - 1 ,
                                        name=f"y_{e1}_{t1}_{e2}_{t2}_gt")
                        model.addConstr(y[e1,t1,e2,t2]<=x[e1,t1] ,
                                        name=f"y_{e1}_{t1}_{e2}_{t2}_lt_x_{e1}_{t1}")
                        model.addConstr(y[e1,t1,e2,t2]<=x[e2,t2] ,
                                        name=f"y_{e1}_{t1}_{e2}_{t2}_lt_x_{e2}_{t2}")
    
    if equity_measure==2:
        for s in S:
            for t in range(1, len(T)):
                for e1 in stu[s]:
                    for e2 in stu[s]:
                        #exams in which the student is enrolled are conflicting by definition
                        if e1!=e2:
                        #if n[e1,e2]>0:
                            model.addConstr(
                                #if at least two conflicting exams scheduled back back for student s, then z[s] is forced to get value 1
                                #otherwise it will be pushed to 0 by minimization function where possible
                                z[s] >=  y[e1,t,e2,t+1]
                                , name=f"back_to_back_{s}_{e1}_{e2}_{t}"
                                )
    if equity_measure==3:
        for e1 in E:
            for e2 in range(e1+1, len(E)+1):
                if n[e1,e2]>0:
                        #w is minimum distance between any two conflicting exams
                        model.addConstr(w <=  gp.quicksum(abs(t1-t2) * y[e1,t1,e2,t2] for t1 in T for t2 in T))

    if equity_measure==4:
        for e1 in E:
            for e2 in range(e1+1, len(E)+1):
                if n[e1,e2]>0:
                        #use linearization of absolute value
                        # d = |f(x)| => d >= x, d >= -x and minimize d
                        model.addConstr(d[e1,e2] >= gp.quicksum(abs(t1-t2) * y[e1,t1,e2,t2] for t1 in T for t2 in T) - 
                                        gp.quicksum(abs(t3-t4) * y[e3,t3,e4,t4] / num_conflicting_exm for e3 in E for e4 in range (e3+1, len(E)+1) for t3 in T for t4 in T if n[e3,e4]>0))
                        
                        model.addConstr(d[e1,e2] >= gp.quicksum(abs(t3-t4) * y[e3,t3,e4,t4] / num_conflicting_exm for e3 in E for e4 in range (e3+1, len(E)+1) for t3 in T for t4 in T if n[e3,e4]>0) - 
                                        gp.quicksum(abs(t1-t2) * y[e1,t1,e2,t2] for t1 in T for t2 in T))
        
        #insert desired minimum average distance between exams 
        # model.addConstr(gp.quicksum(abs(t1-t2) * y[e1,t1,e2,t2] for e1 in E for e2 in range (e1+1, len(E)+1) for t1 in T for t2 in T if n[e1,e2]>0)/num_conflicting_exm >= 10)
    
    
    
    match additional_restriction:
        #case 1: at most 3 consecutive time slots with conflicting exams
        #case 2: if two conflicting exams in consecutive time slots, 
                    #then no conflicting exams can be scheduled in the next 3 time slots (every 4 time slots, only two of them can have conflicting exams)
        #case 3: bonus introduced when 6 consecutive timeslots do not have conflicting exams (allows overlapping of slots)
        case 1|2|3:
            k_dict = dict()
            for t1 in range(1, len(T)):
                for t2 in range(t1+1, len(T)+1):
                    #k[t1,t2]=1 if t1 and t2 have conflicting exams scheduled
                    k_dict[t1, t2] = 0
        
            k = model.addVars(k_dict, vtype=GRB.BINARY, name="k")
        
            for t1 in range(1, len(T)):
                for t2 in range(t1+1, len(T)+1):
                    #if we have zero conflicting exams in t1, t2, then k is forced to 0
                    model.addConstr(k[t1,t2] <= gp.quicksum(y[e1,t1,e2,t2] for e1 in E for e2 in E if n[e1,e2]>0)
                                    , name=f"{t}_{t+1}_conflicting_exams_gt")
                    
                    #if we have m conflicting exams in t1,t2, we force k[t1,t2] to go to 1 dividing by the total number of 
                    #conflicting exams
                    model.addConstr(k[t1,t2] >= gp.quicksum(y[e1,t1,e2,t2] for e1 in E for e2 in E if n[e1,e2]>0)/num_conflicting_exm
                                    , name=f"{t}_{t+1}_conflicting_exams_lt")

        
            
                
            if additional_restriction==1:
                for t in range(1, len(T)-2):
                    #can't have more than 3 consecutive time slots with conflicting exams
                    # if t=1 and have conflicts in 1 and 2, 2 and 3, then I will force time slot 4 to not have conflicts with neither timeslot 3 nor 5
                    model.addConstr(k[t+2, t+3] <= 2 - (k[t,t+1] + k[t+1, t+2]))
                    if t <= len(T)-4:
                        model.addConstr(k[t+3, t+4] <= 2 - (k[t,t+1] + k[t+1, t+2]))
            
            elif additional_restriction==2:
                for t in range(1, len(T)-1):
                    #only two consecutive time slots can have conflicting exams, then 3 timeslots without conflicting exams (between consecutive time slots)
                    #if t=1 and have conflicts in 1 and 2, then I will have no consecutive conflicts in 2 and 3, 3 and 4, 4 and 5, 5 and 6 
                    #to ensure that 3-4-5 are without conflicts
                    model.addConstr(k[t+1, t+2] <= 1 - k[t, t+1])
                    if t <= len(T)-3:
                        model.addConstr(k[t+2, t+3] <= 1 - k[t, t+1])
                    if t <= len(T)-4:
                        model.addConstr(k[t+3, t+4] <= 1 - k[t, t+1])
                    if t <= len(T)-5:
                        model.addConstr(k[t+4, t+5] <= 1 - k[t, t+1])
            
            elif additional_restriction==3:
                for t in range(1, len(T)-4):                    
                        #if the sum is 0, bonus is 1
                        #check between all possible combination between exams from t to t+6
                        model.addConstr(b[t]>= 1 - (gp.quicksum(k[t3,t4] for t3 in range(t,t+5) for t4 in range(t3+1, t+6))))
                        
                        # #if the sum is greater than 0, force the bonus to 0
                        model.addConstr(b[t]<= 1 - (gp.quicksum(k[t3,t4] for t3 in range(t,t+5) for t4 in range(t3+1, t+6))/num_conflicting_exm))

        
##############################################################################################
#################################   SHOW FOUND SOLUTION   ####################################
##############################################################################################
 
    model.setParam(GRB.param.TimeLimit, 600)
    # model.setParam('Presolve', 2)
    #model.setParam('MIPFocus', 3)

    #Optimize the model
    model.optimize()



    #Display results
    if model.SolCount > 0:
        for e, t in x.keys():
                if x[e, t].x > 0.5:
                    print(f"Exam {e} is scheduled in time-slot {t}")

        #implemented only for base constraints
        if check_feasibility(x, n, T, additional_restriction):
            
            #save found solution
            if equity_measure == 0:
                model.write(f"solutions/equity_measure_{equity_measure}/{files_name}_additional_restriction_{additional_restriction}.sol")
            else:
                model.write(f"solutions/equity_measure_{equity_measure}/{files_name}.sol")

            if model.status == GRB.OPTIMAL:
                print("Found optimal solution")
            
        else:
            print("Solution not feasible")
    else:
        print("No solution found.")

    model.dispose()








def check_feasibility(solution, n : dict, T : list, additional_restriction : int) -> bool:
    """Method used to check the feasibility of a solution

    Args:
        solution {gurobi solution}: Gurobi solution to be verified
        n {dict}: Dictionary containing number of students enrolled in pairs of conflicting exams
    Returns:
        feasibile {bool} : True if the solution respects the constraints
    """
    dict_schedules = dict()
    dict_time_slot_schedule = dict()
    for t in T:
        dict_time_slot_schedule[t] = list()

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
            #build dictionary containing all exams scheduled in a certain time slot
            dict_time_slot_schedule[t].append(e)

    #2nd constraint: can't have conflicting exam in same time slot
    match additional_restriction:
        case 1|2|3:
            for e1 in dict_schedules.keys():
                for e2 in dict_schedules.keys():
                    if n[e1, e2] > 0:
                        if dict_schedules[e1] == dict_schedules[e2]:
                            return False
        case 4:
            for t in T:
                #print (dict_time_slot_schedule[t])
                count = np.sum(1 for e1 in dict_time_slot_schedule[t] for e2 in dict_time_slot_schedule[t] if n[e1,e2]>0 and e2>e1)
                #print(f"timeslot {t} count {count}")
                if count >3:
                    return False
    
            
    return True






def solve_instance_diff_formulation(files_name : str):
    """Method used to solve the instance of an ETP problem contained in 3 files with name {files_name}

    Args:
        files_name (str): name of the files containing the instance to solve
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
    y_dict = dict()

    #initialize variable x[e,t]
    for e in E:
        for t in T:
            x_dict[e, t] = 0
    
    #linearization: new binary decision var y[e1,e2,t1, t2] = 1 if x[e1, t1]=1 and x[e2,t2]=1
    for e1 in E:
        for e2 in E:
            if n[e1,e2]>0:
                for i in range (1, 6):
                        if i < len(T):
                            for t in range(1, len(T)-i+1):
                                y_dict[e1, e2, i] = 0



    # Create Gurobi model
    model = gp.Model("Exam_Timetabling")

    # Create gurobi binary decision variable x[e, t]
    x = model.addVars(x_dict, vtype=GRB.BINARY, name="x")

    #linearization
    y = model.addVars(y_dict, vtype=GRB.BINARY, name="y")


##############################################################################################
####################################    OBJ FUNCTION   #######################################
##############################################################################################

    model.setObjective(
        gp.quicksum((2**(5 - i) )* n[e1, e2] *  y[e1,e2,i] / len(S) 
                    for e1 in E for e2 in E for i in range(1, 6) if n[e1,e2]>0),
        GRB.MINIMIZE
    )

##############################################################################################
#####################################   CONSTRAINTS   ########################################
##############################################################################################

    #each exam can be scheduled only once in examination period
    for e in E:
        model.addConstr(gp.quicksum(x[e, t] for t in T) == 1, name=f"{e}_scheduled_once")

    for e1 in E:
        for e2 in range(e1+1, len(E)+1):
            if n[e1, e2] > 0:
                for t in T:
                    model.addConstr(x[e1,t] + x[e2,t] <= 1)
       
    #link y to x (linearize x[e1,t1] * x[e2,t2])
    for e1 in E:
        for e2 in E:
            if n[e1,e2]>0:
                for i in range (1, 6):
                        if i < len(T):
                            for t in range(1, len(T)-i+1):
                                model.addConstr(y[e1,e2,i] >= x[e1, t] + x[e2, t+i] -1 )
    
 
##############################################################################################
#################################   SHOW FOUND SOLUTION   ####################################
##############################################################################################
 
    model.setParam(GRB.param.TimeLimit, 1200)
    #model.setParam('Presolve', 2)
    #model.setParam('MIPFocus', 3)

    #Optimize the model
    model.optimize()

    #Display results
    if model.SolCount > 0:
        for e, t in x.keys():
                if x[e, t].x > 0.5:
                    print(f"Exam {e} is scheduled in time-slot {t}")

        #implemented only for base constraints
        if check_feasibility(x, n, T, 0):

            #save found solution
            model.write(f"solutions/new_formulation/{files_name}.sol")
           
            if model.status == GRB.OPTIMAL:
                print("Found optimal solution")
            
        else:
            print("Solution not feasible")
    else:
        print("No solution found.")

    model.dispose()