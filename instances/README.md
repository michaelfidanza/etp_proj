# Benchmark instances for the DODMproject2023

This folder contains a set of banchmark instances of the problem faced with this project.


## Format of the input files

Each instance, say instanceXX with XX a natural number, is defined by 3 plain text (ASCII) files, with the same name and different extensions:

instanceXX.slo: defines the length of the examination period for instance XX.
Format: a single value in the format
INT
where INT is the number of available time-slots.
This integer number is always positive and is called tmax; the time-slots will have numeric IDs  1, 2, ... , tmax.

instanceXX.stu: defines the exams in which each student is enrolled.
Format: One line for each enrollment. Each line has the format
INT1 INT2
where INT1 is the student ID and INT2 is the ID of the exam in which student INT1 is enrolled.
The student ID is a string that starts with an 's' character and continues with the decimal encoding of the number identifying the student.
The exam ID is the decimal encoding of the number identifying the exam, but it is left-padded with 0's in order to standardize its length to 4 characters.

instanceXX.exm: is an auxiliary file that summarizes the total number of students enrolled per exam. It can be used to check the data integrity of the instance files.
Format: One line per exam. Each line has the format
INT1 INT2
where INT1 is the exam ID and INT2 is the number of students enrolled in exam INT1.


## Format of the output file
An output file is meant to contain either the string "UNFEASIBLE" or a fesasible solution of smallest possible penality.
In case it offers a feasible solution, then it should contain as many lines as the number of exams. Each line has the format
INT1 INT2
where INT1 is the ID of an exam and INT2 is the ID of the time-slot where that exam is allocated.
