CC=gcc 
LIBS=-lm
FLAGS= -O3 -fopenmp

Serial_Bucket: Serial_Bucket.c
	$(CC) $(FLAGS) -o Serial_Bucket Serial_Bucket.c $(LIBS)
	
Pragma_Schedule: Pragma_Schedule.c
	$(CC) $(FLAGS) -o Pragma_Schedule Pragma_Schedule.c $(LIBS)

Bin_Packing: Bin_Packing.c
	$(CC) $(FLAGS) -o Bin_Packing Bin_Packing.c $(LIBS)
	
Pragma_Task: Pragma_Task.c
	$(CC) $(FLAGS) -o Pragma_Task Pragma_Task.c $(LIBS)
	
Run_Time_System: Run_Time_System.c
	$(CC) $(FLAGS) -o Run_Time_System Run_Time_System.c $(LIBS)
	
clean:
	rm -f *.o *~ Serial_Bucket Pragma_Schedule Run_Time_System Pragma_Task Bin_Packing
	