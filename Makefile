#OBJS specifies which files to compile as part of the project
#TARGET = fileio

CC = g++

COMPILER_FLAGS = -Wall
LINKER_FLAGS = -lgsl -lgslcblas -lm
DFALGS = -DHAVE_INLINE -DGSL_C99_INLINE

INC_PATH = -I/usr/local/include
LIB_PATH = -L/usr/local/lib

#This is the target that compiles our executable
all : $(TARGET).cpp
	$(CC) $(COMPILER_FLAGS) $(DFALGS) $(INC_PATH) $(LIB_PATH) -o$(TARGET) $(TARGET).cpp $(LINKER_FLAGS)


clean :
	rm $(TARGET)

