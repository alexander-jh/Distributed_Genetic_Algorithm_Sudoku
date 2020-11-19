#
# Makefile
#	Author: Alexander Hoke
#	Date:	11/19/20
#

CC		= g++
STD		= -std=c++11
LIBS	= OptimizationProblem.hpp OptimizationProblem.cpp \
		  SudokuBoard.hpp SudokuBoard.cpp
HILL	= HillClimbing
GEN		= Genetic
TEST    = Testing/

all: HillClimber Genetic

HillClimber:
	$(CC) $(STD) -o $(HILL) $(HILL).cpp $(LIBS)

Genetic:
	$(CC) $(STD) -D GENETIC -o $(GEN) $(GEN).cpp $(LIBS) -pthread

TestHill:
	$(CC) $(STD) -o $(HILL) $(HILL).cpp $(LIBS)
	./$(HILL) $(TEST)BaseTest.txt
	./$(HILL) $(TEST)Test2.txt
	./$(HILL) $(TEST)Test3.txt
	./$(HILL) $(TEST)Test4.txt
	./$(HILL) $(TEST)Test5.txt
	./$(HILL) $(TEST)Test6.txt
	rm -rf $(HILL)

TestGenetic:
	$(CC) $(STD) -D GENETIC -o $(GEN) $(GEN).cpp $(LIBS) -pthread
	./$(GEN) $(TEST)BaseTest.txt
	./$(GEN) $(TEST)Test2.txt
	./$(GEN) $(TEST)Test3.txt
	./$(GEN) $(TEST)Test4.txt
	./$(GEN) $(TEST)Test5.txt
	./$(GEN) $(TEST)Test6.txt
	rm -rf $(GEN)

clean:
	rm -rf $(GEN) $(HILL)