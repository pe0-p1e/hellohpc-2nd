CXX := g++
OPT := -O3 -std=gnu++17 -march=native -mtune=native -flto -fipa-pta -funroll-loops -DNDEBUG -fno-semantic-interposition -fno-plt
CORE_OBJS := block0.o block1.o block1stevens00.o block1stevens01.o block1stevens10.o block1stevens11.o block1wang.o md5.o
CLI_OBJS := main.o timer.o
LIBS := -pthread

.PHONY := clean all

all: md5fastcoll run

md5fastcoll: $(CORE_OBJS) $(CLI_OBJS)
	$(CXX) $(OPT) $^ $(LIBS) -o $@

run: runner.o $(CORE_OBJS)
	$(CXX) $(OPT) $^ $(LIBS) -o $@

clean:
	rm -f $(CORE_OBJS) $(CLI_OBJS) runner.o md5fastcoll run

%.o: %.cpp
	$(CXX) $(OPT) -pthread -c $< -o $@
