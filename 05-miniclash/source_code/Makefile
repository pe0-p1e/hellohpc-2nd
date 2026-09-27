CXX := g++
CPPFLAGS := -O3 -std=gnu++17 -march=native -flto -fwhole-program -fipa-pta -funroll-loops -DNDEBUG -fno-semantic-interposition -fno-plt
OBJS := block0.o block1.o block1stevens00.o block1stevens01.o block1stevens10.o block1stevens11.o block1wang.o main.o md5.o timer.o
LIBS :=

.PHONY := clean all

all: md5fastcoll run

md5fastcoll: $(OBJS)
	$(CXX) $(CPPFLAGS) $(LIBS) $^ -o $@

run: run.sh
	cp -f run.sh run
	chmod +x run

clean:
	rm -f $(OBJS) md5fastcoll run

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@
