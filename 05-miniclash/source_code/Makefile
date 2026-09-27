CXX ?= g++
CXXFLAGS := -O3 -march=native -mtune=native -flto -fomit-frame-pointer -fno-semantic-interposition -DNDEBUG -std=gnu++17
COLLISION_OBJS := block0.o block1.o block1stevens00.o block1stevens01.o block1stevens10.o block1stevens11.o block1wang.o main.o md5.o timer.o

.PHONY: clean all

all: md5fastcoll run

md5fastcoll: $(COLLISION_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

run: runner.o
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(COLLISION_OBJS) runner.o md5fastcoll run

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
