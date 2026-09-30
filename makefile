CXX = g++

CXXFLAGS = -std=c++20 \
           -I"$(CONDA_PREFIX)/include"

LDFLAGS = -L"$(CONDA_PREFIX)/lib" \
          -Wl,-rpath,"$(CONDA_PREFIX)/lib"

LDLIBS = -lopenblas -llapack -lblas

.PHONY: all bruteforce debug clean test

bruteforce:
	$(CXX) $(CXXFLAGS) bruteforce.cpp ./lib/*.cpp \
		$(LDFLAGS) $(LDLIBS) \
		-o brute
	@echo -e "\n"
	@./brute

test:
	$(CXX) $(CXXFLAGS) test.cpp ./lib/*.cpp \
		$(LDFLAGS) $(LDLIBS) \
		-o test
	@echo -e "\n"
	@./test

debug:
	$(CXX) $(CXXFLAGS) -g bruteforce.cpp ./lib/*.cpp \
		$(LDFLAGS) $(LDLIBS) \
		-o deb
	@echo -e "\n"
	@gdb ./deb

clean:
	rm -f main test deb
