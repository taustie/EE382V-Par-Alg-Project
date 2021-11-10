linux: test_qh.cpp ../quick_hull/parallel_funcs.cpp
	g++ -O3 -std=c++17 -fopenmp test_qh.cpp \
	../common/verify_conv_hull.cpp  ../quick_hull/quick_hull.cpp ../quick_hull/parallel_funcs.cpp \
	-o test_qh -I. -I../quick_hull/ -I../common/

mac: test_qh.cpp ../quick_hull/parallel_funcs.cpp
	g++ -O3 -std=c++17 -Xclang -fopenmp test_qh.cpp \
	../common/verify_conv_hull.cpp ../quick_hull/quick_hull.cpp ../quick_hull/parallel_funcs.cpp \
	-o test_qh -I. -I../quick_hull/ -I../common/ -lomp

clean:
	rm test_qh
