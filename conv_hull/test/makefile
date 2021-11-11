linux: test_qh.cpp ../quick_hull/parallel_funcs.cpp test_mh.cpp ../merge_hull/merge_hull_par.cpp
	g++ -O3 -std=c++17 -fopenmp test_qh.cpp \
	../common/verify_conv_hull.cpp  ../quick_hull/quick_hull.cpp ../quick_hull/parallel_funcs.cpp \
	-o test_qh -I. -I../quick_hull/ -I../common/

	g++ -O3 -std=c++17 -fopenmp test_mh.cpp \
	../merge_hull/merge_hull_par.cpp \
	-o test_mh -I. -I../merge_hull/ -I../common/ -I../quick_hull/

mac: test_qh.cpp ../quick_hull/parallel_funcs.cpp test_mh.cpp ../merge_hull/merge_hull_par.cpp
	g++ -O3 -std=c++17 -Xclang -fopenmp test_qh.cpp \
	../common/verify_conv_hull.cpp ../quick_hull/quick_hull.cpp ../quick_hull/parallel_funcs.cpp \
	-o test_qh -I. -I../quick_hull/ -I../common/ -lomp

	g++ -O3 -std=c++17 -Xclang -fopenmp test_mh.cpp \
	../merge_hull/merge_hull_par.cpp \
	-o test_mh -I. -I../merge_hull/ -I../common/ -I../quick_hull/ -lomp

	#../common/verify_conv_hull.cpp ../quick_hull/quick_hull.cpp ../quick_hull/parallel_funcs.cpp

clean:
	rm test_qh
	rm test_mh
