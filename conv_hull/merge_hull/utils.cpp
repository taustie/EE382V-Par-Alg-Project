#include <iostream>
#include <vector>
#include <omp.h>
#include "./utils.hpp"

using namespace std;


int quad(pair<int, int> p)
{
	if (p.first >= 0 && p.second >= 0)
		return 1;
	if (p.first <= 0 && p.second >= 0)
		return 2;
	if (p.first <= 0 && p.second <= 0)
		return 3;
	return 4;
}

int orientation(pair<int, int> a, pair<int, int> b,
				pair<int, int> c)
{
	int res = (b.second-a.second)*(c.first-b.first) -
			(c.second-b.second)*(b.first-a.first);
			

	if (res == 0)
		return 0;
	if (res > 0)
		return 1;
	return -1;
}

double dis(pair<int, int> p1, pair<int, int> p2)
{
	return (p1.first - p2.first)*(p1.first - p2.first) + (p1.second - p2.second)*(p1.second - p2.second);
}

int calc_right_m(vector<pair<int, int>> a)
{
	int n = a.size();
	int ia = 0;
	for (int i=1; i<n; i++)
		if (a[i].first > a[ia].first)
			ia = i;
	return ia;
}


int calc_left_m(vector<pair<int, int>> a)
{
	int n = a.size();
	int ia = 0;
	for (int i=1; i<n; i++)
		if (a[i].first < a[ia].first)
			ia = i;
	return ia;
}


void fill_vector(vector<pair<int, int>> src, vector<pair<int, int>> des)
{
    size_t *prefix;
    #pragma omp parallel
    {
        int ithread  = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        #pragma omp single
        {
            prefix = new size_t[nthreads+1];
            prefix[0] = 0;
        }
        vector<pair<int, int>> vec_private;
        #pragma omp for schedule(static) nowait
        for(int i=0; i<src.size(); i++) {
            vec_private.push_back(src[i]);
        }

        prefix[ithread+1] = vec_private.size();
        #pragma omp barrier
        #pragma omp single 
        {
            for(int i=1; i<(nthreads+1); i++)
            {
                prefix[i] += prefix[i-1];
            }
            des.resize(des.size() + prefix[nthreads]);
        }
        copy(vec_private.begin(), vec_private.end(), des.begin() + prefix[ithread]);
    }
    delete[] prefix;
}