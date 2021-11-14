#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <random>
#include <fstream>

// #include "./utils.hpp"

using namespace std;
using namespace std::chrono;

pair<int, int> mid;









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


vector<pair<int, int>> fill_vector(vector<pair<int, int>> src)
{
    vector<pair<int, int>> des;
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

    return des;
}












bool compare(pair<int, int> p1, pair<int, int> q1)
{
	pair<int, int> p = make_pair(p1.first - mid.first, p1.second - mid.second);
	pair<int, int> q = make_pair(q1.first - mid.first, q1.second - mid.second);

	int one = quad(p);
	int two = quad(q);

	if (one != two)
    {
        return (one < two);
    }
	return (p.second*q.first < q.second*p.first);
}

pair<pair<int, int>, pair<int, int>> calculate_upper_tangent(vector<pair<int, int>> a, vector<pair<int, int>> b)
{
	int n1 = a.size(), n2 = b.size();
	int ia = calc_right_m(a), ib = calc_left_m(b);

	int inda = ia, indb = ib;
	bool done = 0;
	while (!done)
	{
        // cout << "in while" << endl;
		done = 1;
		while (orientation(b[indb], a[inda], a[(inda+1)%n1]) >= 0)
		{
            // cout << "in while" << endl;
			if(orientation(b[indb], a[inda], a[(inda+1)%n1]) == 0)
			{
				if(n1 == 1)
				{
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[(inda+1)%n1], b[indb]))
					{
						a.erase(a.begin() + ((inda+1)%n1));
					} else {
                        a.erase(a.begin() + inda);
                    }
                    n1 = a.size();
                    inda = calc_right_m(a);
                    indb = calc_left_m(b);
				}
			} else
			{
				inda = (inda + 1) % n1;
			}
		}
			
		while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) <= 0)
		{
            // cout << "in while" << endl;
			if(orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) == 0)
			{
				if(n2 == 1)
				{
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[inda], b[(n2+indb-1)%n2]))
					{
						b.erase(b.begin() + ((n2+indb-1)%n2));
					} else
                    {
                        b.erase(b.begin() + indb);
                    }
                    n2 = b.size();
                    inda = calc_right_m(a);
                    indb = calc_left_m(b);
				}
			} else
			{
				indb = (n2+indb-1)%n2;
			}

			done = 0;
		}
	}
	pair<int, int> point1 = {a[inda].first, a[inda].second};
	pair<int, int> point2 = {b[indb].first, b[indb].second};
	pair<pair<int, int>, pair<int, int>> tangent = {point1, point2};
	return tangent;
}

pair<pair<int, int>, pair<int, int>> calculate_lower_tangent(vector<pair<int, int>> a, vector<pair<int, int>> b)
{
	int n1 = a.size(), n2 = b.size();
	int ia = calc_right_m(a), ib = calc_left_m(b);

	int inda = ia, indb = ib;
	bool done = 0;
	while (!done)
	{
        // cout << "in while" << endl;
		done = 1;
		while (orientation(a[inda], b[indb], b[(indb+1)%n2])>= 0)
		{
            // cout << "in while" << endl;
			if(orientation(a[inda], b[indb], b[(indb+1)%n2]) == 0)
			{
				if(n2 == 1)
				{
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[inda], b[(indb+1)%n2]))
					{
						b.erase(b.begin() + ((indb+1)%n2));
					} else
                    {
                        b.erase(b.begin() + (indb));
                    }
					n2 = b.size();
                    inda = calc_right_m(a);
                    indb = calc_left_m(b);
				}
			} else
			{
				indb=(indb+1)%n2;
			}
		}
			

		while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1]) <= 0)
		{
            // cout << "in while" << endl;
			if(orientation(b[indb], a[inda], a[(n1+inda-1)%n1]) == 0)
			{
				if(n1 == 1)
				{
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[(n1+inda-1)%n1], b[indb]))
					{
						a.erase(a.begin() + ((n1+inda-1)%n1));
					} else
                    {
                        a.erase(a.begin() + inda);
                    }
					n1 = a.size();
                    inda = calc_right_m(a);
                    indb = calc_left_m(b);
				}
			} else
			{
				inda=(n1+inda-1)%n1;
				
			}
			done=0;
		}
	}
	pair<int, int> point1 = {a[inda].first, a[inda].second};
	pair<int, int> point2 = {b[indb].first, b[indb].second};
	pair<pair<int, int>, pair<int, int>> tangent = {point1, point2};
	return tangent;
}

vector<pair<int, int>> merge_hulls(vector<pair<int, int>> a, vector<pair<int, int>> b)
{
	int n1 = a.size(), n2 = b.size();
	vector<pair<int, int>> ua, ub, la, lb;
    ua = fill_vector(a);
    la = fill_vector(a);
    ub = fill_vector(b);
    lb = fill_vector(b);
    // cout << la.size() << endl;
    pair<pair<int, int>, pair<int, int>> tu;
	pair<pair<int, int>, pair<int, int>> tl;

    #pragma omp parallel sections
    {
        #pragma omp section
        tu = calculate_upper_tangent(ua, ub);

        #pragma omp section
        tl = calculate_lower_tangent(la, lb);
    }
	
	vector<pair<int, int>> ret;
	
	pair<int, int> uppera = tu.first, upperb = tu.second;
	pair<int, int> lowera = tl.first, lowerb = tl.second;

    // #pragma omp parallel sections
    // {
    //     #pragma omp section
    //     {
            int ind = 0;
            while(a[ind].first != uppera.first && a[ind].second != uppera.second)
            {
                ind = (ind + 1) % n1;
            }
            ret.push_back(a[ind]);

            while(a[ind].first != lowera.first && a[ind].second != lowera.second)
            {
                ind = (ind + 1) % n1;
                ret.push_back(a[ind]);
            }
            ret.push_back(a[ind]);
        // }

        // #pragma omp section
        // {
            ind = 0;

            while(b[ind].first != lowerb.first && b[ind].second != lowerb.second)
            {
                ind = (ind + 1) % n2;
            }
            ret.push_back(b[ind]);

            while(b[ind].first != upperb.first && b[ind].second != upperb.second)
            {
                ind = (ind + 1) % n2;
                ret.push_back(b[ind]);
            }
            ret.push_back(b[ind]);
    //     }

    // }
	return ret;

}

vector<pair<int, int>> find_convex_hull_brute(vector<pair<int, int>> a)
{
	vector<pair<int, int>> convex_hull;

	for (int i=0; i<a.size(); i++)
	{
		for (int j=i+1; j<a.size(); j++)
		{
			int x1 = a[i].first, x2 = a[j].first;
			int y1 = a[i].second, y2 = a[j].second;

			int a1 = y1-y2;
			int b1 = x2-x1;
			int c1 = x1*y2-y1*x2;
			int pos = 0, neg = 0;
			for (int k=0; k<a.size(); k++)
			{
				if (a1*a[k].first+b1*a[k].second+c1 <= 0)
					neg++;
				if (a1*a[k].first+b1*a[k].second+c1 >= 0)
					pos++;
			}
			if (pos == a.size() || neg == a.size())
			{
				convex_hull.push_back(a[i]);
				convex_hull.push_back(a[j]);
			}
		}
	}

	mid = {0, 0};
	int n = convex_hull.size();
	for (int i=0; i<n; i++)
	{
		mid.first += convex_hull[i].first;
		mid.second += convex_hull[i].second;
		convex_hull[i].first *= n;
		convex_hull[i].second *= n;
	}
	sort(convex_hull.begin(), convex_hull.end(), compare);
	for (int i=0; i<n; i++)
    {
        convex_hull[i] = make_pair(convex_hull[i].first/n, convex_hull[i].second/n);
    }

	return convex_hull;
}

vector<pair<int, int>> find_convex_hull(vector<pair<int, int>> a)
{

	if (a.size() < 6)
    {
        return find_convex_hull_brute(a);
    }
		

	vector<pair<int, int>>left, right;
	for (int i=0; i<a.size()/2; i++)
    {
        left.push_back(a[i]);
    }
		
	for (int i=a.size()/2; i<a.size(); i++)
    {
        right.push_back(a[i]);
    }
		

    vector<pair<int, int>>left_hull;
    vector<pair<int, int>>right_hull;
    
    #pragma omp parallel sections
    {
        #pragma omp section
        left_hull = find_convex_hull(left);

        #pragma omp section
        right_hull = find_convex_hull(right);
    }
	return merge_hulls(left_hull, right_hull);
}

int main( int argc, char* argv[] )
{
    // omp_set_nested(1);
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(-100, 100);
    
    vector<int> n = {100, 1000, 10000};
    for(int i = 0; i < n.size(); i++)
    {
        vector<pair<int, int>> points;
        ofstream myfile;
        string name = "test" + to_string(i) + ".txt";
        myfile.open (name);
        for(int j = 0; j < n[i]; j++)
        {
            auto random_integer = uni(rng);
            int x = (int)random_integer;
            random_integer = uni(rng);
            int y = (int)random_integer;
            points.push_back(make_pair(x, y));
            myfile << points[j].first << ", "<< points[j].second << std::endl;
        }
        
        sort(points.begin(), points.end());
        omp_set_num_threads(4);
        auto startp = high_resolution_clock::now();
        vector<pair<int, int>> final_hull = find_convex_hull(points);
        auto endp = high_resolution_clock::now();
        myfile << "parallel: " << duration_cast<microseconds>(endp - startp).count() << endl;

        omp_set_num_threads(1);
        auto starts = high_resolution_clock::now();
        final_hull = find_convex_hull(points);
        auto ends = high_resolution_clock::now();
        myfile << "sequential: " << duration_cast<microseconds>(ends - starts).count() << endl;
        
        myfile.close();
        // cout << duration.count() << endl;
    }




    return 0;
}
