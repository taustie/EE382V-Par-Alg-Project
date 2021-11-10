// #include <iostream>
// #include <cmath>
// #include <vector>
// #include <omp.h>
// #include<bits/stdc++.h>
// #include <chrono>
// #include <random>

// using namespace std::chrono;
// using namespace std;






// pair<int, int> mid;
// // determines the quadrant of a point
// // (used in compare())
// int quad(pair<int, int> p)
// {
//     if (p.first >= 0 && p.second >= 0)
//         return 1;
//     if (p.first <= 0 && p.second >= 0)
//         return 2;
//     if (p.first <= 0 && p.second <= 0)
//         return 3;
//     return 4;
// }
// // Checks whether the line is crossing the polygon
// int orientation(pair<int, int> a, pair<int, int> b,
//                 pair<int, int> c)
// {
//     int res = (b.second-a.second)*(c.first-b.first) -
//               (c.second-b.second)*(b.first-a.first);
//     if (res == 0)
//         return 0;
//     if (res > 0)
//         return 1;
//     return -1;
// }
// // compare function for sorting
// bool compare(pair<int, int> p1, pair<int, int> q1)
// {
//     pair<int, int> p = make_pair(p1.first - mid.first,
//                                  p1.second - mid.second);
//     pair<int, int> q = make_pair(q1.first - mid.first,
//                                  q1.second - mid.second);
//     int one = quad(p);
//     int two = quad(q);
//     if (one != two)
//         return (one < two);
//     return (p.second*q.first < q.second*p.first);
// }
// // Finds upper tangent of two polygons 'a' and 'b'
// // represented as two vectors.
// vector<pair<int, int>> merger(vector<pair<int, int> > a,
//                               vector<pair<int, int> > b)
// {
//     // n1 -> number of points in polygon a
//     // n2 -> number of points in polygon b
//     cout << "merge hull" << endl;
//     int n1 = a.size(), n2 = b.size();
//     int ia = 0, ib = 0;
//     for (int i=1; i<n1; i++)
//         if (a[i].first > a[ia].first)
//             ia = i;
//     // ib -> leftmost point of b
//     for (int i=1; i<n2; i++)
//         if (b[i].first < b[ib].first)
//             ib=i;
//     // finding the upper tangent
//     int inda = ia, indb = ib;
//     bool done = 0;
//     while (!done)
//     {
//         done = 1;
//         while (orientation(b[indb], a[inda], a[(inda+1)%n1]) >=0)
//             inda = (inda + 1) % n1;
//         while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) <=0)
//         {
//             indb = (n2+indb-1)%n2;
//             done = 0;
//         }
//     }
//     int uppera = inda, upperb = indb;
//     inda = ia, indb=ib;
//     done = 0;
//     int g = 0;
//     while (!done)//finding the lower tangent
//     {
//         done = 1;
//         while (orientation(a[inda], b[indb], b[(indb+1)%n2])>=0)
//             indb=(indb+1)%n2;
//         while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1])<=0)
//         {
//             inda=(n1+inda-1)%n1;
//             done=0;
//         }
//     }
//     int lowera = inda, lowerb = indb;
//     vector<pair<int, int>> ret;
//     //ret contains the convex hull after merging the two convex hulls
//     //with the points sorted in anti-clockwise order
//     int ind = uppera;
//     ret.push_back(a[uppera]);
//     while (ind != lowera)
//     {
//         ind = (ind+1)%n1;
//         ret.push_back(a[ind]);
//     }
//     ind = lowerb;
//     ret.push_back(b[lowerb]);
//     while (ind != upperb)
//     {
//         ind = (ind+1)%n2;
//         ret.push_back(b[ind]);
//     }
//     return ret;
// }
// // Brute force algorithm to find convex hull for a set
// // of less than 6 points
// vector<pair<int, int>> bruteHull(vector<pair<int, int>> a)
// {
//     // Take any pair of points from the set and check
//     // whether it is the edge of the convex hull or not.
//     // if all the remaining points are on the same side
//     // of the line then the line is the edge of convex
//     // hull otherwise not
//     set<pair<int, int> >s;
//     for (int i=0; i<a.size(); i++)
//     {
//         for (int j=i+1; j<a.size(); j++)
//         {
//             int x1 = a[i].first, x2 = a[j].first;
//             int y1 = a[i].second, y2 = a[j].second;
//             int a1 = y1-y2;
//             int b1 = x2-x1;
//             int c1 = x1*y2-y1*x2;
//             int pos = 0, neg = 0;
//             for (int k=0; k<a.size(); k++)
//             {
//                 if (a1*a[k].first+b1*a[k].second+c1 <= 0)
//                     neg++;
//                 if (a1*a[k].first+b1*a[k].second+c1 >= 0)
//                     pos++;
//             }
//             if (pos == a.size() || neg == a.size())
//             {
//                 s.insert(a[i]);
//                 s.insert(a[j]);
//             }
//         }
//     }
//     vector<pair<int, int>>ret;
//     for (auto e:s)
//         ret.push_back(e);
//     // Sorting the points in the anti-clockwise order
//     mid = {0, 0};
//     int n = ret.size();
//     for (int i=0; i<n; i++)
//     {
//         mid.first += ret[i].first;
//         mid.second += ret[i].second;
//         ret[i].first *= n;
//         ret[i].second *= n;
//     }
//     sort(ret.begin(), ret.end(), compare);
//     for (int i=0; i<n; i++)
//         ret[i] = make_pair(ret[i].first/n, ret[i].second/n);
//     return ret;
// }

// vector<pair<int, int>> divide_seq(vector<pair<int, int>> a, int threads)
// {
//     // If the number of points is less than 6 then the
//     // function uses the brute algorithm to find the
//     // convex hull

//     if (a.size() <= 5)
//         return bruteHull(a);
//     // left contains the left half points
//     // right contains the right half points
//     vector<pair<int, int>>left, right;
//     for (int i=0; i<a.size()/2; i++)
//         left.push_back(a[i]);
//     for (int i=a.size()/2; i<a.size(); i++)
//         right.push_back(a[i]);
//     // convex hull for the left and right sets
// 	vector<pair<int, int>>left_hull;
// 	vector<pair<int, int>>right_hull;
// 	left_hull = divide_seq(left, threads/2);
// 	right_hull = divide_seq(right, threads - threads/2);
//     return merger(left_hull, right_hull);
// }

// // Returns the convex hull for the given set of points
// vector<pair<int, int>> divide_par(vector<pair<int, int>> a, int threads)
// {
//     if(threads == 1) {
//         cout << "going to seq" << endl;
//         return divide_seq(a, threads);
//     }
//     if (a.size() <= 5)
//         return bruteHull(a);
//     // left contains the left half points
//     // right contains the right half points
//     vector<pair<int, int>>left, right;
//     for (int i=0; i<a.size()/2; i++)
//         left.push_back(a[i]);
//     for (int i=a.size()/2; i<a.size(); i++)
//         right.push_back(a[i]);
//     // convex hull for the left and right sets
// 	vector<pair<int, int>>left_hull;
// 	vector<pair<int, int>>right_hull;
// 	#pragma omp parallel sections
// 	{
// 		#pragma omp section
//         {
//             cout << omp_get_num_threads() << endl;
// 			left_hull = divide_par(left, threads/2);
//         }
            
// 		#pragma omp section
// 			right_hull = divide_par(right, threads - threads/2);
// 	}
//     return merger(left_hull, right_hull);
// }








// int main( int argc, char* argv[] )
// {
// 	omp_set_num_threads(4);
//     omp_set_nested(1);

// 	std::random_device rd;     // only used once to initialise (seed) engine
// 	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
// 	std::uniform_int_distribution<int> uni(-100,100);
//     vector<pair<int, int> > a;
// 	for(int i = 0; i < 20; i++) {
// 		 // guaranteed unbiased
// 		auto random_integer = uni(rng);
// 		int x = (int)random_integer;
// 		random_integer = uni(rng);
// 		int y = (int)random_integer;
// 		// cout << x << endl;
// 		// cout << y << endl;
// 		a.push_back(make_pair(x, y));
// 	}
//     int n = a.size();
//     // sorting the set of points according
//     // to the x-coordinate
//     sort(a.begin(), a.end());
	
	
// 	vector<pair<int, int> >ans;
// 	auto start = high_resolution_clock::now();
// 	ans = divide_par(a, 4);
// 	auto end = high_resolution_clock::now();
// 	auto duration = duration_cast<microseconds>(end - start);
// 	cout << duration.count() << endl;
// 	// cout << "convex hull:\n";
// 	// for (auto e:ans)
// 	// 	cout << e.first << " "
// 	// 		<< e.second << endl;
//     return 0;
// }




// A divide and conquer program to find convex
// hull of a given set of points.
#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include<bits/stdc++.h>
#include <chrono>
#include <random>
using namespace std::chrono;
using namespace std;

pair<int, int> mid;

int find_quadrant(int x, int y)
{
    if (x >= 0 && y >= 0)
    {
        return 1;
    }
    else if (x <= 0 && y >= 0)
    {
        return 2;
    }
    else if (x <= 0 && y <= 0)
    {
        return 3;
    }
    else 
    {
        return 4;
    }
}

int orientation(pair<int, int> a, pair<int, int> b, pair<int, int> c)
{
    int res = (b.second-a.second)*(c.first-b.first) - (c.second-b.second)*(b.first-a.first);
    
    if (res == 0)
        return 0;
    if (res > 0)
        return 1;
    return -1;
}

bool compare_points(pair<int, int> point1, pair<int, int> point2)
{
    int q1 = find_quadrant(point1.first - mid.first, point1.second - mid.second);
    int q2 = find_quadrant(point2.first - mid.first, point2.second - mid.second);
    if (q1 != q2) 
    {
        return (q1 < q2);
    }
        
    return (point1.second * point2.first < point1.second * point2.first);
}

vector<pair<int, int>> merge_hull(vector<pair<int, int>> left, vector<pair<int, int>> right)
{
    cout << "merge hull" << endl;

    int l_size = left.size();
    int r_size = right.size();

    int right_m = 0;
    int left_m = 0;

    for (int i = 1; i < l_size; i++)
    {
        if (left[i].first > left[right_m].first) 
        {
            right_m = i;
        }
    }

    for (int i = 1; i < r_size; i++)
    {
        if (right[i].first < right[right_m].first) 
        {
            left_m = i;
        }
    }

    int inda = right_m;
    int indb = left_m;
    bool done = 0;
    while (!done)
    {
        done = 1;
        while (orientation(right[indb], left[inda], left[(inda+1)%l_size]) >=0)
        {
            inda=(inda+1)%l_size;
        }
            
        while (orientation(left[inda], right[indb], right[(r_size+indb-1)%r_size]) <=0)
        {
            indb = (r_size+indb-1)%r_size;
            done = 0;
        }
    }

    int uppera = inda, upperb = indb;
    inda = right_m;
    indb=left_m;
    done = 0;
    int g = 0;
    while (!done)
    {
        done = 1;
        while (orientation(left[inda], right[indb], right[(indb+1)%r_size])>=0)
        {
            indb=(indb+1)%r_size;
        }
        while (orientation(right[indb], left[inda], left[(l_size+inda-1)%l_size])<=0)
        {
            inda=(l_size+inda-1)%l_size;
            done=0;
        }
    }
    int lowera = inda, lowerb = indb;

    vector<pair<int, int>> ret;
    //ret contains the convex hull after merging the two convex hulls
    //with the points sorted in anti-clockwise order
    int ind = uppera;
    ret.push_back(left[uppera]);
    while (ind != lowera)
    {
        ind = (ind+1)%l_size;
        ret.push_back(left[ind]);
    }
    ind = lowerb;
    ret.push_back(right[lowerb]);
    while (ind != upperb)
    {
        ind = (ind+1)%r_size;
        ret.push_back(right[ind]);
    }
    return ret;
}

vector<pair<int, int>> brute_force(vector<pair<int, int>> points, int start, int end)
{
    vector<pair<int, int>> selected_points;

    for(int i = start; i <= end; i++)
    {
        selected_points.push_back(points[i]);
    }

    set<pair<int, int>> convex_hull;

    for (int i=0; i < selected_points.size(); i++)
    {
        for (int j=i+1; j < selected_points.size(); j++)
        {
            int px1 = selected_points[i].first;
            int py1 = selected_points[i].second;
            int px2 = selected_points[j].first;
            int py2 = selected_points[j].second;
            int a = py1 - py2;
            int b = px2 - px1;
            int c = (px1 * py2) - (py1 * px2);

            int pos_points = 0;
            int neg_points = 0;
            for (int k = 0; k < selected_points.size(); k++)
            {
                int side = (a*selected_points[k].first) + (b*selected_points[k].second) + c;
                if (side <= 0)
                {
                    neg_points++;
                }
                    
                if (side >= 0)
                {
                    pos_points++;
                }
            }
            if (pos_points == selected_points.size() || neg_points == selected_points.size())
            {
                convex_hull.insert(selected_points[i]);
                convex_hull.insert(selected_points[j]);
            }
        }
    }

    vector<pair<int, int>> sorted_hull;
    for (auto e: convex_hull)
    {
        sorted_hull.push_back(e);
    }
        
    mid = {0, 0};
    for (int i = 0; i < sorted_hull.size(); i++)
    {
        mid.first += sorted_hull[i].first;
        mid.second += sorted_hull[i].second;
        sorted_hull[i].first *= sorted_hull.size();
        sorted_hull[i].second *= sorted_hull.size();
    }
    sort(sorted_hull.begin(), sorted_hull.end(), compare_points);
    for (int i = 0; i < sorted_hull.size(); i++)
    {
        sorted_hull[i] = make_pair(sorted_hull[i].first/sorted_hull.size(), sorted_hull[i].second/sorted_hull.size());
    }
        
    return sorted_hull;
}

vector<pair<int, int>> find_convex_hull_seq(vector<pair<int, int>> points, int start, int end)
{
    if (end - start + 1 <= 5)
    {
        return brute_force(points, start, end);
    }

    int mid_ind = points.size() / 2;

    vector<pair<int, int>>right_hull = find_convex_hull_seq(points, mid_ind, end);
    vector<pair<int, int>>left_hull = find_convex_hull_seq(points, start, mid_ind-1);   

    vector<pair<int, int>> ans =  merge_hull(left_hull, right_hull);
    return ans;
}

vector<pair<int, int>> find_convex_hull_par(vector<pair<int, int>> points, int start, int end, int threads)
{

    if(threads == 1)
    {
        cout << "going to seq" << endl;
        return find_convex_hull_seq(points, start, end);
    }

    if (end - start + 1 <= 5)
    {
        return brute_force(points, start, end);
    }

    int mid_ind = points.size() / 2;

    vector<pair<int, int>>right_hull;
    vector<pair<int, int>>left_hull;
    #pragma omp parallel sections
	{
		#pragma omp section
            right_hull = find_convex_hull_par(points, mid_ind, end, threads/2);

        #pragma omp section
            left_hull = find_convex_hull_par(points, start, mid_ind-1, threads - (threads/2));    
    }

    vector<pair<int, int>> ans =  merge_hull(left_hull, right_hull);
	for (auto e:ans)
		cout << "hello" << endl;
    return ans;
}



// Driver code
int main( int argc, char* argv[] )
{
	omp_set_num_threads(4);
    omp_set_nested(1);
    
    
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(-100,100);
    vector<pair<int, int>> points;
	for(int i = 0; i < 15; i++)
    {
		auto random_integer = uni(rng);
		int x = (int)random_integer;
		random_integer = uni(rng);
		int y = (int)random_integer;
		points.push_back(make_pair(x, y));
	}
    
    sort(points.begin(), points.end());
    auto start = high_resolution_clock::now();
    cout << "hi" << endl;
    vector<pair<int, int>> final_hull = find_convex_hull_par(points, 0, points.size() - 1, 4);
    auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << duration.count() << endl;
    return 0;
}
