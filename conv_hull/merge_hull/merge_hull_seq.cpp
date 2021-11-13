#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include<bits/stdc++.h>
#include <chrono>
#include <random>
using namespace std::chrono;
using namespace std;

// stores the centre of polygon (It is made
// global because it is used in compare function)
pair<int, int> mid;

// determines the quadrant of a point
// (used in compare())
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

// Checks whether the line is crossing the polygon
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

// compare function for sorting
bool compare(pair<int, int> p1, pair<int, int> q1)
{
	pair<int, int> p = make_pair(p1.first - mid.first,
								p1.second - mid.second);
	pair<int, int> q = make_pair(q1.first - mid.first,
								q1.second - mid.second);

	int one = quad(p);
	int two = quad(q);

	if (one != two)
		return (one < two);
	return (p.second*q.first < q.second*p.first);
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

double dis(pair<int, int> p1, pair<int, int> p2)
{
	return (p1.first - p2.first)*(p1.first - p2.first)
	 + (p1.second - p2.second)*(p1.second - p2.second);
}

pair<pair<int, int>, pair<int, int>> calculate_upper_tangent(vector<pair<int, int>> a, vector<pair<int, int>> b)
{
	int n1 = a.size(), n2 = b.size();
	int ia = calc_right_m(a), ib = calc_left_m(b);

	int inda = ia, indb = ib;
	bool done = 0;
	while (!done)
	{
		cout << "break while 1" << endl;
		done = 1;
		while (orientation(b[indb], a[inda], a[(inda+1)%n1]) >= 0)
		{
			if(orientation(b[indb], a[inda], a[(inda+1)%n1]) == 0)
			{
				if(n1 == 1)
				{
					cout << "break while 1" << endl;
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[(inda+1)%n1], b[indb]))
					{
						
						// a.erase(a.begin() + ((inda+1)%n1));
						break;
					}
					inda = (inda + 1) % n1;
					// else
					// {
					// 	a.erase(a.begin() + inda);
						
					// }
					
					// inda = calc_right_m(a);
					// indb = calc_left_m(b);
					// n1 = a.size();
				}
			} else
			{
				inda = (inda + 1) % n1;
			}
			// cout << "while 1" << endl;
			// cout << inda << endl;
			// cout << b[indb].first;
			// cout << a[inda].first << endl;
		}
			
		while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) <= 0)
		{

			if(orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) == 0)
			{
				if(n2 == 1)
				{
					cout << "break while 2" << endl;
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[inda], b[(n2+indb-1)%n2]))
					{
						
						// b.erase(b.begin() + ((n2+indb-1)%n2));
						break;
						
					}
					// else
					// {
					// 	b.erase(b.begin() + indb);
						
					// }
					
					// indb = calc_left_m(b);
					// inda = calc_right_m(a);
					// n2 = b.size();
					indb = (n2+indb-1)%n2;
					done = 0;
					// break;
				}
			} else
			{
				indb = (n2+indb-1)%n2;
			}

			done = 0;
			// cout << "while 2" << endl;
			// cout << b[indb].first;
			// cout << a[inda].first << endl;
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
	cout << "in the merger" << endl;
	int ia = calc_right_m(a), ib = calc_left_m(b);

	// finding the upper tangent
	int inda = ia, indb = ib;
	bool done = 0;
	while (!done)//finding the lower tangent
	{
		cout << "outer while 1" << endl;
		done = 1;
		while (orientation(a[inda], b[indb], b[(indb+1)%n2])>= 0)
		{

			if(orientation(a[inda], b[indb], b[(indb+1)%n2]) == 0)
			{
				if(n2 == 1)
				{
					cout << "break while 3" << endl;
					
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[inda], b[(indb+1)%n2]))
					{
						// b.erase(b.begin() + ((indb+1)%n2));
						break;
					}
					//  else
					// {
					// 	// b.erase(b.begin() + indb);
					// }
					indb=(indb+1)%n2;
					// indb = calc_left_m(b);
					// inda = calc_right_m(a);
					// n2 = b.size();
				}
			} else
			{
				indb=(indb+1)%n2;
			}

			
			// cout << "while 3" << endl;
			// cout << b[indb].first;
			// cout << a[inda].first << endl;
		}
			

		while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1]) <= 0)
		{
			cout << "outer while 2" << endl;
			if(orientation(b[indb], a[inda], a[(n1+inda-1)%n1]) == 0)
			{
				if(n1 == 1)
				{
					cout << "break while 4" << endl;
					
					break;
				} else
				{
					if(dis(a[inda], b[indb]) > dis(a[(n1+inda-1)%n1], b[indb]))
					{
						// a.erase(a.begin() + ((n1+inda-1)%n1));
						break;
					}
					// else
					// {
					// 	a.erase(a.begin() + inda);
					// }
					inda=(n1+inda-1)%n1;
					// inda = calc_right_m(a);
					// indb = calc_left_m(b);
					// n1 = a.size();
				}
			} else
			{
				inda=(n1+inda-1)%n1;
				
			}
			done=0;
			// cout << "while 4" << endl;
			// cout << b[indb].first;
			// cout << a[inda].first << endl;
		}
	}
	pair<int, int> point1 = {a[inda].first, a[inda].second};
	pair<int, int> point2 = {b[indb].first, b[indb].second};
	pair<pair<int, int>, pair<int, int>> tangent = {point1, point2};
	return tangent;
}

vector<pair<int, int>> merger(vector<pair<int, int> > a, vector<pair<int, int> > b)
{
	int n1 = a.size(), n2 = b.size();
	vector<pair<int, int>> ua, ub, la, lb;
	for(int i = 0; i < a.size(); i++)
	{
		la.push_back(a[i]);
		ua.push_back(a[i]);
	}

	for(int i = 0; i < b.size(); i++)
	{
		lb.push_back(b[i]);
		ub.push_back(b[i]);
	}
	pair<pair<int, int>, pair<int, int>> tu = calculate_upper_tangent(ua, ub);
	pair<pair<int, int>, pair<int, int>> tl = calculate_lower_tangent(la, lb);
	
	vector<pair<int, int>> ret;
	
	pair<int, int> uppera = tu.first, upperb = tu.second;
	pair<int, int> lowera = tl.first, lowerb = tl.second;
	int ind = 0;
	// ret.push_back(uppera);
	// ret.push_back(lowera);
	// ret.push_back(lowerb);
	// ret.push_back(upperb);
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
	return ret;

}

// Brute force algorithm to find convex hull for a set
// of less than 6 points
vector<pair<int, int>> bruteHull(vector<pair<int, int>> a)
{
	// Take any pair of points from the set and check
	// whether it is the edge of the convex hull or not.
	// if all the remaining points are on the same side
	// of the line then the line is the edge of convex
	// hull otherwise not
	set<pair<int, int> >s;

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
				s.insert(a[i]);
				s.insert(a[j]);
			}
		}
	}

	vector<pair<int, int>>ret;
	for (auto e:s)
		ret.push_back(e);

	// Sorting the points in the anti-clockwise order
	mid = {0, 0};
	int n = ret.size();
	for (int i=0; i<n; i++)
	{
		mid.first += ret[i].first;
		mid.second += ret[i].second;
		ret[i].first *= n;
		ret[i].second *= n;
	}
	sort(ret.begin(), ret.end(), compare);
	for (int i=0; i<n; i++)
		ret[i] = make_pair(ret[i].first/n, ret[i].second/n);

	return ret;
}

// Returns the convex hull for the given set of points
vector<pair<int, int>> divide(vector<pair<int, int>> a)
{
	// If the number of points is less than 6 then the
	// function uses the brute algorithm to find the
	// convex hull
	if (a.size() < 6)
		return bruteHull(a);

	// left contains the left half points
	// right contains the right half points
	vector<pair<int, int>>left, right;
	for (int i=0; i<a.size()/2; i++)
		left.push_back(a[i]);
	for (int i=a.size()/2; i<a.size(); i++)
		right.push_back(a[i]);

	// convex hull for the left and right sets
	vector<pair<int, int>>left_hull = divide(left);
	vector<pair<int, int>>right_hull = divide(right);

	// merging the convex hulls
	return merger(left_hull, right_hull);
}




void write_input(vector<pair<int, int>> input_points){
	std::ofstream input_points_file;
    input_points_file.open("../visualize/input.txt");
	for(int i = 0; i < input_points.size(); i++){
		input_points_file << "(" << input_points[i].first << ", "<< input_points[i].second << ")" << std::endl;
	}
    input_points_file.close();
}

void write_output(vector<pair<int, int>> output_points){
	// output_points.sort(compare_angles);
	
	std::ofstream output_points_file;
    output_points_file.open("../visualize/output.txt");
	for (int i = 0; i < output_points.size(); i++){
		output_points_file << "(" << output_points[i].first << ", "<< output_points[i].second << ")" << std::endl;
	}
	output_points_file << "(" << output_points[0].first << ", "<< output_points[0].second << ")" << std::endl;
    output_points_file.close();
}




// Driver code
int main( int argc, char* argv[] )
{ 
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(-100, 100);
    vector<pair<int, int>> points;
	for(int i = 0; i < 500; i++)
    {
		auto random_integer = uni(rng);
		int x = (int)random_integer;
		random_integer = uni(rng);
		int y = (int)random_integer;
		points.push_back(make_pair(x, y));
	}
    //parallel sort

	write_input(points);
    sort(points.begin(), points.end());
    auto start = high_resolution_clock::now();
    vector<pair<int, int>> final_hull = divide(points);
	write_output(final_hull);
    auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << duration.count() << endl;
    return 0;
}
