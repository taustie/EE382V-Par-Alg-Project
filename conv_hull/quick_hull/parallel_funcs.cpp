#include "quick_hull.hpp"
#include "parallel_funcs.hpp"
#include "test_config.hpp"

Point * neutral_point(){
	return NULL;
}

Point * myMax(Point * a, Point * b){
	if(a == NULL){
		return b;
	}
	else if(b == NULL){
		return a;
	}

    if(a->y > b->y){
    	return a;
    }
    else if(a->y < b->y){
    	return b;
    }
    else if(a->y == b->y){
    	if(a->x > b->x){
      		return a;
       	}
       	else{
       		return b;
       	}
    }

    return NULL; // this should never happen
}

Point * myMin(Point * a, Point * b){
	if(a == NULL){
		return b;
	}
	else if(b == NULL){
		return a;
	}

    if(a->y > b->y){
    	return b;
    }
    else if(a->y < b->y){
    	return a;
    }
    else if(a->y == b->y){
    	if(a->x > b->x){
      		return b;
       	}
       	else{
       		return a;
       	}
    }

    return NULL; // this should never happen
}

void get_max_min_points_parallel(std::vector<Point *> &input_points, Point** max, Point** min){
	Point *tmp_min = neutral_point();
	Point *tmp_max = neutral_point();

    #pragma omp declare reduction \
        (maxPoint:Point *:omp_out=myMax(omp_out, omp_in)) \
        initializer(omp_priv = neutral_point())

    #pragma omp declare reduction \
        (minPoint:Point *:omp_out=myMin(omp_out, omp_in)) \
        initializer(omp_priv = neutral_point())

    #pragma omp parallel for reduction(maxPoint : tmp_max) reduction(minPoint:tmp_min)
    for(int i = 0; i < input_points.size(); i++){
		if(tmp_max == NULL){
			tmp_max = input_points.at(i);
		}
        else if(input_points.at(i)->y > tmp_max->y){
        	tmp_max = input_points.at(i);
        }
        else if(input_points.at(i)->y == tmp_max->y){
        	if(input_points.at(i)->x > tmp_max->x){
        		tmp_max = input_points.at(i);
        	}
        }

		if(tmp_min == NULL){
			tmp_min = input_points.at(i);
		}
        else if(input_points.at(i)->y < tmp_min->y){
        	tmp_min = input_points.at(i);
        }
        else if(input_points.at(i)->y == tmp_min->y){
        	if(input_points.at(i)->x < tmp_min->x){
        		tmp_min = input_points.at(i);
        	}
        }
    }

    *max = tmp_max;
    *min = tmp_min;
}

Dist_Info neutral_distance(){
	Dist_Info maxDistResult;
	maxDistResult.index = 0;
	maxDistResult.line_dist = 0.0;
	maxDistResult.p1_dist = 0.0;
	return maxDistResult;
}

Dist_Info Dist_Max_Compare(Dist_Info a, Dist_Info b){
	if(a.line_dist < b.line_dist){
		return b;
	}
	else if(a.line_dist > b.line_dist){
		return a;
	}
	else if(a.line_dist == b.line_dist){
		if(a.p1_dist < b.p1_dist){
			return a;
		}
		else{
			return b;
		}
	}

	return neutral_distance(); // should never reach here
}

void get_max_dist_parallel(std::vector<Point *> &left_points, Point** max, Point *p1, Point *p2){
	Line line = get_line(*p1, *p2);

	Dist_Info maxDistResult = neutral_distance();
	long long int index;

	#pragma omp declare reduction \
	(maxDist:Dist_Info:omp_out=Dist_Max_Compare(omp_out, omp_in)) \
	initializer(omp_priv = neutral_distance())

	#pragma omp parallel for reduction(maxDist:maxDistResult)
	for(index = 0; index < left_points.size(); index++){
		double distance = get_distance(line, *left_points.at(index));
		if(maxDistResult.line_dist < distance){
			maxDistResult.index = index;
			maxDistResult.line_dist = distance;
			maxDistResult.p1_dist = get_distance_point(*left_points.at(index), *p1);
		}
		else if(maxDistResult.line_dist == distance){
			double distance_p1 = get_distance_point(*left_points.at(index), *p1);
			if(distance_p1 < maxDistResult.p1_dist){
				maxDistResult.index = index;
				maxDistResult.p1_dist = distance_p1;
			}
		}
	}

	*max = left_points.at(maxDistResult.index);
}

void get_points_on_left_parallel(std::vector<Point *> &left_points, Point ** input_points, int size, Point* p1, Point* p2){
	Vector ref_vector = form_zero_vector(*p1, *p2);
	// Option 1: See great resource at: https://www.py4u.net/discuss/63446
	#pragma omp declare reduction (merge : std::vector<Point *> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
	#pragma omp parallel for reduction(merge: left_points)
	for(long long int i=0; i < size; i++) {
		Vector test_vector = form_zero_vector(*p2, *input_points[i]);
        if (cross_prod_orientation(ref_vector, test_vector)) {
			left_points.push_back(input_points[i]);
		}
	}

	// Option 2: use private copy and merge at end
	// #pragma omp parallel
	// {
	//     std::vector<Point *> vec_private; // fill vec_private in parallel
	// 	Vector test_vector;
	//     #pragma omp for nowait
	//     for(int i=0; i < size; i++) {
	// 		test_vector = form_zero_vector(*p2, *input_points[i]);
	// 		if (cross_prod_orientation(ref_vector, test_vector)) {
	// 			vec_private.push_back(input_points[i]);
	// 		}
	//     }
	//     #pragma omp critical
	//     left_points.insert(left_points.end(), vec_private.begin(), vec_private.end());
	// }
}


// void get_points_on_left_filter_parallel(std::vector<Point *> &left_points, Point ** input_points, int size, Point* p1, Point* p2){
// 	//prefix sum using array
// 	Vector ref_vector = form_zero_vector(*p1, *p2);
// 	int * left_points_filter = new int[size];
// 	int * left_points_result = new int[size];
//
// 	#pragma omp parallel for schedule(static)
// 	for(int k = 0; k < size; k++){
// 		Vector test_vector = form_zero_vector(*p2, *input_points[k]);
// 		if(cross_prod_orientation(ref_vector, test_vector)){
// 			left_points_filter[k] = 1;
// 		}
// 		else{
// 			left_points_filter[k] = 0;
// 		}
// 	}
//
// 	int x = 0;
// 	#pragma omp parallel for simd reduction(inscan,+: x)
// 	for (int k = 0; k < size; k++) {
// 		x += left_points_filter[k];
// 		#pragma omp scan inclusive(x)
// 		left_points_result[k] = x;
// 	}
//
// 	int new_size = left_points_result[size - 1];
// 	Point ** left_points_arr = new Point *[new_size];
//
// 	#pragma omp parallel for schedule(static)
// 	for(int i = 0; i < size; i++){
// 		if(left_points_filter[i] == 1){
// 			//std::cout << "index left_points_result.at(i): " << left_points_result.at(i) - 1 << std::endl;
// 			left_points_arr[left_points_result[i] - 1] = input_points[i];
// 		}
// 	}
//
// 	delete[] left_points_filter;
// 	delete[] left_points_result;
// 	left_points.assign(left_points_arr, left_points_arr + new_size);
// 	delete[] left_points_arr;
// }
